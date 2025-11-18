#!/usr/bin/env python3
import os
# Set environment variable for headless rendering
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from scipy.sparse import load_npz
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
import logging
from typing import Dict, List, Optional, Tuple
import argparse
import warnings

def setup_logging(output_dir: str) -> None:
    """Configure logging with file and console handlers"""
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    file_handler = logging.FileHandler(f"{output_dir}/pangenome_viz.log")
    file_handler.setFormatter(formatter)
    
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    root_logger.addHandler(file_handler)
    
    # Add console handler only if not already present
    if not any(isinstance(h, logging.StreamHandler) for h in root_logger.handlers):
        console_handler = logging.StreamHandler()
        console_handler.setFormatter(formatter)
        root_logger.addHandler(console_handler)

def analyze_pangenome(matrix: np.ndarray, matrix_type: str = "") -> Dict[str, int]:
    """Analyze pangenome distribution with improved threshold handling"""
    n_rows, n_cols = matrix.shape
    presence_counts = np.sum(matrix, axis=1)
    
    stats = {
        'total': n_rows,
        'core': np.sum(presence_counts >= 0.90 * n_cols),
        'strict_core': np.sum(presence_counts == n_cols),
        'soft_core': np.sum((presence_counts >= 0.85 * n_cols) & (presence_counts < 0.90 * n_cols)),
        'shell': np.sum((presence_counts >= 0.25 * n_cols) & (presence_counts < 0.85 * n_cols)),
        'cloud': np.sum(presence_counts < 0.25 * n_cols),
        'singleton': np.sum(presence_counts == 1),
        'shared': np.sum((presence_counts > 1) & (presence_counts < n_cols))
    }
    
    logging.info(f"\n{matrix_type} Matrix Analysis:")
    logging.info(f"Dimensions: {stats['total']} features x {n_cols} strains")
    logging.info("\nDistribution:")
    for category in ['core', 'strict_core', 'soft_core', 'shell', 'cloud']:
        count = stats[category]
        pct = count/stats['total']*100
        logging.info(f"{category}: {count} ({pct:.1f}%)")
    
    logging.info("\nSharing patterns:")
    logging.info(f"Unique: {stats['singleton']} ({stats['singleton']/stats['total']*100:.1f}%)")
    logging.info(f"Shared: {stats['shared']} ({stats['shared']/stats['total']*100:.1f}%)")
    
    return stats

def get_proportional_sample(matrix: np.ndarray, stats: Dict[str, int], 
                           max_features: int = 500, core_boost: int = 3) -> np.ndarray:
    """Enhanced sampling with core gene prioritization and error handling for sparse matrices"""
    n_rows, n_cols = matrix.shape
    presence_counts = np.sum(matrix, axis=1)
    
    try:
        # Check matrix characteristics
        sparsity = 1.0 - (np.count_nonzero(matrix) / matrix.size)
        logging.info(f"Matrix sparsity: {sparsity:.2%}")
        
        # If extremely sparse, focus on rows with some variation
        if sparsity > 0.99:
            logging.warning(f"Matrix is extremely sparse, focusing on informative features")
            # Find rows with at least some presence
            informative_rows = np.where(presence_counts > 0)[0]
            if len(informative_rows) > 10:
                # If we have many informative rows, prioritize rows with variation
                min_present = max(1, int(0.1 * n_cols))  # At least in 10% of genomes
                max_present = min(n_cols - 1, int(0.9 * n_cols))  # At most in 90% of genomes
                varied_rows = np.where((presence_counts >= min_present) & (presence_counts <= max_present))[0]
                
                if len(varied_rows) > 0:
                    logging.info(f"Found {len(varied_rows)} features with interesting variation patterns")
                    return np.random.choice(varied_rows, min(max_features, len(varied_rows)), replace=False)
                else:
                    logging.info(f"Using {min(max_features, len(informative_rows))} informative features")
                    return np.random.choice(informative_rows, min(max_features, len(informative_rows)), replace=False)
        
        # Try to get core genes
        core_mask = presence_counts >= 0.90 * n_cols
        core_indices = np.where(core_mask)[0]
        
        # Get accessory genes 
        non_core_indices = np.where(~core_mask)[0]
        
        # Handle edge cases
        if len(core_indices) == 0 and len(non_core_indices) == 0:
            logging.warning("No valid features found, using random sampling")
            return np.random.choice(n_rows, min(max_features, n_rows), replace=False)
        
        if len(core_indices) == 0:
            logging.info("No core features found, sampling from accessory features")
            return np.random.choice(non_core_indices, min(max_features, len(non_core_indices)), replace=False)
            
        if len(non_core_indices) == 0:
            logging.info("All features are core, sampling randomly from core")
            return np.random.choice(core_indices, min(max_features, len(core_indices)), replace=False)
        
        # Calculate boosted weights for normal case
        total_weight = (len(core_indices) * core_boost) + len(non_core_indices)
        core_sample_size = min(int(max_features * (len(core_indices) * core_boost) / total_weight), 
                             len(core_indices))
        non_core_sample_size = max_features - core_sample_size
        
        # Sample with error handling
        sampled_core = core_indices if core_sample_size >= len(core_indices) else \
            np.random.choice(core_indices, core_sample_size, replace=False)
            
        sampled_non_core = np.random.choice(non_core_indices, 
                                          min(non_core_sample_size, len(non_core_indices)), 
                                          replace=False)
        
        combined = np.sort(np.concatenate([sampled_core, sampled_non_core]))
        logging.info(f"Sampled {len(sampled_core)} core and {len(sampled_non_core)} accessory features")
        return combined
    
    except Exception as e:
        logging.error(f"Sampling failed: {str(e)}")
        # Very simple fallback
        return np.arange(min(max_features, n_rows))

def plot_heatmap(matrix: np.ndarray, row_labels: List[str], col_labels: List[str],
                title: str, output_path: str) -> bool:
    """Enhanced heatmap visualization with improved clustering and error handling.
    Returns True if heatmap was successfully created, False otherwise."""
    
    try:
        # Check if the matrix is too sparse or has no variation
        sparsity = 1.0 - (np.count_nonzero(matrix) / matrix.size)
        
        # Check if matrix has any rows or columns
        if matrix.shape[0] == 0 or matrix.shape[1] == 0:
            logging.error(f"Matrix has 0 rows or columns: shape {matrix.shape}")
            return False
        
        # Check if matrix has no variation (all values the same)
        if np.all(matrix == matrix[0, 0]):
            logging.warning("Matrix has no variation (all values are identical), creating simple heatmap")
            plt.figure(figsize=(12, 10))
            sns.heatmap(matrix, cmap=sns.color_palette(["#FFFFFF", "#2E86C1"]), 
                       cbar=False, xticklabels=col_labels, yticklabels=False)
            plt.title(f"{title}\n(No variation)")
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            return True
            
        # Try to use fastcluster for better performance
        try:
            import fastcluster
            linkage_fn = fastcluster.linkage
            logging.info("Using fastcluster for faster hierarchical clustering")
        except ImportError:
            from scipy.cluster.hierarchy import linkage as linkage_fn
            logging.info("Using scipy for hierarchical clustering (consider installing fastcluster)")
        
        # Create figure
        plt.figure(figsize=(18, 12))
        
        # Create color palette for binary data
        cmap = sns.color_palette(["#FFFFFF", "#2E86C1"])  # White for 0, blue for 1
        
        # Create DataFrame for better labeling
        df = pd.DataFrame(matrix, index=row_labels, columns=col_labels)
        
        # Handle extremely sparse matrices (>99% zeros)
        if sparsity > 0.99:
            logging.warning(f"Matrix is extremely sparse ({sparsity:.2%}), using simplified visualization approach")
            
            # For extremely sparse data, try a simpler approach without hierarchical clustering
            plt.figure(figsize=(18, 12))
            sns.heatmap(df, cmap=cmap, cbar=False, xticklabels=True, yticklabels=False)
            plt.title(f"{title}\n(Simplified view - no clustering due to high sparsity: {sparsity:.2%})")
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logging.info(f"Created simplified heatmap at {output_path}")
            return True
        
        # Create clustermap with improved parameters
        g = sns.clustermap(
            df,
            cmap=cmap,
            method='ward',
            metric='euclidean',
            figsize=(18, 12),
            dendrogram_ratio=(0.15, 0.1),
            cbar_pos=None,  # Remove color bar
            yticklabels=False,
            xticklabels=True,
            tree_kws={'linewidths': 0.5}
        )
        
        # Improve axis labels
        g.ax_heatmap.set_xlabel("Strains", fontsize=10)
        g.ax_heatmap.set_ylabel("Features", fontsize=10)
        
        # Add presence/absence legend in a better position
        plt.text(1.02, 0.98, "Present (1)", transform=g.ax_heatmap.transAxes,
                color=cmap[1], fontsize=10)
        plt.text(1.02, 0.95, "Absent (0)", transform=g.ax_heatmap.transAxes,
                color=cmap[0], fontsize=10, bbox=dict(facecolor='white', edgecolor='none', alpha=0.7))
        
        # Rotate strain labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=8)
        
        # Save high-resolution output
        plt.savefig(output_path, dpi=300, bbox_inches='tight')
        plt.close()
        logging.info(f"Created clustered heatmap at {output_path}")
        return True
        
    except Exception as e:
        logging.error(f"Heatmap generation failed: {str(e)}")
        
        # Create a simple alternative visualization as fallback
        try:
            logging.info("Attempting to create simplified matrix visualization")
            plt.figure(figsize=(12, 10))
            
            # Take a small subset if matrix is too large
            if matrix.shape[0] > 100:
                idx = np.random.choice(matrix.shape[0], 100, replace=False)
                viz_matrix = matrix[idx]
                viz_row_labels = [row_labels[i] for i in idx]
            else:
                viz_matrix = matrix
                viz_row_labels = row_labels
                
            # Simple heatmap without clustering
            sns.heatmap(viz_matrix, cmap=sns.color_palette(["#FFFFFF", "#2E86C1"]), 
                       cbar=False, xticklabels=col_labels, yticklabels=False)
            
            plt.title(f"{title}\n(Simplified view - no clustering)")
            plt.tight_layout()
            plt.savefig(output_path, dpi=300)
            plt.close()
            logging.info(f"Created simplified visualization at {output_path}")
            return True
            
        except Exception as e2:
            logging.error(f"Even simplified visualization failed: {str(e2)}")
            return False

def load_labels(label_path: str) -> List[str]:
    """Load feature or strain labels from file"""
    try:
        with open(label_path, 'r') as f:
            return [line.strip() for line in f]
    except FileNotFoundError:
        logging.warning(f"Label file {label_path} not found, using generated labels")
        return []

def main(args) -> None:
    """Main analysis workflow with enhanced error handling"""
    setup_logging(args.output_dir)
    
    # Check for additional packages
    try:
        import fastcluster
        logging.info("Using fastcluster for improved clustering performance")
    except ImportError:
        logging.warning("fastcluster not found, using scipy for clustering (will be slower)")
    
    success = False
    for matrix_type in ['gene', 'allele']:
        try:
            # Load data and labels
            matrix_path = f"{args.input_dir}/{args.name}_strain_by_{matrix_type}.npz"
            label_path = f"{args.input_dir}/{args.name}_strain_by_{matrix_type}.npz.labels.txt"
            
            # Check if files exist
            if not os.path.exists(matrix_path):
                logging.error(f"Matrix file not found: {matrix_path}")
                continue
                
            # Load the matrix
            matrix = load_npz(matrix_path).toarray()
            logging.info(f"Loaded {matrix_type} matrix with shape {matrix.shape}")
            
            # Load or generate labels
            labels = load_labels(label_path) or [f"Feature_{i}" for i in range(matrix.shape[0])]
            strain_labels = [f"Strain_{i}" for i in range(matrix.shape[1])]
            
            # Analysis and sampling
            stats = analyze_pangenome(matrix, matrix_type.capitalize())
            sampled_indices = get_proportional_sample(matrix, stats, args.max_features, args.core_boost)
            
            if len(sampled_indices) == 0:
                logging.error(f"No features sampled for {matrix_type} matrix")
                continue
                
            # Visualization
            output_file = f"{args.output_dir}/{args.name}_{matrix_type}_heatmap.png"
            plot_title = f"{args.name} {matrix_type} Pangenome\n{sampled_indices.size} Features × {matrix.shape[1]} Strains"
            
            heatmap_success = plot_heatmap(
                matrix[sampled_indices],
                [labels[i] for i in sampled_indices],
                strain_labels,
                plot_title,
                output_file
            )
            
            if heatmap_success:
                success = True
                logging.info(f"Successfully created {matrix_type} heatmap")
            else:
                logging.error(f"Failed to create {matrix_type} heatmap")
            
        except Exception as e:
            logging.error(f"Failed processing {matrix_type} matrix: {str(e)}")
    
    if not success:
        # Create a dummy file to indicate the process ran, even if it didn't produce a heatmap
        with open(f"{args.output_dir}/{args.name}_no_heatmap.txt", 'w') as f:
            f.write("No heatmap could be generated due to data characteristics or errors.\n")
            f.write("See pangenome_viz.log for details.\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Pangenome Visualization Tool")
    parser.add_argument("--input-dir", required=True, help="Directory containing input matrices")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--name", required=True, help="Base name for files")
    parser.add_argument("--max-features", type=int, default=500, 
                        help="Maximum features to include in heatmap")
    parser.add_argument("--core-boost", type=int, default=3,
                        help="Boost factor for core gene sampling")
    
    args = parser.parse_args()
    main(args)
