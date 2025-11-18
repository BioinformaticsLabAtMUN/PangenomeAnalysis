#!/usr/bin/env python3
"""
Enhanced Pangenome Visualization with Organism Name Support

This script creates improved visualizations of pangenome data, including:
- Reduced gene sampling for better visibility (200 genes)
- Support for organism names instead of accession IDs
- Better visual styling for poster presentations
- Improved proportional sampling of genes by category
- Category color highlighting and legends

Usage:
  python pan_viz.py --input-dir DIR --output-dir DIR --name PREFIX [options]

Options:
  --max-features INT     Maximum number of genes to show (default: 200)
  --core-boost INT       Boost factor for core gene sampling (default: 4)
  --figure-width FLOAT   Width of figure in inches (default: 18)
  --figure-height FLOAT  Height of figure in inches (default: 14)
  --organism-mapping FILE Path to TSV file mapping accession IDs to organism names
"""

import os
# Set environment variable for headless rendering
os.environ['QT_QPA_PLATFORM'] = 'offscreen'

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
from scipy.sparse import load_npz
from scipy.cluster.hierarchy import linkage, dendrogram
import logging
import argparse
import warnings
from pathlib import Path
from matplotlib.patches import Patch
import re

# Improve font rendering
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans', 'Liberation Sans', 'Bitstream Vera Sans', 'sans-serif']
matplotlib.rcParams['axes.linewidth'] = 1.5  # Make axes bolder
matplotlib.rcParams['axes.labelsize'] = 14   # Larger axis labels
matplotlib.rcParams['xtick.major.width'] = 1.5  # Bold tick marks
matplotlib.rcParams['ytick.major.width'] = 1.5

# Setup logger
def setup_logging(output_dir):
    """Configure logging with file and console handlers"""
    log_file = os.path.join(output_dir, "pangenome_viz.log")
    
    # Create formatter
    formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
    
    # Setup file handler
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(formatter)
    
    # Setup console handler
    console_handler = logging.StreamHandler()
    console_handler.setFormatter(formatter)
    
    # Configure root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.INFO)
    
    # Remove any existing handlers
    for handler in root_logger.handlers[:]:
        root_logger.removeHandler(handler)
    
    # Add our handlers
    root_logger.addHandler(file_handler)
    root_logger.addHandler(console_handler)
    
    return logging.getLogger(__name__)

def load_matrix_and_labels(matrix_path, labels_path=None):
    """
    Load a sparse matrix and its labels
    
    Args:
        matrix_path: Path to npz file with scipy.sparse matrix
        labels_path: Path to text file with index and column names
                     If None, uses <matrix_path>.labels.txt
    
    Returns:
        matrix: Dense numpy array
        row_labels: List of row labels
        col_labels: List of column labels
    """
    try:
        # Load sparse matrix and convert to dense
        matrix = load_npz(matrix_path)
        
        # Determine labels path
        if labels_path is None:
            labels_path = f"{matrix_path}.labels.txt"
        
        # Load labels
        if os.path.exists(labels_path):
            with open(labels_path, 'r') as f:
                all_labels = [line.strip() for line in f]
            
            n_rows, n_cols = matrix.shape
            row_labels = all_labels[:n_rows]
            col_labels = all_labels[n_rows:]
            
            logging.info(f"Loaded {len(row_labels)} row labels and {len(col_labels)} column labels")
        else:
            logging.warning(f"Labels file {labels_path} not found. Using generic labels.")
            row_labels = [f"Feature_{i}" for i in range(matrix.shape[0])]
            col_labels = [f"Strain_{i}" for i in range(matrix.shape[1])]
        
        # Convert to dense array
        dense_matrix = matrix.toarray()
        logging.info(f"Loaded matrix with shape {dense_matrix.shape}")
        
        return dense_matrix, row_labels, col_labels
        
    except Exception as e:
        logging.error(f"Error loading matrix and labels: {str(e)}")
        raise

def load_organism_mapping(mapping_file):
    """
    Load organism name mapping from TSV file
    
    Args:
        mapping_file: Path to TSV file with accession to organism mapping
    
    Returns:
        Dictionary mapping accession IDs to organism names
    """
    try:
        if not os.path.exists(mapping_file):
            logging.warning(f"Organism mapping file {mapping_file} not found. Using original strain labels.")
            return {}
        
        # Load the mapping file
        mapping_df = pd.read_csv(mapping_file, sep='\t')
        
        # Check for required columns
        if 'accession_id' not in mapping_df.columns:
            logging.warning("Missing 'accession_id' column in mapping file. Using original strain labels.")
            return {}
        
        # Determine which name column to use
        if 'short_name' in mapping_df.columns:
            name_col = 'short_name'
        elif 'organism_name' in mapping_df.columns:
            name_col = 'organism_name'
        else:
            logging.warning("No name column found in mapping file. Using original strain labels.")
            return {}
        
        # Create mapping dictionary
        mapping = dict(zip(mapping_df['accession_id'], mapping_df[name_col]))
        logging.info(f"Loaded organism mapping for {len(mapping)} accessions")
        
        return mapping
    
    except Exception as e:
        logging.error(f"Error loading organism mapping: {str(e)}")
        return {}

def analyze_pangenome(matrix):
    """
    Analyze pangenome gene distribution with careful categorization
    
    Args:
        matrix: Binary presence/absence matrix (genes x strains)
    
    Returns:
        Dictionary with pangenome statistics and gene categorization
    """
    n_rows, n_cols = matrix.shape
    presence_counts = np.sum(matrix, axis=1)
    
    # Calculate presence percentages
    presence_pct = presence_counts / n_cols * 100
    
    # Define thresholds for gene categories - using standard definitions
    core_threshold = 95  # Present in >= 95% of genomes
    softcore_threshold = 85  # Present in >= 85% but < 95% of genomes
    shell_threshold = 15  # Present in >= 15% but < 85% of genomes
    cloud_threshold = 15  # Present in < 15% of genomes
    
    # Categorize genes
    core_mask = presence_pct >= core_threshold
    softcore_mask = (presence_pct >= softcore_threshold) & (presence_pct < core_threshold)
    shell_mask = (presence_pct >= cloud_threshold) & (presence_pct < softcore_threshold)
    cloud_mask = presence_pct < cloud_threshold
    
    # Check that all genes are assigned to exactly one category
    total_assigned = np.sum(core_mask) + np.sum(softcore_mask) + np.sum(shell_mask) + np.sum(cloud_mask)
    if total_assigned != n_rows:
        logging.error(f"Error in gene categorization! {total_assigned} genes assigned vs {n_rows} total genes")
    
    # Count genes in each category
    core_genes = np.sum(core_mask)
    strict_core_genes = np.sum(presence_counts == n_cols)  # Present in all genomes
    softcore_genes = np.sum(softcore_mask)
    shell_genes = np.sum(shell_mask)
    cloud_genes = np.sum(cloud_mask)
    
    # Calculate singleton and shared
    singleton_genes = np.sum(presence_counts == 1)  # Present in exactly 1 genome
    shared_genes = np.sum((presence_counts > 1) & (presence_counts < n_cols))  # Present in some but not all
    
    # Calculate total unique genes
    total_genes = matrix.shape[0]
    
    # Create category mask arrays for future use
    categories = {
        'core': core_mask,
        'softcore': softcore_mask,
        'shell': shell_mask,
        'cloud': cloud_mask
    }
    
    # Create category to color mapping - maintaining similar colors as requested
    category_colors = {
        'core': '#d62728',     # red for core genes
        'softcore': '#FFFDD0',  # cream for soft-core genes
        'shell': '#ff7f0e',    # orange for shell genes
        'cloud': '#1f77b4'     # blue for cloud genes
    }
    
    # Generate summary statistics
    stats = {
        'total': total_genes,
        'core': core_genes,
        'strict_core': strict_core_genes,
        'softcore': softcore_genes,
        'shell': shell_genes,
        'cloud': cloud_genes,
        'singleton': singleton_genes,
        'shared': shared_genes,
        'categories': categories,
        'category_colors': category_colors,
        'presence_counts': presence_counts,
        'presence_pct': presence_pct
    }
    
    # Log summary statistics
    logging.info("\nPangenome Analysis Summary:")
    logging.info(f"Total genes: {total_genes}")
    logging.info(f"Core genes (≥{core_threshold}%): {core_genes} ({core_genes/total_genes*100:.1f}%)")
    logging.info(f"  - Strict core (100%): {strict_core_genes} ({strict_core_genes/total_genes*100:.1f}%)")
    logging.info(f"Soft-core genes (≥{softcore_threshold}% to <{core_threshold}%): {softcore_genes} ({softcore_genes/total_genes*100:.1f}%)")
    logging.info(f"Shell genes (≥{cloud_threshold}% to <{softcore_threshold}%): {shell_genes} ({shell_genes/total_genes*100:.1f}%)")
    logging.info(f"Cloud genes (<{cloud_threshold}%): {cloud_genes} ({cloud_genes/total_genes*100:.1f}%)")
    logging.info(f"Singleton genes: {singleton_genes} ({singleton_genes/total_genes*100:.1f}%)")
    
    return stats

def sample_genes_for_visualization(matrix, stats, max_features=200, core_boost=4):
    """
    Sample genes for visualization with improved proportional representation
    
    Args:
        matrix: Binary presence/absence matrix (genes x strains)
        stats: Dictionary with pangenome statistics from analyze_pangenome()
        max_features: Maximum number of features to include (reduced to 200)
        core_boost: Weighting factor to boost sampling of core genes (increased to 4)
    
    Returns:
        Indices of selected genes
    """
    categories = stats['categories']
    
    # Count features in each category
    category_counts = {
        category: np.sum(mask) 
        for category, mask in categories.items()
    }
    
    logging.info(f"Category counts for sampling: {category_counts}")
    
    # Check if the categories make sense
    if sum(category_counts.values()) != matrix.shape[0]:
        logging.error(f"Category counts don't match total genes! {sum(category_counts.values())} vs {matrix.shape[0]}")
    
    # Apply core gene boost to weights
    # Increase core and softcore sampling for better visibility in the poster
    weights = {
        'core': category_counts['core'] * core_boost,
        'softcore': category_counts['softcore'] * (core_boost * 0.8),
        'shell': category_counts['shell'] * 0.8,  # Slightly reduce shell
        'cloud': category_counts['cloud'] * 0.4    # Further reduce cloud gene sampling
    }
    
    # Calculate proportional sampling
    total_weight = sum(weights.values())
    proportions = {
        category: weight / total_weight
        for category, weight in weights.items()
    }
    
    # Calculate number of features to sample from each category
    sample_counts = {
        category: min(int(max_features * prop), count)
        for category, prop, count in zip(
            proportions.keys(), 
            proportions.values(), 
            [category_counts[cat] for cat in proportions.keys()]
        )
    }
    
    # Ensure we don't exceed max_features
    total_samples = sum(sample_counts.values())
    if total_samples < max_features:
        # Distribute remaining slots proportionally to categories with capacity
        remaining = max_features - total_samples
        for category in sorted(sample_counts.keys(), 
                              key=lambda x: sample_counts[x] / category_counts[x]):
            available = category_counts[category] - sample_counts[category]
            if available > 0:
                add_count = min(remaining, available)
                sample_counts[category] += add_count
                remaining -= add_count
            if remaining == 0:
                break
    
    # Sample from each category
    selected_indices = []
    for category, count in sample_counts.items():
        if count > 0:
            category_indices = np.where(categories[category])[0]
            if len(category_indices) > 0:
                sampled = np.random.choice(category_indices, size=min(count, len(category_indices)), replace=False)
                selected_indices.extend(sampled)
    
    # Convert to numpy array and sort
    selected_indices = np.array(selected_indices)
    
    logging.info("\nProportional sampling results:")
    for category, count in sample_counts.items():
        logging.info(f"{category}: {count} genes")
    logging.info(f"Total selected: {len(selected_indices)} genes")
    
    return selected_indices

def create_enhanced_heatmap(matrix, row_labels, col_labels, stats, selected_indices, output_file, name, 
                           organism_mapping=None, include_categories=True, fig_width=18, fig_height=14, improved_layout=True):
    """
    Create an enhanced heatmap visualization with improved layout for poster presentation
    """
    try:
        # Extract selected data
        selected_matrix = matrix[selected_indices, :]
        selected_row_labels = [row_labels[i] for i in selected_indices]
        
        # Apply organism mapping to column labels if provided
        mapped_col_labels = col_labels.copy()
        if organism_mapping and len(organism_mapping) > 0:
            logging.info("Applying organism mapping to strain labels")
            for i, label in enumerate(col_labels):
                # Extract accession ID using regex pattern
                pattern = r'(GC[AF]_\d+\.\d+)'
                match = re.search(pattern, label)
                
                if match:
                    accession_id = match.group(1)
                    if accession_id in organism_mapping:
                        mapped_col_labels[i] = organism_mapping[accession_id]
                        logging.debug(f"Mapped {label} to {mapped_col_labels[i]}")
        
        # Get category masks for selected genes
        selected_categories = {}
        for category, mask in stats['categories'].items():
            selected_categories[category] = mask[selected_indices]
        
        # Convert to DataFrame for seaborn
        df = pd.DataFrame(selected_matrix, 
                        index=selected_row_labels, 
                        columns=mapped_col_labels)
        
        # Create custom color palette
        binary_cmap = sns.color_palette(['#FFFFFF', '#0047AB'])  # White and Cobalt Blue
        
        # Try to use fastcluster for better performance if available
        try:
            import fastcluster
            method = 'average'  # UPGMA linkage for cleaner dendrograms
            row_linkage = fastcluster.linkage(selected_matrix, method=method)
            col_linkage = fastcluster.linkage(selected_matrix.T, method=method)
            logging.info("Using fastcluster for hierarchical clustering")
        except ImportError:
            # Fall back to scipy
            method = 'average'  # UPGMA linkage for cleaner dendrograms
            row_linkage = linkage(selected_matrix, method=method)
            col_linkage = linkage(selected_matrix.T, method=method)
            logging.info("Using scipy for hierarchical clustering")

        # Create row colors for gene categories
        row_colors = []
        for i in range(len(selected_indices)):
            # Find which category this gene belongs to
            for category in ['core', 'softcore', 'shell', 'cloud']:
                if selected_categories[category][i]:
                    row_colors.append(stats['category_colors'][category])
                    break
        
        # Create the row_colors dataframe for seaborn
        row_colors_df = pd.DataFrame({'Category': row_colors}, index=selected_row_labels)
        
        # First, create a larger figure with more space at the bottom
        fig = plt.figure(figsize=(fig_width, fig_height + 2))  # Add extra height for legend
        
        # Create clustermap
        g = sns.clustermap(
            df,
            cmap=binary_cmap,
            row_linkage=row_linkage,
            col_linkage=col_linkage,
            figsize=(fig_width, fig_height),
            row_colors=row_colors_df,
            dendrogram_ratio=(0.12, 0.12),  # Smaller dendrograms
            cbar_pos=None,  # No color bar for binary data
            yticklabels=False,  # Too many gene labels to display
            xticklabels=True,  # Keep strain labels 
            tree_kws={'linewidths': 2.0, 'colors': '#333333'},  # Bolder dendrogram lines
            colors_ratio=0.05,  # Thinner category color bar
            row_cluster=True,  
            col_cluster=True  
        )
        
        # Add title
        plt.suptitle(f"{name} Pangenome Analysis\n{len(selected_indices)} genes × {matrix.shape[1]} strains", 
                    fontsize=20, y=1.02, fontweight='bold')
        
        # Improve axis labels
        g.ax_heatmap.set_xlabel("Strains", fontsize=18, fontweight='bold', labelpad=15)
        g.ax_heatmap.set_ylabel("Genes", fontsize=18, fontweight='bold', labelpad=15)
        g.ax_heatmap.tick_params(axis='x', pad=2) 
        
        # Beautify strain labels
        plt.setp(g.ax_heatmap.get_xticklabels(), rotation=45, ha='right', fontsize=10, fontweight='bold')
        
        # Create a dedicated axis for the legend at the bottom
        legend_ax = fig.add_axes([0.2, 0.02, 0.6, 0.03])  # [left, bottom, width, height]
        legend_ax.axis('off')  # Turn off axis
        
        # Add category legend to the dedicated axis
        legend_elements = [
            Patch(facecolor=stats['category_colors']['core'], edgecolor=None, label=f'Core (≥95%)'),
            Patch(facecolor=stats['category_colors']['softcore'], edgecolor=None, label=f'Soft-core (85-95%)'),
            Patch(facecolor=stats['category_colors']['shell'], edgecolor=None, label=f'Shell (15-85%)'),
            Patch(facecolor=stats['category_colors']['cloud'], edgecolor=None, label=f'Cloud (<15%)')
        ]
        
        # Add legend with horizontal layout
        legend_ax.legend(handles=legend_elements, loc='center', ncol=4, fontsize=14, 
                       frameon=False, title="Gene Categories", title_fontsize=16)
        
        # Save with tight layout to ensure everything fits
        plt.tight_layout()
        #plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.savefig(output_file, dpi=300, bbox_inches='tight', transparent=True)

        plt.close()
        
        logging.info(f"Successfully created enhanced heatmap at {output_file}")
        return True
    
    except Exception as e:
        logging.error(f"Error creating enhanced heatmap: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

def create_pangenome_composition_pie(stats, output_file, name):
    """
    Create a simple, clean pie chart showing the composition of the pangenome
    """
    try:
        plt.figure(figsize=(10, 8))
        
        # Data for pie chart
        categories = ['Core', 'Soft-core', 'Shell', 'Cloud']
        sizes = [
            stats['core'], 
            stats['softcore'], 
            stats['shell'], 
            stats['cloud']
        ]
        colors = [
            stats['category_colors']['core'],
            stats['category_colors']['softcore'],
            stats['category_colors']['shell'],
            stats['category_colors']['cloud']
        ]
        
        # Add percentage to labels
        total = sum(sizes)
        labels = [f'{cat}\n{size} ({size/total*100:.1f}%)' for cat, size in zip(categories, sizes)]
        
        # Create a simple pie chart with no borders or effects
        wedges, texts = plt.pie(
            sizes, 
            colors=colors,
            startangle=90,
            # Remove all borders and styling
            wedgeprops={'edgecolor': 'none', 'linewidth': 0}
        )
        
        # Add centered legend
        plt.legend(
            wedges, 
            labels,
            loc="center",
            bbox_to_anchor=(0.5, 0.5),
            fontsize=14,
            frameon=True
        )
        
        plt.title(f'{name} Pangenome Composition', fontsize=18, fontweight='bold')
        
        # Save pie chart
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        logging.info(f"Created pangenome composition pie chart at {output_file}")
        return True
    
    except Exception as e:
        logging.error(f"Error creating pie chart: {str(e)}")
        return False

def extract_mapping_from_data_summary(data_summary_path):
    """
    Extract simplified mapping from data_summary.tsv file - genus initial + species only
    
    Args:
        data_summary_path: Path to data_summary.tsv file
    
    Returns:
        Dictionary mapping accession IDs to simplified organism names
    """
    try:
        if not os.path.exists(data_summary_path):
            logging.warning(f"Data summary file {data_summary_path} not found")
            return {}
        
        # Load the data summary file
        df = pd.read_csv(data_summary_path, sep='\t')
        
        # Check for required columns
        required_cols = ['Organism Scientific Name', 'Assembly Accession']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            logging.warning(f"Missing columns in data summary file: {missing_cols}")
            return {}
        
        # Create mapping dictionary with simplified names
        mapping = {}
        for _, row in df.iterrows():
            org_name = row['Organism Scientific Name']
            accession = row['Assembly Accession']
            
            if pd.isna(org_name) or pd.isna(accession):
                continue
                
            # Create short name (genus initial + species only)
            parts = org_name.split()
            if len(parts) >= 2:
                genus_initial = parts[0][0]
                species = parts[1]
                simplified_name = f"{genus_initial}. {species}"
            else:
                simplified_name = org_name
            
            mapping[accession] = simplified_name
        
        logging.info(f"Extracted {len(mapping)} simplified organism mappings from data summary file")
        return mapping
    
    except Exception as e:
        logging.error(f"Error extracting mapping from data summary: {str(e)}")
        return {}

def main():
    # Parse command line arguments
    parser = argparse.ArgumentParser(description="Enhanced Pangenome Visualization for Poster")
    parser.add_argument("--input-dir", required=True, help="Directory containing input matrices")
    parser.add_argument("--output-dir", required=True, help="Output directory for results")
    parser.add_argument("--name", required=True, help="Base name for files")
    parser.add_argument("--max-features", type=int, default=200, help="Maximum features to include in heatmap")
    parser.add_argument("--core-boost", type=int, default=4, help="Boost factor for core gene sampling")
    parser.add_argument("--figure-width", type=int, default=18, help="Figure width in inches")
    parser.add_argument("--figure-height", type=int, default=14, help="Figure height in inches")
    parser.add_argument("--improved-layout", action="store_true", default=True, help="Use improved layout for poster")
    parser.add_argument("--organism-mapping", help="Path to TSV file mapping accession IDs to organism names")
    parser.add_argument("--data-summary", help="Path to data_summary.tsv file to extract organism mapping")
    
    args = parser.parse_args()
    
    # Setup logging
    logger = setup_logging(args.output_dir)
    
    try:
        # Load organism mapping - try different sources
        organism_mapping = {}
        
        # 1. Try using provided organism mapping file first
        if args.organism_mapping and os.path.exists(args.organism_mapping):
            organism_mapping = load_organism_mapping(args.organism_mapping)
        
        # 2. If no mapping loaded yet, try using data_summary.tsv
        if not organism_mapping and args.data_summary:
            organism_mapping = extract_mapping_from_data_summary(args.data_summary)
        
        # 3. If still no mapping, try default data_summary.tsv location
        if not organism_mapping:
            default_data_summary = "/home/saba/Documents/thesis_project/input/data/data_summary.tsv"
            if os.path.exists(default_data_summary):
                logging.info(f"Using default data_summary.tsv at {default_data_summary}")
                organism_mapping = extract_mapping_from_data_summary(default_data_summary)
        
        # Process only gene matrix (skip allele matrix)
        matrix_type = 'gene'
        try:
            # Matrix path
            matrix_path = os.path.join(args.input_dir, f"{args.name}_strain_by_{matrix_type}.npz")
            
            # Check if matrix exists
            if not os.path.exists(matrix_path):
                logging.error(f"Matrix file not found: {matrix_path}")
                # Create a failure file
                with open(os.path.join(args.output_dir, f"{args.name}_no_heatmap.txt"), 'w') as f:
                    f.write(f"Matrix file not found: {matrix_path}\n")
                return
            
            # Load matrix and labels
            logging.info(f"Processing {matrix_type} matrix: {matrix_path}")
            matrix, row_labels, col_labels = load_matrix_and_labels(matrix_path)
            
            # Analyze pangenome
            stats = analyze_pangenome(matrix)
            
            # Sample genes for visualization
            selected_indices = sample_genes_for_visualization(
                matrix, stats, args.max_features, args.core_boost)
            
            # Create heatmap with categories - main visualization for poster
            heatmap_with_cats_file = os.path.join(args.output_dir, f"{args.name}_{matrix_type}_heatmap_poster.png")
            create_enhanced_heatmap(
                matrix, row_labels, col_labels, stats, 
                selected_indices, heatmap_with_cats_file, 
                f"{args.name}", organism_mapping=organism_mapping,
                include_categories=True,
                fig_width=args.figure_width, fig_height=args.figure_height,
                improved_layout=args.improved_layout
            )
            
            # Create pangenome composition pie chart
            pie_chart_file = os.path.join(args.output_dir, f"{args.name}_{matrix_type}_pangenome_composition.png")
            create_pangenome_composition_pie(stats, pie_chart_file, f"{args.name}")
            
            logging.info(f"Completed {matrix_type} matrix visualization for poster")
            
        except Exception as e:
            logging.error(f"Error processing {matrix_type} matrix: {str(e)}")
            import traceback
            traceback.print_exc()
            # Create a failure file
            with open(os.path.join(args.output_dir, f"{args.name}_no_heatmap.txt"), 'w') as f:
                f.write(f"Visualization failed: {str(e)}\n")
        
    except Exception as e:
        logging.error(f"Visualization process failed: {str(e)}")
        # Create a file indicating failure
        with open(os.path.join(args.output_dir, f"{args.name}_viz_failed.txt"), 'w') as f:
            f.write(f"Visualization failed: {str(e)}\n")

if __name__ == "__main__":
    main()
