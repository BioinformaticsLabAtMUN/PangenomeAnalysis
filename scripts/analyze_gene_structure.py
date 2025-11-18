#!/usr/bin/env python3
"""
Gene-Level Structural Analysis
Focuses on quantitative/structural aspects of pangenome without functional annotations.
Handles: core coverage, unique genes, missing genes, gene frequencies, outlier detection.
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path
import sys
import warnings
from collections import defaultdict, Counter
from scipy import stats
import warnings
warnings.filterwarnings('ignore')

# Optional imports
try:
    import plotly.express as px
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
    PLOTLY_AVAILABLE = True
except ImportError:
    PLOTLY_AVAILABLE = False
    print("Plotly not available - skipping interactive plots")

# Import custom sparse utilities
try:
    import sparse_utils
except ImportError:
    print("ERROR: sparse_utils module not found. Please ensure it's in your PYTHONPATH.")
    sys.exit(1)

def load_data(gene_matrix, gene_labels, core_genes_file, metadata_file=None, annotations_file=None):
    """Load all required data files for structural analysis."""
    print("Loading data files for structural analysis...")
    
    # Load gene matrix
    df_genes = sparse_utils.read_lsdf(gene_matrix, gene_labels)
    print(f"Loaded gene matrix: {df_genes.shape}")
    
    # Load core genes
    with open(core_genes_file, 'r') as f:
        core_genes = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(core_genes)} core genes")
    
    # Load metadata if provided
    metadata = None
    if metadata_file and Path(metadata_file).exists():
        try:
            metadata = pd.read_csv(metadata_file, sep='\t')
            print(f"Loaded metadata for {len(metadata)} proteomes")
        except Exception as e:
            print(f"Warning: Could not load metadata: {e}")
    
    # Load annotations if provided
    annotations = None
    if annotations_file and Path(annotations_file).exists():
        try:
            annotations = pd.read_csv(annotations_file, sep='\t')
            print(f"Loaded gene annotations for {len(annotations)} genes")
        except Exception as e:
            print(f"Warning: Could not load annotations: {e}")

    if annotations is None:
        print("No gene annotations found - will use cluster IDs")
    
    return df_genes, core_genes, metadata, annotations

def get_full_organism_name(proteome_id, metadata):
    """Get complete organism name from metadata - preserves full names."""
    if metadata is None:
        return proteome_id
    
    # Handle case where metadata might be a string (error condition)
    if isinstance(metadata, str):
        return metadata
    
    try:
        matching_row = metadata[metadata['Proteome Id'] == proteome_id]
        if not matching_row.empty:
            organism_name = matching_row.iloc[0]['Organism']
            return organism_name  # Return complete name as-is
    except (KeyError, TypeError, AttributeError):
        # If there's any issue with the metadata lookup, return the proteome_id
        pass
    
    return proteome_id

def analyze_core_genome_structure(df_genes, core_genes, metadata):
    """Comprehensive core genome structural analysis."""
    print("\n=== Core Genome Structural Analysis ===")
    
    # Filter to core genes in matrix
    core_genes_in_matrix = [g for g in core_genes if g in df_genes.index]
    print(f"Analyzing {len(core_genes_in_matrix)} core genes across {len(df_genes.columns)} strains")
    
    # Get core genome subset
    core_df = df_genes.labelslice(indices=core_genes_in_matrix)
    
    # Calculate coverage per strain
    coverage_counts = core_df.sum(axis='columns')
    total_core_genes = len(core_genes_in_matrix)
    coverage_percentages = (coverage_counts / total_core_genes) * 100
    
    # Calculate missing genes per strain
    missing_counts = total_core_genes - coverage_counts
    
    # Create comprehensive results dataframe
    results_df = pd.DataFrame({
        'strain_id': list(df_genes.columns),
        'core_genes_present': coverage_counts.astype(int),
        'core_genes_missing': missing_counts.astype(int),
        'total_core_genes': total_core_genes,
        'core_coverage_percent': coverage_percentages,
        'total_genes_in_strain': df_genes.sum(axis='columns').astype(int)
    })
    
    # Add metadata integration
    if metadata is not None:
        # Create metadata mapping
        metadata_dict = {}
        for _, row in metadata.iterrows():
            proteome_id = row['Proteome Id']
            metadata_dict[proteome_id] = {
                'organism': row['Organism'],
                'protein_count': row.get('Protein count', np.nan)
            }
        
        # Add organism names and protein counts
        organism_names = []
        protein_counts = []
        for strain_id in results_df['strain_id']:
            if strain_id in metadata_dict:
                organism_names.append(metadata_dict[strain_id]['organism'])
                protein_counts.append(metadata_dict[strain_id]['protein_count'])
            else:
                organism_names.append(strain_id)
                protein_counts.append(np.nan)
        
        results_df['organism'] = organism_names
        results_df['protein_count'] = protein_counts
        
        # Calculate genome size efficiency
        results_df['genes_per_protein_ratio'] = results_df['total_genes_in_strain'] / results_df['protein_count']
    
    # Statistical summary
    print(f"Core Coverage Statistics:")
    print(f"  Mean coverage: {np.mean(coverage_percentages):.2f}%")
    print(f"  Median coverage: {np.median(coverage_percentages):.2f}%")
    print(f"  Range: {np.min(coverage_percentages):.1f}% - {np.max(coverage_percentages):.1f}%")
    print(f"  Standard deviation: {np.std(coverage_percentages):.2f}%")
    
    # Identify structural outliers
    low_coverage_threshold = np.percentile(coverage_percentages, 10)  # Bottom 10%
    high_missing_threshold = np.percentile(missing_counts, 90)       # Top 10% missing
    
    structural_outliers = results_df[
        (results_df['core_coverage_percent'] <= low_coverage_threshold) |
        (results_df['core_genes_missing'] >= high_missing_threshold)
    ]
    
    print(f"Structural outliers identified: {len(structural_outliers)} strains")
    
    return results_df, structural_outliers

def analyze_gene_frequency_distribution(df_genes, core_genes):
    """Analyze gene frequency distribution across the pangenome."""
    print("\n=== Gene Frequency Distribution Analysis ===")
    
    n_strains = len(df_genes.columns)
    
    # Calculate gene frequencies
    gene_counts = df_genes.sum(axis='index')
    gene_frequencies = (gene_counts / n_strains) * 100
    
    # Categorize genes by frequency
    core_genes_set = set(core_genes)
    
    # Gene categories
    unique_genes = df_genes.index[gene_counts == 1]
    rare_genes = df_genes.index[(gene_counts > 1) & (gene_counts <= 0.05 * n_strains)]
    accessory_genes = df_genes.index[(gene_counts > 0.05 * n_strains) & (gene_counts < 0.95 * n_strains)]
    core_genes_found = df_genes.index[gene_counts >= 0.95 * n_strains]
    
    # Filter out core genes from accessory for cleaner categorization
    accessory_genes_clean = [g for g in accessory_genes if g not in core_genes_set]
    
    # FIXED: Calculate other_shared_genes (rare + accessory)
    other_shared_genes = len(rare_genes) + len(accessory_genes_clean)
    
    freq_analysis = {
        'total_genes': len(df_genes.index),
        'unique_genes': len(unique_genes),
        'rare_genes': len(rare_genes),
        'accessory_genes': len(accessory_genes_clean),
        'other_shared_genes': other_shared_genes,  # FIXED: Added this key
        'core_genes_by_frequency': len(core_genes_found),
        'core_genes_ml_defined': len(core_genes),
        'core_genes': len(core_genes),  # FIXED: Added this key as well for consistency
        'gene_frequencies': gene_frequencies,
        'gene_counts': gene_counts
    }
    
    print(f"Gene frequency categories:")
    print(f"  Total genes: {freq_analysis['total_genes']:,}")
    print(f"  Unique genes (1 strain): {freq_analysis['unique_genes']:,} ({freq_analysis['unique_genes']/freq_analysis['total_genes']*100:.1f}%)")
    print(f"  Rare genes (2-{int(0.05*n_strains)} strains): {freq_analysis['rare_genes']:,} ({freq_analysis['rare_genes']/freq_analysis['total_genes']*100:.1f}%)")
    print(f"  Accessory genes: {freq_analysis['accessory_genes']:,} ({freq_analysis['accessory_genes']/freq_analysis['total_genes']*100:.1f}%)")
    print(f"  Core genes (≥95% strains): {freq_analysis['core_genes_by_frequency']:,}")
    print(f"  ML-defined core genes: {freq_analysis['core_genes_ml_defined']:,}")
    
    return freq_analysis, unique_genes, rare_genes, accessory_genes_clean

def analyze_unique_genes_per_strain(df_genes, unique_genes, metadata):
    """Analyze unique gene distribution per strain."""
    print("\n=== Unique Genes Per Strain Analysis ===")
    
    unique_stats = []
    gene_matrix_csc = df_genes.data.tocsc()
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        strain_column = gene_matrix_csc.getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(df_genes.index[i] for i in present_gene_indices)
        
        # Count unique genes for this strain
        strain_unique_count = len(present_genes.intersection(set(unique_genes)))
        
        organism = get_full_organism_name(strain_id, metadata)
        
        unique_stats.append({
            'strain_id': strain_id,
            'organism': organism,
            'unique_genes': strain_unique_count,
            'total_genes': len(present_genes),
            'unique_gene_percent': (strain_unique_count / len(present_genes) * 100) if len(present_genes) > 0 else 0
        })
    
    unique_df = pd.DataFrame(unique_stats)
    
    print(f"Unique genes per strain statistics:")
    print(f"  Mean: {unique_df['unique_genes'].mean():.0f} genes")
    print(f"  Range: {unique_df['unique_genes'].min()}-{unique_df['unique_genes'].max()} genes")
    print(f"  Mean percentage: {unique_df['unique_gene_percent'].mean():.1f}%")
    
    return unique_df

def identify_missing_core_genes_structural(df_genes, core_genes):
    """Identify missing core genes with focus on structural patterns."""
    print("\n=== Missing Core Genes Structural Analysis ===")
    
    core_genes_in_matrix = [g for g in core_genes if g in df_genes.index]
    core_df = df_genes.labelslice(indices=core_genes_in_matrix)
    
    missing_genes_per_strain = {}
    core_gene_absence_frequency = Counter()
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        strain_column = core_df.data.tocsc().getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(core_genes_in_matrix[i] for i in present_gene_indices)
        
        missing_genes = set(core_genes_in_matrix) - present_genes
        missing_genes_per_strain[strain_id] = missing_genes
        
        # Count frequency of each missing gene
        for gene in missing_genes:
            core_gene_absence_frequency[gene] += 1
    
    # Calculate missing patterns
    missing_counts = {strain_id: len(missing_genes) for strain_id, missing_genes in missing_genes_per_strain.items()}
    
    print(f"Missing core genes analysis:")
    print(f"  Strains with missing core genes: {len([c for c in missing_counts.values() if c > 0])}")
    print(f"  Mean missing genes: {np.mean(list(missing_counts.values())):.1f}")
    print(f"  Max missing genes: {max(missing_counts.values())}")
    
    return missing_genes_per_strain, missing_counts, core_gene_absence_frequency

def create_structural_visualizations(core_results, freq_analysis, unique_df, missing_counts, 
                                   core_gene_absence_frequency, metadata, annotations, output_dir, max_strains_plot=50):
    """Create essential structural visualizations only."""
    print("\n=== Creating Essential Structural Visualizations ===")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Set style
    plt.style.use('default')
    
    # 1. Core Coverage Analysis with Interactive HTML + Static versions
    create_core_coverage_plots(core_results, output_dir)
    
    # 2. Simplified Pangenome Composition (only pie chart)
    create_simple_pangenome_composition(freq_analysis, output_dir)
    
    # 3. All Strains Missing Core Genes (the important one you wanted)
    create_all_strains_missing_plot(missing_counts, metadata, output_dir, max_strains_plot)
    
    # 4. Detailed Missing Core Genes Plot (enhanced version)
    create_enhanced_missing_genes_detailed_plot(core_gene_absence_frequency, annotations, output_dir)
    
    # 5. NEW: Create missing genes summary table
    create_missing_genes_summary_table(core_gene_absence_frequency, annotations, output_dir)
    
    # 6. NEW: Create unique genes distribution plot
    create_unique_genes_distribution_plot(unique_df, output_dir)

def create_missing_genes_summary_table(core_gene_absence_frequency, annotations, output_dir):
    """Create a comprehensive summary table of missing core genes with detailed annotations."""
    
    if not core_gene_absence_frequency:
        return
    
    # Get top missing genes
    top_missing = core_gene_absence_frequency.most_common(30)
    
    # Create comprehensive table data
    table_data = []
    
    for gene_cluster, missing_count in top_missing:
        # Initialize row with defaults
        row = {
            'Gene Name': gene_cluster,
            'Protein Function': '',
            'Biological Process': '',
            'Missing Count': missing_count,
            'Missing %': round((missing_count / 153) * 100, 1)  # Assuming 153 total strains
        }
        
        # Look up annotations
        if annotations is not None:
            gene_info = annotations[annotations['gene'] == gene_cluster]
            if not gene_info.empty:
                gene_row = gene_info.iloc[0]
                
                # Get protein product name
                if pd.notna(gene_row.get('product', '')):
                    product_name = str(gene_row['product']).strip()
                    if product_name and product_name != 'nan' and product_name.lower() != 'deleted':
                        row['Gene Name'] = product_name
                        row['Protein Function'] = product_name
                
                # Get biological process from go_process
                if pd.notna(gene_row.get('go_process', '')):
                    bio_process = str(gene_row['go_process']).strip()
                    if bio_process and bio_process != 'nan':
                        row['Biological Process'] = bio_process
        
        table_data.append(row)
    
    # Create DataFrame and save
    df_table = pd.DataFrame(table_data)
    
    # Save as TSV
    df_table.to_csv(output_dir / 'missing_core_genes_summary_table.tsv', sep='\t', index=False)
    
    # Create a formatted HTML table for better visualization
    html_content = f"""
    <!DOCTYPE html>
    <html>
    <head>
        <title>Missing Core Genes Summary</title>
        <style>
            body {{ font-family: Arial, sans-serif; margin: 20px; }}
            h1 {{ color: #2c3e50; text-align: center; }}
            table {{ border-collapse: collapse; width: 100%; margin: 20px 0; }}
            th {{ background-color: #3498db; color: white; padding: 12px; text-align: left; font-weight: bold; }}
            td {{ padding: 10px; border-bottom: 1px solid #ddd; }}
            tr:nth-child(even) {{ background-color: #f2f2f2; }}
            tr:hover {{ background-color: #e8f4fd; }}
            .gene-name {{ font-weight: bold; color: #2980b9; }}
            .missing-high {{ background-color: #ffebee; color: #c62828; font-weight: bold; }}
            .missing-medium {{ background-color: #fff3e0; color: #f57c00; font-weight: bold; }}
            .missing-low {{ background-color: #e3f2fd; color: #1976d2; font-weight: bold; }}
        </style>
    </head>
    <body>
        <h1>Most Frequently Missing Core Genes in Streptomyces Pangenome</h1>
        <table>
            <thead>
                <tr>
                    <th>Gene Name</th>
                    <th>Protein Function</th>
                    <th>Biological Process</th>
                    <th>Missing Count</th>
                    <th>Missing %</th>
                </tr>
            </thead>
            <tbody>
    """
    
    for _, row in df_table.iterrows():
        # Color code missing percentage
        missing_pct = row['Missing %']
        if missing_pct >= 10:
            missing_class = 'missing-high'
        elif missing_pct >= 5:
            missing_class = 'missing-medium'
        else:
            missing_class = 'missing-low'
        
        html_content += f"""
                <tr>
                    <td class="gene-name">{row['Gene Name']}</td>
                    <td>{row['Protein Function']}</td>
                    <td>{row['Biological Process']}</td>
                    <td class="{missing_class}">{row['Missing Count']}</td>
                    <td class="{missing_class}">{row['Missing %']}%</td>
                </tr>
        """
    
    html_content += """
            </tbody>
        </table>
    </body>
    </html>
    """
    
    # Save HTML version
    with open(output_dir / 'missing_core_genes_summary_table.html', 'w') as f:
        f.write(html_content)
    
    print(f"✅ Missing genes summary table created (TSV + HTML format)")

def create_unique_genes_distribution_plot(unique_df, output_dir):
    """Create unique genes distribution plot - histogram + top 10 bar chart."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # Left plot: Histogram of unique genes distribution
    unique_counts = unique_df['unique_genes']
    mean_unique = unique_counts.mean()
    
    ax1.hist(unique_counts, bins=25, alpha=0.7, color='lightcoral', edgecolor='darkred', linewidth=1.5)
    ax1.axvline(mean_unique, color='red', linestyle='--', linewidth=2, label=f'Mean: {mean_unique:.0f}')
    ax1.set_xlabel('Number of Unique Genes per Genome', fontsize=12, fontweight='bold')
    ax1.set_ylabel('Number of Streptomyces Strains', fontsize=12, fontweight='bold')
    ax1.set_title('Distribution of Unique Genes Across Streptomyces', fontsize=14, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.legend(fontsize=12)
    
    # Right plot: Top 10 strains with most unique genes
    top_10 = unique_df.nlargest(10, 'unique_genes')
    
    # Truncate long organism names for better display
    display_names = []
    for organism in top_10['organism']:
        if len(organism) > 50:
            display_names.append(organism[:47] + "...")
        else:
            display_names.append(organism)
    
    bars = ax2.barh(range(len(top_10)), top_10['unique_genes'], 
                   color='lightblue', alpha=0.8, edgecolor='darkblue', linewidth=1.5)
    
    ax2.set_yticks(range(len(top_10)))
    ax2.set_yticklabels(display_names, fontsize=10)
    ax2.set_xlabel('Number of Unique Genes', fontsize=12, fontweight='bold')
    ax2.set_title('Top 10 Streptomyces with Most Unique Genes', fontsize=14, fontweight='bold')
    ax2.grid(True, alpha=0.3, axis='x')
    
    # Add value labels on bars
    for i, (bar, count) in enumerate(zip(bars, top_10['unique_genes'])):
        width = bar.get_width()
        ax2.text(width + max(top_10['unique_genes'])*0.01, bar.get_y() + bar.get_height()/2,
                f'{count}', ha='left', va='center', fontsize=10, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'unique_genes_distribution.png', dpi=300, bbox_inches='tight')
    plt.close()
    print(f"✅ Unique genes distribution plot saved")

def create_core_coverage_plots(core_results, output_dir):
    """Create core coverage analysis plots - Interactive HTML + Static versions."""
    
    # 1. Interactive HTML plot with outlier species names (like original)
    if PLOTLY_AVAILABLE and 'protein_count' in core_results.columns:
        plot_data = core_results.dropna(subset=['protein_count'])
        
        if len(plot_data) > 0:
            # Calculate trend line and identify outliers
            z = np.polyfit(plot_data['protein_count'], plot_data['core_coverage_percent'], 1)
            p = np.poly1d(z)
            predicted_coverage = p(plot_data['protein_count'])
            residuals = plot_data['core_coverage_percent'] - predicted_coverage
            residual_threshold = np.percentile(residuals, 15)  # Bottom 15% as outliers
            
            plot_data = plot_data.copy()
            plot_data['is_outlier'] = residuals < residual_threshold
            plot_data['residual'] = residuals
            
            # Create hover text with FULL organism names for ALL points
            hover_text = []
            for idx, row in plot_data.iterrows():
                full_organism_name = row.get('organism', row['strain_id'])
                status = "OUTLIER" if row['is_outlier'] else "Normal"
                hover_text.append(
                    f"<b style='font-size:14px'>{full_organism_name}</b><br><br>" +
                    f"<b>Status:</b> {status}<br>" +
                    f"<b>Protein Count:</b> {int(row['protein_count']):,}<br>" +
                    f"<b>Core Coverage:</b> {row['core_coverage_percent']:.1f}%<br>" +
                    f"<b>Missing Core Genes:</b> {row['core_genes_missing']}<br>" +
                    f"<b>Total Genes:</b> {row['total_genes_in_strain']:,}<br>" +
                    f"<b>Deviation from Expected:</b> {row['residual']:.2f}%"
                )
            
            fig = go.Figure()
            
            # Add normal points (blue)
            normal_data = plot_data[~plot_data['is_outlier']]
            fig.add_trace(go.Scatter(
                x=normal_data['protein_count'],
                y=normal_data['core_coverage_percent'],
                mode='markers',
                marker=dict(size=8, color='steelblue', opacity=0.7, line=dict(width=1, color='darkblue')),
                text=[hover_text[i] for i in normal_data.index],
                hovertemplate='%{text}<extra></extra>',
                name=f'Normal Strains ({len(normal_data)})',
                showlegend=True
            ))
            
            # Add outlier points (red, larger)
            outlier_data = plot_data[plot_data['is_outlier']]
            fig.add_trace(go.Scatter(
                x=outlier_data['protein_count'],
                y=outlier_data['core_coverage_percent'],
                mode='markers',
                marker=dict(size=12, color='red', line=dict(width=2, color='darkred')),
                text=[hover_text[i] for i in outlier_data.index],
                hovertemplate='%{text}<extra></extra>',
                name=f'Outlier Strains ({len(outlier_data)})',
                showlegend=True
            ))
            
            # Add trend line
            correlation = plot_data['protein_count'].corr(plot_data['core_coverage_percent'])
            fig.add_trace(go.Scatter(
                x=plot_data['protein_count'],
                y=p(plot_data['protein_count']),
                mode='lines',
                line=dict(dash='dash', color='red', width=3),
                name=f'Trend Line (r={correlation:.3f})',
                hovertemplate=f'Correlation: {correlation:.3f}<extra></extra>',
                showlegend=True
            ))
            
            fig.update_layout(
                title={
                    'text': '<b>Interactive: Core Coverage vs Protein Count</b><br><span style="font-size:14px">🖱️ Hover over points to see full species names</span>',
                    'x': 0.5,
                    'xanchor': 'center',
                    'font': {'size': 18}
                },
                xaxis_title='<b>Total Protein Count</b>',
                yaxis_title='<b>Core Genome Coverage (%)</b>',
                width=1200,
                height=700,
                hovermode='closest',
                plot_bgcolor='rgba(248,248,248,0.8)',
                paper_bgcolor='white',
                font=dict(size=12),
                legend=dict(
                    yanchor="top", y=0.99, xanchor="left", x=0.01,
                    bgcolor="rgba(255,255,255,0.9)", bordercolor="rgba(0,0,0,0.3)", borderwidth=2
                )
            )
            
            fig.update_xaxes(showgrid=True, gridwidth=1, gridcolor='rgba(200,200,200,0.4)')
            fig.update_yaxes(showgrid=True, gridwidth=1, gridcolor='rgba(200,200,200,0.4)')
            
            # Save interactive plot
            fig.write_html(output_dir / 'interactive_core_coverage_with_outlier_names.html',
                          config={'displayModeBar': True, 'displaylogo': False})
            print(f"✅ Interactive plot saved with {len(outlier_data)} outliers - hover for full species names!")
    
    # 2. Static version (exactly like original blue scatter)
    if 'protein_count' in core_results.columns:
        plot_data = core_results.dropna(subset=['protein_count'])
        
        if len(plot_data) > 0:
            fig, ax = plt.subplots(1, 1, figsize=(10, 8))
            
            # Calculate correlation and trend line
            correlation = plot_data['protein_count'].corr(plot_data['core_coverage_percent'])
            z = np.polyfit(plot_data['protein_count'], plot_data['core_coverage_percent'], 1)
            p = np.poly1d(z)
            
            # Blue scatter plot (exactly like original)
            ax.scatter(plot_data['protein_count'], plot_data['core_coverage_percent'], 
                      alpha=0.6, s=60, color='steelblue', edgecolors='darkblue', linewidth=0.5)
            
            # Add trend line
            ax.plot(plot_data['protein_count'], p(plot_data['protein_count']), 
                   "r--", alpha=0.8, linewidth=2, label=f'Trend line (r={correlation:.3f})')
            
            ax.set_xlabel('Total Protein Count', fontsize=14, fontweight='bold')
            ax.set_ylabel('Core Genome Coverage (%)', fontsize=14, fontweight='bold')
            ax.set_title('Core Coverage vs Protein Count\nStreptomyces Pangenome Analysis', 
                        fontsize=16, fontweight='bold', pad=20)
            ax.grid(True, alpha=0.3)
            ax.legend(fontsize=12)
            
            # Add summary statistics
            ax.text(0.02, 0.98, 
                   f'Strains: {len(plot_data)}\nCorrelation: r={correlation:.3f}\nMean Coverage: {plot_data["core_coverage_percent"].mean():.1f}%',
                   transform=ax.transAxes, fontsize=11, verticalalignment='top',
                   bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8))
            
            plt.tight_layout()
            plt.savefig(output_dir / 'core_coverage_vs_protein_count_blue_scatter.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    # 3. NEW: Static Version with Color-Coded CIRCULAR Legend 
    if 'protein_count' in core_results.columns:
        plot_data = core_results.dropna(subset=['protein_count'])
        
        if len(plot_data) > 0:
            # Calculate outliers
            z = np.polyfit(plot_data['protein_count'], plot_data['core_coverage_percent'], 1)
            p = np.poly1d(z)
            predicted_coverage = p(plot_data['protein_count'])
            residuals = plot_data['core_coverage_percent'] - predicted_coverage
            residual_threshold = np.percentile(residuals, 15)
            
            # Get outliers
            outliers_mask = residuals < residual_threshold
            outliers = plot_data[outliers_mask].copy()
            normal_points = plot_data[~outliers_mask]
            
            # Create the plot with more space for legend
            fig, (ax_main, ax_legend) = plt.subplots(1, 2, figsize=(22, 10), 
                                                     gridspec_kw={'width_ratios': [2.5, 1]})
            
            # Plot normal points
            ax_main.scatter(normal_points['protein_count'], normal_points['core_coverage_percent'], 
                           alpha=0.6, s=60, color='steelblue', label='Normal strains', zorder=2)
            
            # Generate unique colors for each outlier
            colors = plt.cm.Set1(np.linspace(0, 1, len(outliers)))  # Use Set1 colormap for distinct colors
            
            # Plot each outlier with unique color
            legend_info = []
            for i, (idx, row) in enumerate(outliers.iterrows()):
                color = colors[i]
                ax_main.scatter(row['protein_count'], row['core_coverage_percent'], 
                               s=120, color=color, edgecolors='black', linewidth=2, 
                               zorder=3, alpha=0.9)
                
                # Store info for legend
                full_organism_name = row.get('organism', row['strain_id'])
                legend_info.append((color, full_organism_name, row['core_coverage_percent']))
            
            # Add trend line
            correlation = plot_data['protein_count'].corr(plot_data['core_coverage_percent'])
            ax_main.plot(plot_data['protein_count'], p(plot_data['protein_count']), 
                        "r--", alpha=0.8, linewidth=2, label='Trend line', zorder=1)
            
            # Main plot styling
            ax_main.set_xlabel('Total Protein Count', fontsize=14, fontweight='bold')
            ax_main.set_ylabel('Core Genome Coverage (%)', fontsize=14, fontweight='bold')
            ax_main.set_title('Core Coverage vs Protein Count in Streptomyces\n(Each outlier has unique color - see legend)', 
                             fontsize=16, fontweight='bold', pad=20)
            ax_main.grid(True, alpha=0.3)
            ax_main.legend(loc='lower right', fontsize=12)
            
            # Add correlation info
            ax_main.text(0.02, 0.98, f'Pearson r = {correlation:.3f}', transform=ax_main.transAxes,
                        bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8),
                        fontsize=12, fontweight='bold', va='top')
            
            # Create detailed legend on the right with PROPER COLORED CIRCLES
            ax_legend.set_xlim(0, 1)
            ax_legend.set_ylim(0, len(legend_info) + 1)
            ax_legend.axis('off')
            
            # Add legend title
            ax_legend.text(0.5, len(legend_info) + 0.3, 'OUTLIER SPECIES', 
                          ha='center', va='center', fontsize=14, fontweight='bold',
                          bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgray', alpha=0.8))
            
            # Add each outlier to legend with PROPER COLORED CIRCLES using scatter
            for i, (color, full_name, coverage) in enumerate(reversed(legend_info)):
                y_pos = len(legend_info) - i - 0.5
                
                # Use scatter instead of Circle for proper color display
                ax_legend.scatter(0.08, y_pos, s=300, color=color, edgecolors='black', 
                                 linewidth=2, zorder=10, alpha=0.9)
                
                # Add full species name (wrap if too long)
                if len(full_name) > 45:
                    # Split long names into two lines at a sensible point
                    words = full_name.split()
                    if len(words) > 3:
                        # Try to split at middle, but keep species name together if possible
                        if 'strain' in full_name.lower():
                            # Split before strain info
                            strain_idx = next((j for j, word in enumerate(words) if 'strain' in word.lower()), len(words)//2)
                            line1 = ' '.join(words[:strain_idx])
                            line2 = ' '.join(words[strain_idx:])
                        else:
                            mid_point = len(words) // 2
                            line1 = ' '.join(words[:mid_point])
                            line2 = ' '.join(words[mid_point:])
                        text = f"{line1}\n{line2}"
                    else:
                        text = full_name
                else:
                    text = full_name
                
                # Position text next to colored circle
                ax_legend.text(0.18, y_pos, f"{text}\n({coverage:.1f}%)", 
                              ha='left', va='center', fontsize=9, 
                              bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='gray'))
            
            plt.tight_layout()
            plt.savefig(output_dir / 'color_coded_full_species_names.png', 
                       dpi=300, bbox_inches='tight', facecolor='white')
            plt.close()
            
            print(f"✅ Color-coded static plot saved with {len(outliers)} outliers and COLORED CIRCLE legend")

def create_simple_pangenome_composition(freq_analysis, output_dir):
    """Create simple pangenome composition - ONLY pie chart with readable text."""
    
    fig, ax = plt.subplots(1, 1, figsize=(10, 10))
    
    # Simplified categories: Unique, Core, Other Shared
    categories = ['Unique\n(1 strain only)', 'Core\n(ML-defined)', 'Other Shared\n(2+ strains, not core)']
    
    # FIXED: Use the correct keys from the freq_analysis
    sizes = [
        freq_analysis['unique_genes'], 
        freq_analysis.get('core_genes', freq_analysis.get('core_genes_ml_defined', 0)), 
        freq_analysis['other_shared_genes']  # FIXED: This key now exists!
    ]
    
    # IMPROVED: More distinguishable and eye-friendly colors
    colors = ['#DDA0DD', '#FF6347', '#20B2AA']
    
    # Create pie chart with better text formatting
    wedges, texts, autotexts = ax.pie(sizes, labels=categories, autopct='%1.1f%%', 
                                     colors=colors, startangle=90, 
                                     textprops={'fontsize': 14, 'fontweight': 'bold'},
                                     pctdistance=0.85)
    
    # Make percentage text larger and bold
    for autotext in autotexts:
        autotext.set_color('white')
        autotext.set_fontweight('bold')
        autotext.set_fontsize(16)
    
    # Make category labels more readable
    for text in texts:
        text.set_fontsize(13)
        text.set_fontweight('bold')
    
    ax.set_title('Streptomyces Pangenome Composition\nSimplified Gene Categories', 
                fontsize=18, fontweight='bold', pad=30)
    
    # Add gene counts as text annotations with improved colors
    total = sum(sizes)
    for i, (size, color) in enumerate(zip(sizes, colors)):
        percentage = (size / total) * 100
        ax.text(0.02, 0.98 - i*0.08, f'{categories[i].split()[0]}: {size:,} genes ({percentage:.1f}%)', 
               transform=ax.transAxes, fontsize=12, fontweight='bold',
               bbox=dict(boxstyle='round,pad=0.3', facecolor=color, alpha=0.8, edgecolor='darkgray'))
    
    plt.tight_layout()
    plt.savefig(output_dir / 'pangenome_composition_simplified.png', dpi=300, bbox_inches='tight')
    plt.close()


def create_all_strains_missing_plot(missing_counts, metadata, output_dir, max_strains=50):
    """Create the ALL strains missing core genes plot - with flexible display options.
    
    Parameters:
    -----------
    max_strains : int
        Maximum number of strains to display in plot. 
        0 = show all strains (might be very wide)
        >0 = show only top N strains with most missing genes
    """
    
    # Filter to strains with missing genes
    strains_with_missing = {k: v for k, v in missing_counts.items() if v > 0}
    
    if len(strains_with_missing) == 0:
        print("No strains have missing core genes")
        return
    
    # Sort by missing count (descending)
    sorted_strains = sorted(strains_with_missing.items(), key=lambda x: x[1], reverse=True)
    
    # ALWAYS save complete data to TSV first (before any filtering)
    print(f"Saving complete missing genes data for {len(sorted_strains)} strains to TSV...")
    complete_data = []
    for strain_id, count in sorted_strains:
        full_name = get_full_organism_name(strain_id, metadata)
        complete_data.append({
            'strain_id': strain_id,
            'organism': full_name,
            'missing_core_genes': count
        })
    
    df_complete = pd.DataFrame(complete_data)
    df_complete.to_csv(output_dir / 'all_strains_missing_core_genes_complete.tsv', 
                      sep='\t', index=False)
    print(f"Complete data saved: {len(sorted_strains)} strains in TSV")
    
    # Determine how many strains to plot
    if max_strains > 0 and len(sorted_strains) > max_strains:
        plot_strains = sorted_strains[:max_strains]
        plot_title = f'Top {max_strains} Strains with Most Missing Core Genes\n(out of {len(sorted_strains)} total strains with missing genes)'
        print(f"Limiting plot to top {max_strains} strains for readability")
    else:
        plot_strains = sorted_strains
        plot_title = f'All {len(plot_strains)} Streptomyces Strains with Missing Core Genes'
        if len(plot_strains) > 100:
            print(f"WARNING: Plotting {len(plot_strains)} strains - plot may be very wide!")
    
    # Get full organism names for plotting
    strain_names = []
    missing_values = []
    for strain_id, count in plot_strains:
        full_name = get_full_organism_name(strain_id, metadata)
        # Truncate very long names for plot readability
        if len(full_name) > 60:
            full_name = full_name[:57] + "..."
        strain_names.append(full_name)
        missing_values.append(count)
    
    # ---- Horizontal bar layout for readability ----
    n = len(plot_strains)

    # Dynamic figure HEIGHT based on number of strains (portrait)
    fig_height = max(8, min(0.35 * n, 20))   # e.g. ~17.5 for 50 strains
    fig, ax = plt.subplots(1, 1, figsize=(10, fig_height))

    # Color coding by severity (same logic)
    colors = []
    for count in missing_values:
        if count >= 50:
            colors.append('#d32f2f')  # Dark red
        elif count >= 20:
            colors.append('#f57c00')  # Dark orange  
        elif count >= 14:
            colors.append('#fbc02d')  # Yellow
        else:
            colors.append('#1976d2')  # Blue

    # Horizontal bars
    y_positions = range(len(strain_names))
    bars = ax.barh(y_positions, missing_values,
                   color=colors, alpha=0.8,
                   edgecolor='black', linewidth=0.5)

    # Species names on y-axis, italic + larger
    ax.set_yticks(y_positions)
    ax.set_yticklabels(strain_names, fontsize=12, fontstyle='italic')  # 🔼 bigger

    # Put the worst offenders at the TOP
    ax.invert_yaxis()

    ax.set_xlabel('Number of Missing Core Genes', fontsize=14, fontweight='bold')
    ax.set_title(plot_title, fontsize=16, fontweight='bold', pad=20)
    ax.grid(True, alpha=0.3, axis='x')

    # Threshold line at 14 missing genes (now vertical)
    ax.axvline(x=14, color='purple', linestyle='--', linewidth=2,
               alpha=0.7, label='14-gene threshold')

    # Add value labels to the right of bars (slightly bigger)
    max_val = max(missing_values)
    for i, (bar, count) in enumerate(zip(bars, missing_values)):
        show_label = (n <= 30 or
                      i % (2 if n <= 60 else 3) == 0 or
                      count >= 20)
        if show_label:
            width = bar.get_width()
            ax.text(width + max_val * 0.01,
                    bar.get_y() + bar.get_height() / 2,
                    f'{count}',
                    ha='left', va='center',
                    fontsize=10, fontweight='bold')  # 🔼 from 8 -> 10
    
    # Add summary box – moved to bottom-right
    summary_text = (f"Total strains with missing genes: {len(strains_with_missing)}\n"
                   f"Shown in plot: {len(plot_strains)}\n"
                   f"Mean missing: {np.mean(list(strains_with_missing.values())):.1f}\n"
                   f"Max missing: {max(strains_with_missing.values())}")
    
    ax.text(0.98, 0.02, summary_text, transform=ax.transAxes,
            bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.8),
            fontsize=11, verticalalignment='bottom',  # ⬅ bottom
            horizontalalignment='right', fontweight='bold')
    
    ax.legend(loc='upper center')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'all_strains_missing_core_genes.png', dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Missing core genes plot saved:")
    print(f"   - Plot shows: {len(plot_strains)} strains")
    print(f"   - Full data: {len(sorted_strains)} strains in TSV file")
    
    # Create a second summary plot if we're limiting the display
    if max_strains > 0 and len(sorted_strains) > max_strains:
        fig2, ax2 = plt.subplots(1, 1, figsize=(10, 6))
        
        bins = [0, 5, 10, 14, 20, 50, 100, max(strains_with_missing.values())+1]
        labels = ['1-5', '6-10', '11-14', '15-20', '21-50', '51-100', '100+']
        
        counts_per_bin = []
        for i in range(len(bins)-1):
            count = sum(1 for v in strains_with_missing.values() 
                       if bins[i] < v <= bins[i+1])
            counts_per_bin.append(count)
        
        used_labels = []
        used_counts = []
        for label, count in zip(labels, counts_per_bin):
            if count > 0:
                used_labels.append(label)
                used_counts.append(count)
        
        bars2 = ax2.bar(range(len(used_labels)), used_counts, color='coral', 
                       edgecolor='darkred', linewidth=1.5)
        ax2.set_xticks(range(len(used_labels)))
        ax2.set_xticklabels(used_labels)
        ax2.set_xlabel('Number of Missing Core Genes', fontsize=12, fontweight='bold')
        ax2.set_ylabel('Number of Strains', fontsize=12, fontweight='bold')
        ax2.set_title(f'Distribution of Missing Core Genes Across All {len(strains_with_missing)} Strains',
                     fontsize=14, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        
        for bar, count in zip(bars2, used_counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{count}', ha='center', va='bottom', fontsize=10, fontweight='bold')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'missing_core_genes_distribution.png', dpi=300, bbox_inches='tight')
        plt.close()
        print(f"Additional distribution plot created for all {len(strains_with_missing)} strains")

def create_enhanced_missing_genes_detailed_plot(core_gene_absence_frequency, annotations, output_dir):
    """Create enhanced detailed missing core genes plot with ACTUAL gene names and protein functions."""
    
    most_frequent_missing = core_gene_absence_frequency.most_common(25)
    
    if len(most_frequent_missing) == 0:
        return
    
    # Get actual gene names and protein functions if annotations are available
    gene_display_names = []
    frequencies = []
    
    for gene_cluster, frequency in most_frequent_missing:
        frequencies.append(frequency)
        
        # Try to get actual gene name and protein function from annotations
        display_name = gene_cluster  # Default fallback
        
        if annotations is not None:
            gene_info = annotations[annotations['gene'] == gene_cluster]
            if not gene_info.empty:
                row = gene_info.iloc[0]
                
                # Get protein product name from 'product' column
                if pd.notna(row.get('product', '')):
                    protein_product = str(row['product']).strip()
                    if protein_product and protein_product != 'nan' and protein_product.lower() != 'deleted':
                        if len(protein_product) > 60:
                            display_name = protein_product[:57] + "..."
                        else:
                            display_name = protein_product
                    else:
                        display_name = gene_cluster
                else:
                    display_name = gene_cluster
        
        gene_display_names.append(display_name)
    
    # Create enhanced horizontal plot
    fig, ax = plt.subplots(1, 1, figsize=(16, max(14, len(gene_display_names) * 0.8)))
    
    # Create gradient colors (professional blue-red scale)
    max_freq = max(frequencies)
    colors = []
    for freq in frequencies:
        intensity = freq / max_freq
        if intensity >= 0.8:
            colors.append('#1565C0')  # Dark blue - most frequent
        elif intensity >= 0.6:
            colors.append('#1976D2')  # Medium blue
        elif intensity >= 0.4:
            colors.append('#42A5F5')  # Light blue
        elif intensity >= 0.2:
            colors.append('#81C784')  # Light green
        else:
            colors.append('#A5D6A7')  # Very light green - least frequent
    
    # Create horizontal bars
    bars = ax.barh(range(len(gene_display_names)), frequencies, color=colors, alpha=0.9, 
                   edgecolor='darkgray', linewidth=1.5)
    
    # Set labels and formatting
    ax.set_yticks(range(len(gene_display_names)))
    ax.set_yticklabels(gene_display_names, fontsize=10, fontweight='bold')
    ax.set_xlabel('Number of Strains Missing This Gene', fontsize=16, fontweight='bold')
    ax.set_title('Most Frequently Missing Core Genes\n(With Gene Names and Protein Functions)', 
                fontsize=16, fontweight='bold', pad=25)
    ax.grid(True, alpha=0.3, axis='x', color='gray', linestyle='-', linewidth=0.5)
    
    # Add frequency labels with styled boxes
    for i, (bar, freq) in enumerate(zip(bars, frequencies)):
        width = bar.get_width()
        ax.text(width + max(frequencies)*0.01, bar.get_y() + bar.get_height()/2,
                f'{freq}', ha='left', va='center', fontsize=11, fontweight='bold',
                bbox=dict(boxstyle='round,pad=0.3', facecolor='white', 
                         edgecolor='darkblue', alpha=0.9, linewidth=1))
    
    # Improve plot appearance
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_linewidth(2)
    ax.spines['bottom'].set_linewidth(2)
    ax.set_facecolor('#f8f9fa')
    fig.patch.set_facecolor('white')
    
    # Adjust layout for better gene name visibility
    plt.tight_layout()
    plt.subplots_adjust(left=0.5)  # More space for gene names and protein functions
    plt.savefig(output_dir / 'most_frequently_missing_core_genes_detailed.png', 
               dpi=300, bbox_inches='tight', facecolor='white')
    plt.close()
    print(f"✅ Enhanced detailed missing genes plot saved with gene names and protein functions")

# Removed old complex visualization functions - keeping only essential plots as requested

def save_structural_results(core_results, freq_analysis, unique_df, missing_counts, 
                          core_gene_absence_frequency, df_genes, output_dir, annotations=None):
    """Save all structural analysis results."""
    print("\n=== Saving Structural Analysis Results ===")
    
    output_dir = Path(output_dir)
    
    # 1. Core genome structural analysis
    core_results_sorted = core_results.sort_values('core_coverage_percent', ascending=False)
    core_results_sorted.to_csv(output_dir / 'core_genome_structural_analysis.tsv', sep='\t', index=False)
    
    # 2. Unique genes per strain
    unique_df_sorted = unique_df.sort_values('unique_genes', ascending=False)
    unique_df_sorted.to_csv(output_dir / 'unique_genes_per_strain_analysis.tsv', sep='\t', index=False)
    
    # 3. Missing genes structural data
    missing_data = []
    for strain_id, count in missing_counts.items():
        if count > 0:
            # Get organism name from core_results DataFrame directly
            organism_row = core_results[core_results['strain_id'] == strain_id]
            if not organism_row.empty and 'organism' in organism_row.columns:
                organism = organism_row.iloc[0]['organism']
            else:
                organism = strain_id
            missing_data.append({
                'strain_id': strain_id,
                'organism': organism,
                'missing_core_genes_count': count
            })
    
    if missing_data:
        df_missing = pd.DataFrame(missing_data)
        df_missing = df_missing.sort_values('missing_core_genes_count', ascending=False)
        df_missing.to_csv(output_dir / 'missing_core_genes_structural.tsv', sep='\t', index=False)
    
    # 4. Gene frequency analysis
    freq_data = []
    gene_counts = freq_analysis['gene_counts']
    gene_names = df_genes.index
    
    for gene, count in zip(gene_names, gene_counts):
        freq_data.append({
            'gene': gene,
            'strain_count': int(count),
            'frequency_percent': (int(count) / len(core_results)) * 100
        })
    
    df_freq = pd.DataFrame(freq_data)
    df_freq = df_freq.sort_values('strain_count', ascending=False)
    df_freq.to_csv(output_dir / 'gene_frequency_analysis.tsv', sep='\t', index=False)
    
    # 5. Frequently missing core genes
    if core_gene_absence_frequency:
        missing_freq_data = []
        for gene, absence_count in core_gene_absence_frequency.most_common():
            row = {
                'gene': gene,
                'product': '',
                'strains_missing': absence_count,
                'missing_percentage': (absence_count / len(core_results)) * 100
            }
            
            # Add product name if available
            if annotations is not None:
                gene_info = annotations[annotations['gene'] == gene]
                if not gene_info.empty:
                    if pd.notna(gene_info.iloc[0].get('product', '')):
                        product = str(gene_info.iloc[0]['product']).strip()
                        if product and product != 'nan' and product.lower() != 'deleted':
                            row['product'] = product
            
            missing_freq_data.append(row)
        
        df_missing_freq = pd.DataFrame(missing_freq_data)
        df_missing_freq.to_csv(output_dir / 'frequently_missing_core_genes_structural.tsv', sep='\t', index=False)
    # 6. Summary report
    summary_lines = [
        "GENE-LEVEL STRUCTURAL ANALYSIS SUMMARY",
        "=" * 45,
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "DATASET OVERVIEW:",
        f"  Total strains analyzed: {len(core_results)}",
        f"  Total genes in pangenome: {freq_analysis['total_genes']:,}",
        f"  Core genes (ML-defined): {freq_analysis['core_genes_ml_defined']:,}",
        "",
        "CORE GENOME STRUCTURE:",
        f"  Mean core coverage: {core_results['core_coverage_percent'].mean():.2f}%",
        f"  Coverage range: {core_results['core_coverage_percent'].min():.1f}% - {core_results['core_coverage_percent'].max():.1f}%",
        f"  Standard deviation: {core_results['core_coverage_percent'].std():.2f}%",
        f"  Strains with <90% coverage: {len(core_results[core_results['core_coverage_percent'] < 90])}",
        "",
        "GENE FREQUENCY DISTRIBUTION:",
        f"  Unique genes (1 strain): {freq_analysis['unique_genes']:,} ({freq_analysis['unique_genes']/freq_analysis['total_genes']*100:.1f}%)",
        f"  Rare genes (2-5% strains): {freq_analysis['rare_genes']:,} ({freq_analysis['rare_genes']/freq_analysis['total_genes']*100:.1f}%)",
        f"  Accessory genes (5-95%): {freq_analysis['accessory_genes']:,} ({freq_analysis['accessory_genes']/freq_analysis['total_genes']*100:.1f}%)",
        f"  Core genes (≥95% strains): {freq_analysis['core_genes_by_frequency']:,} ({freq_analysis['core_genes_by_frequency']/freq_analysis['total_genes']*100:.1f}%)",
        "",
        "UNIQUE GENES:",
        f"  Mean per strain: {unique_df['unique_genes'].mean():.0f} genes",
        f"  Range: {unique_df['unique_genes'].min()}-{unique_df['unique_genes'].max()} genes",
        f"  Total unique genes: {freq_analysis['unique_genes']:,}",
        "",
        "MISSING CORE GENES:",
        f"  Strains with missing core genes: {len([c for c in missing_counts.values() if c > 0])}",
        f"  Mean missing per strain: {np.mean(list(missing_counts.values())):.1f}",
        f"  Maximum missing: {max(missing_counts.values()) if missing_counts else 0}",
        "",
        "FILES GENERATED:",
        "  • core_genome_structural_analysis.tsv - Core coverage per strain",
        "  • unique_genes_per_strain_analysis.tsv - Unique genes analysis",
        "  • gene_frequency_analysis.tsv - Frequency of all genes",
        "  • missing_core_genes_structural.tsv - Missing core genes per strain",
        "  • frequently_missing_core_genes_structural.tsv - Most frequently missing genes",
        "  • Multiple visualization PNG files",
        "  • interactive_core_coverage_analysis.html - Interactive plot (if available)"
    ]
    
    with open(output_dir / 'structural_analysis_summary.txt', 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"Structural analysis results saved to: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Gene-Level Structural Analysis')
    parser.add_argument('--gene-matrix', required=True, help='Path to gene matrix NPZ file')
    parser.add_argument('--gene-labels', required=True, help='Path to gene labels file')
    parser.add_argument('--core-genes', required=True, help='Path to core genes list file')
    parser.add_argument('--metadata', help='Path to proteome metadata TSV file')
    parser.add_argument('--annotations', help='Path to core gene annotations file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--max-strains-plot', type=int, default=50, 
                       help='Maximum strains to show in missing genes plot (0=all, default=50)')
    
    args = parser.parse_args()
    

    print("🧬 GENE-LEVEL STRUCTURAL ANALYSIS")
    print("=" * 40)
    
    try:
        # Load data
        df_genes, core_genes, metadata, annotations = load_data(
            args.gene_matrix, args.gene_labels, args.core_genes, args.metadata, args.annotations
        )
        
        print(f"Analyzing {df_genes.shape[0]} genes across {df_genes.shape[1]} strains")
        
        # Run structural analyses
        core_results, structural_outliers = analyze_core_genome_structure(df_genes, core_genes, metadata)
        freq_analysis, unique_genes, rare_genes, accessory_genes_clean = analyze_gene_frequency_distribution(df_genes, core_genes)
        unique_df = analyze_unique_genes_per_strain(df_genes, unique_genes, metadata)
        missing_genes_per_strain, missing_counts, core_gene_absence_frequency = identify_missing_core_genes_structural(df_genes, core_genes)
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate essential visualizations only
        create_structural_visualizations(
            core_results, freq_analysis, unique_df, missing_counts, 
            core_gene_absence_frequency, metadata, annotations, output_dir,
            args.max_strains_plot
        )
        
        # Save results
        save_structural_results(
            core_results, freq_analysis, unique_df, missing_counts, 
            core_gene_absence_frequency, df_genes, output_dir, annotations
        )
        
        print(f"\nStructural analysis completed!")
        print(f"Results saved to {output_dir}")
        
        # Print key findings
        print(f"\nKey structural findings:")
        print(f"  • Mean core coverage: {core_results['core_coverage_percent'].mean():.1f}%")
        print(f"  • {freq_analysis['unique_genes']:,} unique genes ({freq_analysis['unique_genes']/freq_analysis['total_genes']*100:.1f}% of pangenome)")
        print(f"  • {freq_analysis['core_genes']:,} core genes ({freq_analysis['core_genes']/freq_analysis['total_genes']*100:.1f}% of pangenome)")
        print(f"  • {freq_analysis['other_shared_genes']:,} other shared genes ({freq_analysis['other_shared_genes']/freq_analysis['total_genes']*100:.1f}% of pangenome)")
        print(f"  • {len([c for c in missing_counts.values() if c > 0])} strains have missing core genes")
        print(f"  • {len(structural_outliers)} structural outliers identified")
        
        if 'protein_count' in core_results.columns and len(core_results.dropna(subset=['protein_count'])) > 0:
            plot_data = core_results.dropna(subset=['protein_count'])
            correlation = plot_data['protein_count'].corr(plot_data['core_coverage_percent'])
            print(f"  • Protein count vs core coverage correlation: r={correlation:.3f}")
        
        # Print summary of generated files
        print(f"\n📊 Generated Files:")
        print(f"  • Core coverage plots (interactive HTML + static PNG + color-coded)")
        print(f"  • Pangenome composition pie chart (improved colors)")
        print(f"  • All strains missing genes (vertical bars)")
        print(f"  • Missing genes detailed plot (with gene names + protein functions)")
        print(f"  • Missing genes summary table (TSV + HTML)")
        print(f"  • Comprehensive TSV data files for further analysis")
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
