#!/usr/bin/env python3
import os
import sys
import argparse
import pandas as pd
import sparse_utils

def load_core_genes(core_genes_file):
    """Load list of core genes from file."""
    core_genes = set()
    with open(core_genes_file, 'r') as f:
        for line in f:
            gene = line.strip()
            if gene:
                core_genes.add(gene)
    return core_genes

def check_core_coverage(gene_matrix, core_genes):
    """Calculate core genome coverage for each strain."""
    print(f"Analyzing {len(core_genes)} core genes across {len(gene_matrix.columns)} strains...")
    
    # Get indices of core genes in the matrix
    core_gene_indices = []
    missing_core_genes = []
    
    for gene in core_genes:
        if gene in gene_matrix.index_map:
            core_gene_indices.append(gene_matrix.index_map[gene])
        else:
            missing_core_genes.append(gene)
    
    if missing_core_genes:
        print(f"Warning: {len(missing_core_genes)} core genes not found in matrix")
    
    # Calculate coverage per strain
    coverage_data = []
    for i, strain in enumerate(gene_matrix.columns):
        # Count present core genes (non-zero values)
        present_genes = 0
        for core_idx in core_gene_indices:
            if gene_matrix.data.tocsc()[core_idx, i] > 0:
                present_genes += 1
        
        total_core = len(core_gene_indices)
        coverage_pct = (present_genes / total_core) * 100 if total_core > 0 else 0
        
        coverage_data.append({
            'strain': strain,
            'core_genes_present': present_genes,
            'total_core_genes': total_core,
            'coverage_percentage': round(coverage_pct, 2)
        })
    
    coverage_df = pd.DataFrame(coverage_data)
    coverage_df = coverage_df.sort_values('coverage_percentage', ascending=False)
    
    return coverage_df

def main():
    parser = argparse.ArgumentParser(description='Check core genome coverage for each strain')
    parser.add_argument('--gene-matrix', required=True, help='Path to strain_by_gene.npz file')
    parser.add_argument('--core-genes', required=True, help='Path to core genes list file')
    parser.add_argument('--output', required=True, help='Output TSV file')
    parser.add_argument('--threshold', type=float, default=90.0, help='Coverage threshold for flagging (default: 90.0)')
    
    args = parser.parse_args()
    
    print("=== CORE GENOME COVERAGE CHECK ===")
    print(f"Gene matrix: {args.gene_matrix}")
    print(f"Core genes: {args.core_genes}")
    print(f"Output: {args.output}")
    print(f"Threshold: {args.threshold}%")
    
    # Load data
    print("\nLoading data...")
    gene_matrix = sparse_utils.read_lsdf(args.gene_matrix)
    core_genes = load_core_genes(args.core_genes)
    
    print(f"Loaded: {gene_matrix.shape[0]} genes × {gene_matrix.shape[1]} strains")
    print(f"Core genes: {len(core_genes)}")
    
    # Calculate coverage
    coverage_df = check_core_coverage(gene_matrix, core_genes)
    
    # Save results
    coverage_df.to_csv(args.output, sep='\t', index=False)
    print(f"\nResults saved to: {args.output}")
    
    # Summary statistics
    print(f"\n=== SUMMARY ===")
    print(f"Mean coverage: {coverage_df['coverage_percentage'].mean():.2f}%")
    print(f"Median coverage: {coverage_df['coverage_percentage'].median():.2f}%")
    print(f"Min coverage: {coverage_df['coverage_percentage'].min():.2f}%")
    print(f"Max coverage: {coverage_df['coverage_percentage'].max():.2f}%")
    
    # Problematic strains
    problematic = coverage_df[coverage_df['coverage_percentage'] < args.threshold]
    if len(problematic) > 0:
        print(f"\nStrains with <{args.threshold}% core coverage ({len(problematic)} strains):")
        for _, row in problematic.iterrows():
            print(f"  {row['strain']}: {row['coverage_percentage']}% ({row['core_genes_present']}/{row['total_core_genes']})")
    else:
        print(f"\nAll strains have ≥{args.threshold}% core coverage! ✓")
    
    # Top and bottom performers
    print(f"\nTop 5 strains (highest coverage):")
    for _, row in coverage_df.head(5).iterrows():
        print(f"  {row['strain']}: {row['coverage_percentage']}%")
    
    print(f"\nBottom 5 strains (lowest coverage):")
    for _, row in coverage_df.tail(5).iterrows():
        print(f"  {row['strain']}: {row['coverage_percentage']}%")

if __name__ == "__main__":
    main()
