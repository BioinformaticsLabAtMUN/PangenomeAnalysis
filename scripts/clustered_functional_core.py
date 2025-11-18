#!/usr/bin/env python3
"""
GOATools Clustering-Based Functional Core Genome Analysis
Focuses on functional annotation analysis using GOATools GO term clustering.
Handles: GOATools cluster analysis, missing functional clusters, strain-wise patterns.

Adapted from Revigo analysis to work with GOATools clustering outputs.
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
import re
warnings.filterwarnings('ignore')

# Import custom sparse utilities
try:
    import sparse_utils
except ImportError:
    print("ERROR: sparse_utils module not found. Please ensure it's in your PYTHONPATH.")
    sys.exit(1)

def load_goatools_data(goatools_bp_file, goatools_mf_file):
    """Load and parse GOATools clustering results."""
    print("Loading GOATools clustering data...")
    
    # Load GOATools BP data
    goatools_bp = pd.read_csv(goatools_bp_file, sep='\t')
    print(f"Loaded {len(goatools_bp)} GO Biological Process terms from GOATools")
    
    # Load GOATools MF data
    goatools_mf = pd.read_csv(goatools_mf_file, sep='\t')
    print(f"Loaded {len(goatools_mf)} GO Molecular Function terms from GOATools")
    
    # Clean the data - remove quotes if present and handle whitespace
    def clean_column(col):
        if col.dtype == 'object':
            return col.astype(str).str.replace('"', '').str.strip()
        return col
    
    # GOATools format: Cluster, TermID, Name, Representative_Name, Value
    goatools_bp['TermID'] = clean_column(goatools_bp['TermID'])
    goatools_bp['Name'] = clean_column(goatools_bp['Name'])
    goatools_bp['Representative_Name'] = clean_column(goatools_bp['Representative_Name'])
    
    goatools_mf['TermID'] = clean_column(goatools_mf['TermID'])
    goatools_mf['Name'] = clean_column(goatools_mf['Name'])
    goatools_mf['Representative_Name'] = clean_column(goatools_mf['Representative_Name'])
    
    # Create term name to representative name mappings
    bp_term_to_rep = dict(zip(goatools_bp['Name'], goatools_bp['Representative_Name']))
    mf_term_to_rep = dict(zip(goatools_mf['Name'], goatools_mf['Representative_Name']))
    
    print("Processing BP representatives...")
    
    # Process BP data - build representatives structure
    bp_representatives = {}
    for _, row in goatools_bp.iterrows():
        go_term = row['Name']
        rep_name = row['Representative_Name']
        
        print(f"  {go_term} -> {rep_name}")
        
        # Build representative clusters
        if rep_name not in bp_representatives:
            bp_representatives[rep_name] = {
                'name': rep_name,
                'members': [go_term],
                'total_value': row['Value']
            }
        else:
            bp_representatives[rep_name]['members'].append(go_term)
            bp_representatives[rep_name]['total_value'] += row['Value']
    
    print("Processing MF representatives...")
    
    # Process MF data - build representatives structure
    mf_representatives = {}
    for _, row in goatools_mf.iterrows():
        go_term = row['Name']
        rep_name = row['Representative_Name']
        
        print(f"  {go_term} -> {rep_name}")
        
        # Build representative clusters
        if rep_name not in mf_representatives:
            mf_representatives[rep_name] = {
                'name': rep_name,
                'members': [go_term],
                'total_value': row['Value']
            }
        else:
            mf_representatives[rep_name]['members'].append(go_term)
            mf_representatives[rep_name]['total_value'] += row['Value']
    
    print(f"Created {len(bp_representatives)} BP clusters and {len(mf_representatives)} MF clusters")
    
    # Print some examples for debugging
    print("\nExample BP clusters:")
    for i, (rep_name, cluster_info) in enumerate(list(bp_representatives.items())[:3]):
        print(f"  Cluster {i+1}: {rep_name}")
        print(f"    Members: {cluster_info['members'][:3]}{'...' if len(cluster_info['members']) > 3 else ''}")
        print(f"    Total value: {cluster_info['total_value']}")
    
    return {
        'bp_term_to_rep': bp_term_to_rep,
        'mf_term_to_rep': mf_term_to_rep,
        'bp_representatives': bp_representatives,
        'mf_representatives': mf_representatives,
        'missing_representatives_bp': [],  # GOATools doesn't have missing reps like Revigo
        'missing_representatives_mf': []
    }

def load_allele_cluster_mapping(allele_names_file):
    """Load pre-computed gene to cluster mapping from allele names file.
    
    Expected format:
    Original_ID    New_ID    Cluster_ID    Allele_Number    representative
    WP_329167035.1    Strep_C0A0    Strep_C0    0    True
    """
    print("Loading pre-computed cluster mappings...")
    
    allele_names = pd.read_csv(allele_names_file, sep='\t')
    print(f"Loaded {len(allele_names)} gene-to-cluster mappings")
    
    # Create mapping from Original_ID to Cluster_ID
    gene_to_cluster = dict(zip(allele_names['Original_ID'], allele_names['Cluster_ID']))
    
    # Also create mapping from New_ID to Cluster_ID (in case genes are referenced by new names)
    new_id_to_cluster = dict(zip(allele_names['New_ID'], allele_names['Cluster_ID']))
    gene_to_cluster.update(new_id_to_cluster)
    
    print(f"Created mapping for {len(gene_to_cluster)} unique gene identifiers")
    print(f"Example mapping: {list(gene_to_cluster.items())[:3]}")
    
    return gene_to_cluster

def load_data(gene_matrix, gene_labels, core_genes_file, annotations_file, allele_names_file, metadata_file=None):
    """Load all required data files."""
    print("Loading data files...")
    
    # Load gene matrix
    df_genes = sparse_utils.read_lsdf(gene_matrix, gene_labels)
    print(f"Loaded gene matrix: {df_genes.shape}")
    
    # Load core genes
    with open(core_genes_file, 'r') as f:
        core_genes = [line.strip() for line in f if line.strip()]
    print(f"Loaded {len(core_genes)} core genes")
    
    # Load annotations
    annotations = pd.read_csv(annotations_file, sep='\t')
    print(f"Loaded annotations for {len(annotations)} genes")
    
    # Load pre-computed cluster mappings
    gene_to_cluster = load_allele_cluster_mapping(allele_names_file)
    
    # Debug: Print column information
    print(f"Annotation columns: {list(annotations.columns)}")
    
    # Check which GO columns are available
    go_cols_found = []
    if 'Gene Ontology (biological process)' in annotations.columns:
        go_cols_found.append('Gene Ontology (biological process)')
    if 'Gene Ontology (molecular function)' in annotations.columns:
        go_cols_found.append('Gene Ontology (molecular function)')
    if 'go_process' in annotations.columns:
        go_cols_found.append('go_process')
    if 'go_function' in annotations.columns:
        go_cols_found.append('go_function')
    if 'go_component' in annotations.columns:
        go_cols_found.append('go_component')
        
    print(f"GO-related columns found: {go_cols_found}")
    
    # Load metadata if provided
    metadata = None
    if metadata_file and Path(metadata_file).exists():
        try:
            metadata = pd.read_csv(metadata_file, sep='\t')
            print(f"Loaded metadata for {len(metadata)} proteomes")
        except Exception as e:
            print(f"Warning: Could not load metadata: {e}")
    
    return df_genes, core_genes, annotations, gene_to_cluster, metadata

def get_full_organism_name(proteome_id, metadata):
    """Get complete organism name from metadata."""
    if metadata is None:
        return proteome_id
    
    matching_row = metadata[metadata['Proteome Id'] == proteome_id]
    if not matching_row.empty:
        organism_name = matching_row.iloc[0]['Organism']
        return organism_name
    
    return proteome_id

def parse_go_terms(go_terms_str, debug_count=[0]):
    """Parse GO terms from string format with improved error handling.
    
    Handles two formats:
    1. Legacy format: term [GO:1234567];term2 [GO:2345678]
    2. New annotateSmartDominantAlleles format: term_name|GO_ID||evidence_code|term2|GO_ID2||evidence_code
    """
    if pd.isna(go_terms_str) or str(go_terms_str).strip() == '' or str(go_terms_str) == 'nan':
        return []
    
    go_terms = []
    go_str = str(go_terms_str).strip()
    
    # Check for new pipe-separated format first
    if '|' in go_str:
        # New format: term_name|GO_ID||evidence_code separated by pipes
        # Split by pipes and extract every 4th element (term names)
        parts = go_str.split('|')
        for i in range(0, len(parts), 4):
            if i < len(parts):
                term = parts[i].strip()
                if term and term != '':
                    go_terms.append(term)
        # Debug: Show parsing result for first 5 entries only
        if len(go_terms) > 0 and debug_count[0] < 5:
            print(f"    Parsed pipe format: {go_str[:100]}... -> {go_terms[:3]}{'...' if len(go_terms) > 3 else ''}")
            debug_count[0] += 1
    # Handle legacy semicolon-separated format
    elif ';' in go_str:
        terms = go_str.split(';')
        for term in terms:
            term = term.strip()
            if '[GO:' in term:
                desc = term.split('[GO:')[0].strip()
                if desc and desc != '':
                    go_terms.append(desc)
            elif term and term != '':
                # Handle cases where GO terms don't have [GO:] format
                go_terms.append(term)
    else:
        if '[GO:' in go_str:
            desc = go_str.split('[GO:')[0].strip()
            if desc and desc != '':
                go_terms.append(desc)
        elif go_str and go_str != '':
            # Handle single term without [GO:] format
            go_terms.append(go_str)
    
    return go_terms

def map_genes_to_goatools_clusters(annotations, goatools_data):
    """Map genes to GOATools functional clusters with debugging."""
    print("Mapping genes to GOATools functional clusters...")
    
    gene_to_bp_clusters = {}
    gene_to_mf_clusters = {}
    unmapped_bp_terms = set()
    unmapped_mf_terms = set()
    
    bp_term_to_rep = goatools_data['bp_term_to_rep']
    mf_term_to_rep = goatools_data['mf_term_to_rep']
    
    for _, row in annotations.iterrows():
        gene_id = row['gene']
        
        # Map BP terms to clusters - handle both old and new column names
        bp_col = None
        if 'Gene Ontology (biological process)' in row.index:
            bp_col = 'Gene Ontology (biological process)'
        elif 'go_process' in row.index:
            bp_col = 'go_process'
            
        bp_terms = parse_go_terms(row.get(bp_col, '')) if bp_col else []
        bp_clusters = set()
        for term in bp_terms:
            if term in bp_term_to_rep:
                bp_clusters.add(bp_term_to_rep[term])
            else:
                unmapped_bp_terms.add(term)
        
        if bp_clusters:
            gene_to_bp_clusters[gene_id] = list(bp_clusters)
        
        # Map MF terms to clusters - handle both old and new column names
        mf_col = None
        if 'Gene Ontology (molecular function)' in row.index:
            mf_col = 'Gene Ontology (molecular function)'
        elif 'go_function' in row.index:
            mf_col = 'go_function'
            
        mf_terms = parse_go_terms(row.get(mf_col, '')) if mf_col else []
        mf_clusters = set()
        for term in mf_terms:
            if term in mf_term_to_rep:
                mf_clusters.add(mf_term_to_rep[term])
            else:
                unmapped_mf_terms.add(term)
        
        if mf_clusters:
            gene_to_mf_clusters[gene_id] = list(mf_clusters)
    
    print(f"Mapped {len(gene_to_bp_clusters)} genes to BP clusters")
    print(f"Mapped {len(gene_to_mf_clusters)} genes to MF clusters")
    
    # DEBUGGING: Report unmapped terms
    if unmapped_bp_terms:
        print(f"\n⚠️  WARNING: {len(unmapped_bp_terms)} BP terms not found in GOATools mapping:")
        for term in list(unmapped_bp_terms)[:5]:  # Show first 5
            print(f"    '{term}'")
        if len(unmapped_bp_terms) > 5:
            print(f"    ... and {len(unmapped_bp_terms) - 5} more")
    
    if unmapped_mf_terms:
        print(f"\n⚠️  WARNING: {len(unmapped_mf_terms)} MF terms not found in GOATools mapping:")
        for term in list(unmapped_mf_terms)[:5]:  # Show first 5
            print(f"    '{term}'")
        if len(unmapped_mf_terms) > 5:
            print(f"    ... and {len(unmapped_mf_terms) - 5} more")
    
    return gene_to_bp_clusters, gene_to_mf_clusters

def analyze_core_goatools_clusters(core_genes, gene_to_bp_clusters, gene_to_mf_clusters, goatools_data):
    """Analyze GOATools functional clusters in core genes with debugging."""
    print("\n=== Core Genes GOATools Cluster Analysis ===")
    
    # Count core genes in each BP cluster
    core_bp_clusters = Counter()
    core_bp_genes = {}
    core_genes_with_bp = 0
    
    for gene in core_genes:
        if gene in gene_to_bp_clusters:
            core_genes_with_bp += 1
            for cluster in gene_to_bp_clusters[gene]:
                core_bp_clusters[cluster] += 1
                if cluster not in core_bp_genes:
                    core_bp_genes[cluster] = []
                core_bp_genes[cluster].append(gene)
    
    # Count core genes in each MF cluster
    core_mf_clusters = Counter()
    core_mf_genes = {}
    core_genes_with_mf = 0
    
    for gene in core_genes:
        if gene in gene_to_mf_clusters:
            core_genes_with_mf += 1
            for cluster in gene_to_mf_clusters[gene]:
                core_mf_clusters[cluster] += 1
                if cluster not in core_mf_genes:
                    core_mf_genes[cluster] = []
                core_mf_genes[cluster].append(gene)
    
    print(f"Core genes with BP annotations: {core_genes_with_bp}/{len(core_genes)} ({core_genes_with_bp/len(core_genes)*100:.1f}%)")
    print(f"Core genes with MF annotations: {core_genes_with_mf}/{len(core_genes)} ({core_genes_with_mf/len(core_genes)*100:.1f}%)")
    
    # Create summary DataFrames
    bp_cluster_summary = []
    for cluster, count in core_bp_clusters.most_common():
        cluster_info = goatools_data['bp_representatives'][cluster]
        bp_cluster_summary.append({
            'cluster': cluster,
            'gene_count': count,
            'total_cluster_value': cluster_info['total_value'],
            'member_terms': '; '.join(cluster_info['members'][:5]),  # First 5 terms
            'core_genes': '; '.join(core_bp_genes[cluster][:5])  # First 5 genes
        })
    
    mf_cluster_summary = []
    for cluster, count in core_mf_clusters.most_common():
        cluster_info = goatools_data['mf_representatives'][cluster]
        mf_cluster_summary.append({
            'cluster': cluster,
            'gene_count': count,
            'total_cluster_value': cluster_info['total_value'],
            'member_terms': '; '.join(cluster_info['members'][:5]),
            'core_genes': '; '.join(core_mf_genes[cluster][:5])
        })
    
    print(f"Core genes span {len(core_bp_clusters)} BP clusters and {len(core_mf_clusters)} MF clusters")
    
    return {
        'bp_cluster_summary': bp_cluster_summary,
        'mf_cluster_summary': mf_cluster_summary,
        'core_bp_genes': core_bp_genes,
        'core_mf_genes': core_mf_genes,
        'stats': {
            'core_genes_with_bp': core_genes_with_bp,
            'core_genes_with_mf': core_genes_with_mf,
            'total_core_genes': len(core_genes)
        }
    }

def analyze_missing_goatools_clusters(df_genes, core_genes, gene_to_bp_clusters, gene_to_mf_clusters, 
                                    core_analysis, metadata):
    """Analyze missing GOATools functional clusters per strain."""
    print("\n=== Missing GOATools Clusters Analysis ===")
    
    core_genes_in_matrix = [g for g in core_genes if g in df_genes.index]
    core_df = df_genes.labelslice(indices=core_genes_in_matrix)
    
    core_bp_genes = core_analysis['core_bp_genes']
    core_mf_genes = core_analysis['core_mf_genes']
    
    strain_missing_analysis = {}
    missing_bp_clusters = Counter()
    missing_mf_clusters = Counter()
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        strain_column = core_df.data.tocsc().getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(core_genes_in_matrix[i] for i in present_gene_indices)
        missing_genes = set(core_genes_in_matrix) - present_genes
        
        if len(missing_genes) == 0:
            continue
        
        # Find missing BP clusters
        missing_bp_clusters_strain = set()
        for cluster, cluster_genes in core_bp_genes.items():
            cluster_genes_in_matrix = [g for g in cluster_genes if g in core_genes_in_matrix]
            if cluster_genes_in_matrix:
                missing_from_cluster = [g for g in cluster_genes_in_matrix if g in missing_genes]
                if missing_from_cluster:
                    missing_bp_clusters_strain.add(cluster)
                    missing_bp_clusters[cluster] += 1
        
        # Find missing MF clusters
        missing_mf_clusters_strain = set()
        for cluster, cluster_genes in core_mf_genes.items():
            cluster_genes_in_matrix = [g for g in cluster_genes if g in core_genes_in_matrix]
            if cluster_genes_in_matrix:
                missing_from_cluster = [g for g in cluster_genes_in_matrix if g in missing_genes]
                if missing_from_cluster:
                    missing_mf_clusters_strain.add(cluster)
                    missing_mf_clusters[cluster] += 1
        
        organism = get_full_organism_name(strain_id, metadata)
        
        strain_missing_analysis[strain_id] = {
            'organism': organism,
            'missing_gene_count': len(missing_genes),
            'missing_bp_clusters': missing_bp_clusters_strain,
            'missing_mf_clusters': missing_mf_clusters_strain
        }
    
    print(f"Strains with missing functional clusters: {len(strain_missing_analysis)}")
    print(f"Most frequently missing BP clusters: {len(missing_bp_clusters)}")
    print(f"Most frequently missing MF clusters: {len(missing_mf_clusters)}")
    
    return {
        'strain_analysis': strain_missing_analysis,
        'missing_bp_clusters': missing_bp_clusters,
        'missing_mf_clusters': missing_mf_clusters
    }

def analyze_completely_absent_clusters_optimized(df_genes, gene_to_bp_clusters, gene_to_mf_clusters, 
                                               core_analysis, metadata, gene_to_cluster):
    """OPTIMIZED: Analyze strains that completely lack any genes from functional clusters.
    
    Uses pre-computed cluster mappings to avoid expensive nested loops.
    """
    print("\n=== Completely Absent Functional Clusters Analysis (OPTIMIZED) ===")
    
    core_bp_genes = core_analysis['core_bp_genes']
    core_mf_genes = core_analysis['core_mf_genes']
    
    # Pre-compute which genes belong to each functional cluster
    print("Pre-computing functional cluster memberships...")
    bp_cluster_to_genes = {}
    mf_cluster_to_genes = {}
    
    # Build reverse mapping from functional clusters to genes
    for gene, bp_clusters in gene_to_bp_clusters.items():
        for cluster in bp_clusters:
            if cluster not in bp_cluster_to_genes:
                bp_cluster_to_genes[cluster] = set()
            bp_cluster_to_genes[cluster].add(gene)
    
    for gene, mf_clusters in gene_to_mf_clusters.items():
        for cluster in mf_clusters:
            if cluster not in mf_cluster_to_genes:
                mf_cluster_to_genes[cluster] = set()
            mf_cluster_to_genes[cluster].add(gene)
    
    print(f"Built mappings for {len(bp_cluster_to_genes)} BP clusters and {len(mf_cluster_to_genes)} MF clusters")
    
    strain_absent_analysis = {}
    absent_bp_clusters = Counter()
    absent_mf_clusters = Counter()
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        # Get ALL genes present in this strain
        strain_column = df_genes.data.tocsc().getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(df_genes.index[i] for i in present_gene_indices)
        
        # Check which BP clusters are completely absent
        absent_bp_clusters_strain = set()
        for cluster in core_bp_genes.keys():
            if cluster in bp_cluster_to_genes:
                # Check if ANY genes from this functional cluster are present
                cluster_genes = bp_cluster_to_genes[cluster]
                genes_in_strain = cluster_genes.intersection(present_genes)
                
                if len(genes_in_strain) == 0:
                    absent_bp_clusters_strain.add(cluster)
                    absent_bp_clusters[cluster] += 1
        
        # Check which MF clusters are completely absent  
        absent_mf_clusters_strain = set()
        for cluster in core_mf_genes.keys():
            if cluster in mf_cluster_to_genes:
                cluster_genes = mf_cluster_to_genes[cluster]
                genes_in_strain = cluster_genes.intersection(present_genes)
                
                if len(genes_in_strain) == 0:
                    absent_mf_clusters_strain.add(cluster)
                    absent_mf_clusters[cluster] += 1
        
        organism = get_full_organism_name(strain_id, metadata)
        
        # Only record strains that have some completely absent clusters
        if absent_bp_clusters_strain or absent_mf_clusters_strain:
            strain_absent_analysis[strain_id] = {
                'organism': organism,
                'absent_bp_clusters': absent_bp_clusters_strain,
                'absent_mf_clusters': absent_mf_clusters_strain
            }
    
    print(f"Strains with completely absent functional clusters: {len(strain_absent_analysis)}")
    print(f"BP clusters completely absent in some strains: {len(absent_bp_clusters)}")
    print(f"MF clusters completely absent in some strains: {len(absent_mf_clusters)}")
    
    return {
        'strain_analysis': strain_absent_analysis,
        'absent_bp_clusters': absent_bp_clusters,
        'absent_mf_clusters': absent_mf_clusters
    }

# Keep the old function as backup
def analyze_completely_absent_clusters(df_genes, gene_to_bp_clusters, gene_to_mf_clusters, 
                                      core_analysis, metadata):
    """Analyze strains that completely lack any genes from functional clusters."""
    print("\n=== Completely Absent Functional Clusters Analysis ===")
    
    core_bp_genes = core_analysis['core_bp_genes']
    core_mf_genes = core_analysis['core_mf_genes']
    
    strain_absent_analysis = {}
    absent_bp_clusters = Counter()
    absent_mf_clusters = Counter()
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        # Get ALL genes present in this strain (not just core genes)
        strain_column = df_genes.data.tocsc().getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(df_genes.index[i] for i in present_gene_indices)
        
        # Check which BP clusters are completely absent
        absent_bp_clusters_strain = set()
        for cluster in core_bp_genes.keys():
            # Find ALL genes in the genome that belong to this cluster (not just core genes)
            cluster_genes_in_strain = []
            for gene in present_genes:
                if gene in gene_to_bp_clusters and cluster in gene_to_bp_clusters[gene]:
                    cluster_genes_in_strain.append(gene)
            
            # If no genes found for this cluster, it's completely absent
            if len(cluster_genes_in_strain) == 0:
                absent_bp_clusters_strain.add(cluster)
                absent_bp_clusters[cluster] += 1
        
        # Check which MF clusters are completely absent
        absent_mf_clusters_strain = set()
        for cluster in core_mf_genes.keys():
            # Find ALL genes in the genome that belong to this cluster (not just core genes)
            cluster_genes_in_strain = []
            for gene in present_genes:
                if gene in gene_to_mf_clusters and cluster in gene_to_mf_clusters[gene]:
                    cluster_genes_in_strain.append(gene)
            
            # If no genes found for this cluster, it's completely absent
            if len(cluster_genes_in_strain) == 0:
                absent_mf_clusters_strain.add(cluster)
                absent_mf_clusters[cluster] += 1
        
        organism = get_full_organism_name(strain_id, metadata)
        
        # Only record strains that have some completely absent clusters
        if absent_bp_clusters_strain or absent_mf_clusters_strain:
            strain_absent_analysis[strain_id] = {
                'organism': organism,
                'absent_bp_clusters': absent_bp_clusters_strain,
                'absent_mf_clusters': absent_mf_clusters_strain
            }
    
    print(f"Strains with completely absent functional clusters: {len(strain_absent_analysis)}")
    print(f"BP clusters completely absent in some strains: {len(absent_bp_clusters)}")
    print(f"MF clusters completely absent in some strains: {len(absent_mf_clusters)}")
    
    return {
        'strain_analysis': strain_absent_analysis,
        'absent_bp_clusters': absent_bp_clusters,
        'absent_mf_clusters': absent_mf_clusters
    }

def create_goatools_visualizations(core_analysis, missing_analysis, absent_analysis, output_dir):
    """Create GOATools-based visualizations."""
    print("\n=== Creating GOATools Visualizations ===")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # 1. Core functional clusters
    create_core_cluster_plots(core_analysis, output_dir)
    
    # 2. Missing functional clusters (incomplete core clusters)
    create_missing_cluster_plots(missing_analysis, output_dir)
    
    # 3. Completely absent functional clusters
    create_absent_cluster_plots(absent_analysis, output_dir)
    
    # 4. Strain-wise incomplete patterns
    create_strain_missing_plots(missing_analysis, output_dir)
    
    # 5. Strain-wise completely absent patterns
    create_strain_absent_plots(absent_analysis, output_dir)

def create_core_cluster_plots(core_analysis, output_dir):
    """Create core functional cluster plots."""
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 1. Top BP clusters in core genes
    bp_summary = core_analysis['bp_cluster_summary'][:15]  # Top 15
    if bp_summary:
        clusters = [item['cluster'][:50] + '...' if len(item['cluster']) > 50 else item['cluster'] 
                   for item in bp_summary]
        counts = [item['gene_count'] for item in bp_summary]
        
        bars1 = ax1.barh(range(len(clusters)), counts, color='lightblue', alpha=0.8)
        ax1.set_yticks(range(len(clusters)))
        ax1.set_yticklabels(clusters, fontsize=9)
        ax1.set_xlabel('Number of Core Genes')
        ax1.set_title('Core Genes: Top GOATools BP Clusters', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars1, counts):
            width = bar.get_width()
            ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    # 2. Top MF clusters in core genes
    mf_summary = core_analysis['mf_cluster_summary'][:15]  # Top 15
    if mf_summary:
        clusters = [item['cluster'][:50] + '...' if len(item['cluster']) > 50 else item['cluster'] 
                   for item in mf_summary]
        counts = [item['gene_count'] for item in mf_summary]
        
        bars2 = ax2.barh(range(len(clusters)), counts, color='lightgreen', alpha=0.8)
        ax2.set_yticks(range(len(clusters)))
        ax2.set_yticklabels(clusters, fontsize=9)
        ax2.set_xlabel('Number of Core Genes')
        ax2.set_title('Core Genes: Top GOATools MF Clusters', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars2, counts):
            width = bar.get_width()
            ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'clustered_core_functional_clusters.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_missing_cluster_plots(missing_analysis, output_dir):
    """Create missing functional cluster plots."""
    
    missing_bp = missing_analysis['missing_bp_clusters'].most_common(15)
    missing_mf = missing_analysis['missing_mf_clusters'].most_common(15)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 1. Most frequently missing BP clusters
    if missing_bp:
        clusters = [item[0][:45] + '...' if len(item[0]) > 45 else item[0] for item in missing_bp]
        counts = [item[1] for item in missing_bp]
        
        bars1 = ax1.barh(range(len(clusters)), counts, color='lightcoral', alpha=0.8)
        ax1.set_yticks(range(len(clusters)))
        ax1.set_yticklabels(clusters, fontsize=9)
        ax1.set_xlabel('Number of Strains')
        ax1.set_title('Incomplete Core BP Clusters\n(Strains missing some core genes from cluster)', 
                     fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars1, counts):
            width = bar.get_width()
            ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    # 2. Most frequently missing MF clusters
    if missing_mf:
        clusters = [item[0][:45] + '...' if len(item[0]) > 45 else item[0] for item in missing_mf]
        counts = [item[1] for item in missing_mf]
        
        bars2 = ax2.barh(range(len(clusters)), counts, color='lightsteelblue', alpha=0.8)
        ax2.set_yticks(range(len(clusters)))
        ax2.set_yticklabels(clusters, fontsize=9)
        ax2.set_xlabel('Number of Strains')
        ax2.set_title('Incomplete Core MF Clusters\n(Strains missing some core genes from cluster)', 
                     fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars2, counts):
            width = bar.get_width()
            ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'clustered_incomplete_functional_clusters.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_absent_cluster_plots(absent_analysis, output_dir):
    """Create completely absent functional cluster plots."""
    
    absent_bp = absent_analysis['absent_bp_clusters'].most_common(15)
    absent_mf = absent_analysis['absent_mf_clusters'].most_common(15)
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 8))
    
    # 1. Most frequently completely absent BP clusters
    if absent_bp:
        clusters = [item[0][:45] + '...' if len(item[0]) > 45 else item[0] for item in absent_bp]
        counts = [item[1] for item in absent_bp]
        
        bars1 = ax1.barh(range(len(clusters)), counts, color='darkred', alpha=0.8)
        ax1.set_yticks(range(len(clusters)))
        ax1.set_yticklabels(clusters, fontsize=9)
        ax1.set_xlabel('Number of Strains')
        ax1.set_title('Completely Absent BP Clusters\n(Strains with ZERO genes of this function)', 
                     fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars1, counts):
            width = bar.get_width()
            ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    else:
        ax1.text(0.5, 0.5, 'No completely absent\nBP clusters found', 
                ha='center', va='center', transform=ax1.transAxes, fontsize=14)
        ax1.set_title('Completely Absent BP Clusters\n(Strains with ZERO genes of this function)', 
                     fontsize=12, fontweight='bold')
    
    # 2. Most frequently completely absent MF clusters
    if absent_mf:
        clusters = [item[0][:45] + '...' if len(item[0]) > 45 else item[0] for item in absent_mf]
        counts = [item[1] for item in absent_mf]
        
        bars2 = ax2.barh(range(len(clusters)), counts, color='darkblue', alpha=0.8)
        ax2.set_yticks(range(len(clusters)))
        ax2.set_yticklabels(clusters, fontsize=9)
        ax2.set_xlabel('Number of Strains')
        ax2.set_title('Completely Absent MF Clusters\n(Strains with ZERO genes of this function)', 
                     fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars2, counts):
            width = bar.get_width()
            ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    else:
        ax2.text(0.5, 0.5, 'No completely absent\nMF clusters found', 
                ha='center', va='center', transform=ax2.transAxes, fontsize=14)
        ax2.set_title('Completely Absent MF Clusters\n(Strains with ZERO genes of this function)', 
                     fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'clustered_completely_absent_clusters.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_strain_missing_plots(missing_analysis, output_dir):
    """Create strain-wise incomplete functional patterns plot."""
    
    strain_analysis = missing_analysis['strain_analysis']
    
    # Get strains with most incomplete clusters
    strain_missing_data = []
    for strain_id, data in strain_analysis.items():
        total_missing_clusters = len(data['missing_bp_clusters']) + len(data['missing_mf_clusters'])
        strain_missing_data.append({
            'strain_id': strain_id,
            'organism': data['organism'],
            'missing_genes': data['missing_gene_count'],
            'missing_bp_clusters': len(data['missing_bp_clusters']),
            'missing_mf_clusters': len(data['missing_mf_clusters']),
            'total_missing_clusters': total_missing_clusters
        })
    
    # Sort by total missing clusters and take top 15
    strain_missing_data.sort(key=lambda x: x['total_missing_clusters'], reverse=True)
    top_strains = strain_missing_data[:15]
    
    if top_strains:
        fig, ax = plt.subplots(1, 1, figsize=(14, 8))
        
        # Missing clusters per strain
        organism_names = []
        bp_counts = []
        mf_counts = []
        
        for strain_data in top_strains:
            organism = strain_data['organism']
            if len(organism) > 40:
                organism_names.append(organism[:37] + '...')
            else:
                organism_names.append(organism)
            bp_counts.append(strain_data['missing_bp_clusters'])
            mf_counts.append(strain_data['missing_mf_clusters'])
        
        x_pos = np.arange(len(organism_names))
        width = 0.35
        
        bars1 = ax.bar(x_pos - width/2, bp_counts, width, label='Incomplete BP Clusters', 
                       color='lightcoral', alpha=0.7)
        bars2 = ax.bar(x_pos + width/2, mf_counts, width, label='Incomplete MF Clusters', 
                       color='lightsteelblue', alpha=0.7)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(organism_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Number of Incomplete Functional Clusters')
        ax.set_title('Strains with Most Incomplete Functional Clusters\n(Missing some core genes from clusters)', 
                    fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'clustered_strain_incomplete_patterns.png', dpi=300, bbox_inches='tight')
        plt.close()

def create_strain_absent_plots(absent_analysis, output_dir):
    """Create strain-wise completely absent functional patterns plot."""
    
    strain_analysis = absent_analysis['strain_analysis']
    
    if not strain_analysis:
        # Create empty plot with message
        fig, ax = plt.subplots(1, 1, figsize=(14, 8))
        ax.text(0.5, 0.5, 'No strains with completely absent\nfunctional clusters found', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        ax.set_title('Strains with Completely Absent Functional Clusters\n(Zero genes of these functions)', 
                    fontsize=12, fontweight='bold')
        ax.set_xlabel('Organisms')
        ax.set_ylabel('Number of Completely Absent Functional Clusters')
        plt.tight_layout()
        plt.savefig(output_dir / 'clustered_strain_completely_absent_patterns.png', dpi=300, bbox_inches='tight')
        plt.close()
        return
    
    # Get strains with most completely absent clusters
    strain_absent_data = []
    for strain_id, data in strain_analysis.items():
        total_absent_clusters = len(data['absent_bp_clusters']) + len(data['absent_mf_clusters'])
        strain_absent_data.append({
            'strain_id': strain_id,
            'organism': data['organism'],
            'absent_bp_clusters': len(data['absent_bp_clusters']),
            'absent_mf_clusters': len(data['absent_mf_clusters']),
            'total_absent_clusters': total_absent_clusters
        })
    
    # Sort by total absent clusters and take top 15
    strain_absent_data.sort(key=lambda x: x['total_absent_clusters'], reverse=True)
    top_strains = strain_absent_data[:15]
    
    if top_strains:
        fig, ax = plt.subplots(1, 1, figsize=(14, 8))
        
        # Absent clusters per strain
        organism_names = []
        bp_counts = []
        mf_counts = []
        
        for strain_data in top_strains:
            organism = strain_data['organism']
            if len(organism) > 40:
                organism_names.append(organism[:37] + '...')
            else:
                organism_names.append(organism)
            bp_counts.append(strain_data['absent_bp_clusters'])
            mf_counts.append(strain_data['absent_mf_clusters'])
        
        x_pos = np.arange(len(organism_names))
        width = 0.35
        
        bars1 = ax.bar(x_pos - width/2, bp_counts, width, label='Completely Absent BP Clusters', 
                       color='darkred', alpha=0.7)
        bars2 = ax.bar(x_pos + width/2, mf_counts, width, label='Completely Absent MF Clusters', 
                       color='darkblue', alpha=0.7)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(organism_names, rotation=45, ha='right', fontsize=8)
        ax.set_ylabel('Number of Completely Absent Functional Clusters')
        ax.set_title('Strains with Completely Absent Functional Clusters\n(Zero genes of these functions)', 
                    fontsize=12, fontweight='bold')
        ax.legend()
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        plt.savefig(output_dir / 'clustered_strain_completely_absent_patterns.png', dpi=300, bbox_inches='tight')
        plt.close()

def create_absent_clusters_html(strain_absent_df, absent_analysis, output_dir):
    """Create an HTML version of the completely absent clusters analysis with improved readability."""
    
    output_file = output_dir / 'strain_completely_absent_clustered_clusters.html'
    
    # Calculate summary statistics
    total_strains = len(strain_absent_df)
    strains_with_absent_bp = len(strain_absent_df[strain_absent_df['completely_absent_bp_clusters_count'] > 0]) if total_strains > 0 else 0
    strains_with_absent_mf = len(strain_absent_df[strain_absent_df['completely_absent_mf_clusters_count'] > 0]) if total_strains > 0 else 0
    max_absent_bp = strain_absent_df['completely_absent_bp_clusters_count'].max() if total_strains > 0 else 0
    max_absent_mf = strain_absent_df['completely_absent_mf_clusters_count'].max() if total_strains > 0 else 0
    
    # HTML template with embedded CSS and JavaScript (same structure as Revigo version)
    html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Completely Absent Functional Clusters Analysis (GOATools)</title>
    <style>
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            margin: 0;
            padding: 20px;
            background-color: #f5f5f5;
            line-height: 1.6;
        }}
        .container {{
            max-width: 1400px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            text-align: center;
            margin-bottom: 30px;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        .summary {{
            background-color: #ecf0f1;
            padding: 20px;
            border-radius: 8px;
            margin-bottom: 30px;
            border-left: 5px solid #3498db;
        }}
        .summary h2 {{
            color: #2c3e50;
            margin-top: 0;
        }}
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
            margin-top: 15px;
        }}
        .summary-item {{
            background-color: white;
            padding: 15px;
            border-radius: 6px;
            text-align: center;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        .summary-number {{
            font-size: 24px;
            font-weight: bold;
            color: #e74c3c;
        }}
        .summary-label {{
            font-size: 14px;
            color: #7f8c8d;
            margin-top: 5px;
        }}
        .explanation {{
            background-color: #e3f2fd;
            padding: 15px;
            border-radius: 6px;
            margin-bottom: 20px;
            border-left: 4px solid #2196f3;
        }}
        .explanation h3 {{
            margin-top: 0;
            color: #1976d2;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background-color: white;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        th {{
            background-color: #34495e;
            color: white;
            padding: 15px 10px;
            text-align: left;
            font-weight: 600;
        }}
        td {{
            padding: 12px 10px;
            border-bottom: 1px solid #ecf0f1;
            vertical-align: top;
        }}
        tr:nth-child(even) {{
            background-color: #f8f9fa;
        }}
        tr:hover {{
            background-color: #e3f2fd;
        }}
        .strain-id {{
            font-family: 'Courier New', monospace;
            font-weight: bold;
            color: #2c3e50;
        }}
        .organism-name {{
            max-width: 300px;
            word-wrap: break-word;
        }}
        .count-cell {{
            text-align: center;
            font-weight: bold;
            padding: 8px;
            border-radius: 4px;
        }}
        .count-high {{
            background-color: #ffebee;
            color: #c62828;
        }}
        .count-medium {{
            background-color: #fff3e0;
            color: #ef6c00;
        }}
        .count-low {{
            background-color: #e8f5e8;
            color: #2e7d32;
        }}
        .count-zero {{
            background-color: #f5f5f5;
            color: #757575;
        }}
        .function-list {{
            max-width: 400px;
            word-wrap: break-word;
        }}
        .function-item {{
            display: inline-block;
            background-color: #ffcdd2;
            color: #b71c1c;
            padding: 3px 8px;
            margin: 2px;
            border-radius: 15px;
            font-size: 12px;
            border: 1px solid #ef5350;
        }}
        .no-functions {{
            color: #757575;
            font-style: italic;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>🧬 Completely Absent Functional Clusters Analysis (GOATools)</h1>
        
        <div class="explanation">
            <h3>📋 What This Analysis Shows</h3>
            <p><strong>Completely Absent Clusters:</strong> Strains that have <em>zero genes</em> for entire functional categories anywhere in their genome. This analysis uses GOATools clustering instead of Revigo.</p>
            <p><strong>BP (Biological Process):</strong> What the cell does (e.g., DNA repair, metabolism)</p>
            <p><strong>MF (Molecular Function):</strong> What proteins do (e.g., enzyme activity, binding)</p>
        </div>
        
        <div class="summary">
            <h2>📊 Summary Statistics</h2>
            <div class="summary-grid">
                <div class="summary-item">
                    <div class="summary-number">{total_strains}</div>
                    <div class="summary-label">Total Strains Analyzed</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{strains_with_absent_bp}</div>
                    <div class="summary-label">Strains Missing BP Functions</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{strains_with_absent_mf}</div>
                    <div class="summary-label">Strains Missing MF Functions</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{max_absent_bp}</div>
                    <div class="summary-label">Max BP Functions Lost</div>
                </div>
                <div class="summary-item">
                    <div class="summary-number">{max_absent_mf}</div>
                    <div class="summary-label">Max MF Functions Lost</div>
                </div>
            </div>
        </div>
        
        <table>
            <thead>
                <tr>
                    <th>Strain ID</th>
                    <th>Organism</th>
                    <th>Absent BP Count</th>
                    <th>Absent MF Count</th>
                    <th>Absent Biological Process Functions</th>
                    <th>Absent Molecular Function Functions</th>
                </tr>
            </thead>
            <tbody>
"""
    
    # Add table rows
    if total_strains == 0:
        html_content += """
                <tr>
                    <td colspan="6" style="text-align: center; padding: 40px; color: #757575; font-style: italic;">
                        🎉 No strains with completely absent functional clusters found!<br>
                        This suggests all strains retain at least some genes for all core functional categories.
                    </td>
                </tr>
"""
    else:
        for _, row in strain_absent_df.iterrows():
            bp_count = row['completely_absent_bp_clusters_count']
            mf_count = row['completely_absent_mf_clusters_count']
            
            # Color coding for counts
            def get_count_class(count):
                if count == 0:
                    return 'count-zero'
                elif count <= 2:
                    return 'count-low'
                elif count <= 5:
                    return 'count-medium'
                else:
                    return 'count-high'
            
            bp_class = get_count_class(bp_count)
            mf_class = get_count_class(mf_count)
            
            # Format function lists
            bp_functions = row['completely_absent_bp_clusters'].split('; ') if row['completely_absent_bp_clusters'] else []
            mf_functions = row['completely_absent_mf_clusters'].split('; ') if row['completely_absent_mf_clusters'] else []
            
            bp_functions_html = ''
            if bp_functions and bp_functions[0]:
                bp_functions_html = ' '.join([f'<span class="function-item">{func}</span>' for func in bp_functions])
            else:
                bp_functions_html = '<span class="no-functions">None</span>'
            
            mf_functions_html = ''
            if mf_functions and mf_functions[0]:
                mf_functions_html = ' '.join([f'<span class="function-item">{func}</span>' for func in mf_functions])
            else:
                mf_functions_html = '<span class="no-functions">None</span>'
            
            html_content += f"""
                <tr>
                    <td class="strain-id">{row['strain_id']}</td>
                    <td class="organism-name">{row['organism']}</td>
                    <td class="count-cell {bp_class}">{bp_count}</td>
                    <td class="count-cell {mf_class}">{mf_count}</td>
                    <td class="function-list">{bp_functions_html}</td>
                    <td class="function-list">{mf_functions_html}</td>
                </tr>
"""
    
    # Close HTML
    html_content += """
            </tbody>
        </table>
    </div>
</body>
</html>
"""
    
    # Write HTML file
    with open(output_file, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"Created interactive HTML report: {output_file}")

def save_goatools_results(core_analysis, missing_analysis, absent_analysis, goatools_data, output_dir):
    """Save all GOATools analysis results."""
    print("\n=== Saving GOATools Analysis Results ===")
    
    output_dir = Path(output_dir)
    
    # 1. Core BP clusters
    bp_df = pd.DataFrame(core_analysis['bp_cluster_summary'])
    bp_df.to_csv(output_dir / 'clustered_core_bp_clusters.tsv', sep='\t', index=False)
    
    # 2. Core MF clusters
    mf_df = pd.DataFrame(core_analysis['mf_cluster_summary'])
    mf_df.to_csv(output_dir / 'clustered_core_mf_clusters.tsv', sep='\t', index=False)
    
    # 3. Missing BP analysis (incomplete core clusters)
    missing_bp_data = []
    for cluster, count in missing_analysis['missing_bp_clusters'].most_common():
        missing_bp_data.append({
            'cluster': cluster,
            'strains_with_incomplete_cluster': count
        })
    
    missing_bp_df = pd.DataFrame(missing_bp_data)
    missing_bp_df.to_csv(output_dir / 'clustered_incomplete_bp_analysis.tsv', sep='\t', index=False)
    
    # 4. Missing MF analysis (incomplete core clusters)
    missing_mf_data = []
    for cluster, count in missing_analysis['missing_mf_clusters'].most_common():
        missing_mf_data.append({
            'cluster': cluster,
            'strains_with_incomplete_cluster': count
        })
    
    missing_mf_df = pd.DataFrame(missing_mf_data)
    missing_mf_df.to_csv(output_dir / 'clustered_incomplete_mf_analysis.tsv', sep='\t', index=False)
    
    # 5. Completely absent BP analysis
    absent_bp_data = []
    for cluster, count in absent_analysis['absent_bp_clusters'].most_common():
        absent_bp_data.append({
            'cluster': cluster,
            'strains_completely_lacking_function': count
        })
    
    absent_bp_df = pd.DataFrame(absent_bp_data)
    absent_bp_df.to_csv(output_dir / 'clustered_completely_absent_bp_analysis.tsv', sep='\t', index=False)
    
    # 6. Completely absent MF analysis
    absent_mf_data = []
    for cluster, count in absent_analysis['absent_mf_clusters'].most_common():
        absent_mf_data.append({
            'cluster': cluster,
            'strains_completely_lacking_function': count
        })
    
    absent_mf_df = pd.DataFrame(absent_mf_data)
    absent_mf_df.to_csv(output_dir / 'clustered_completely_absent_mf_analysis.tsv', sep='\t', index=False)
    
    # 7. Strain-wise missing clusters
    strain_data = []
    for strain_id, data in missing_analysis['strain_analysis'].items():
        strain_data.append({
            'strain_id': strain_id,
            'organism': data['organism'],
            'missing_gene_count': data['missing_gene_count'],
            'incomplete_bp_clusters_count': len(data['missing_bp_clusters']),
            'incomplete_mf_clusters_count': len(data['missing_mf_clusters']),
            'incomplete_bp_clusters': '; '.join(data['missing_bp_clusters']),
            'incomplete_mf_clusters': '; '.join(data['missing_mf_clusters'])
        })
    
    strain_df = pd.DataFrame(strain_data)
    strain_df = strain_df.sort_values('missing_gene_count', ascending=False)
    strain_df.to_csv(output_dir / 'strain_incomplete_clustered_clusters.tsv', sep='\t', index=False)
    
    # 8. Strain-wise completely absent clusters
    strain_absent_data = []
    for strain_id, data in absent_analysis['strain_analysis'].items():
        strain_absent_data.append({
            'strain_id': strain_id,
            'organism': data['organism'],
            'completely_absent_bp_clusters_count': len(data['absent_bp_clusters']),
            'completely_absent_mf_clusters_count': len(data['absent_mf_clusters']),
            'completely_absent_bp_clusters': '; '.join(data['absent_bp_clusters']),
            'completely_absent_mf_clusters': '; '.join(data['absent_mf_clusters'])
        })
    
    strain_absent_df = pd.DataFrame(strain_absent_data)
    strain_absent_df = strain_absent_df.sort_values('completely_absent_bp_clusters_count', ascending=False)
    strain_absent_df.to_csv(output_dir / 'strain_completely_absent_clustered_clusters.tsv', sep='\t', index=False)
    
    # 8b. Create HTML version of completely absent clusters
    create_absent_clusters_html(strain_absent_df, absent_analysis, output_dir)
    
    # 9. Summary report
    summary_lines = [
        "GOATOOLS CLUSTERING-BASED FUNCTIONAL CORE GENOME ANALYSIS SUMMARY",
        "=" * 65,
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "CORE FUNCTIONAL CLUSTERS:",
        f"  GOATools BP clusters represented: {len(core_analysis['bp_cluster_summary'])}",
        f"  GOATools MF clusters represented: {len(core_analysis['mf_cluster_summary'])}",
        "",
        "ANNOTATION COVERAGE:",
        f"  Core genes with BP annotations: {core_analysis['stats']['core_genes_with_bp']}/{core_analysis['stats']['total_core_genes']} ({core_analysis['stats']['core_genes_with_bp']/core_analysis['stats']['total_core_genes']*100:.1f}%)",
        f"  Core genes with MF annotations: {core_analysis['stats']['core_genes_with_mf']}/{core_analysis['stats']['total_core_genes']} ({core_analysis['stats']['core_genes_with_mf']/core_analysis['stats']['total_core_genes']*100:.1f}%)",
        "",
        "INCOMPLETE FUNCTIONAL CLUSTERS (missing some core genes):",
        f"  Strains with incomplete clusters: {len(missing_analysis['strain_analysis'])}",
        f"  BP clusters with missing core genes: {len(missing_analysis['missing_bp_clusters'])}",
        f"  MF clusters with missing core genes: {len(missing_analysis['missing_mf_clusters'])}",
        "",
        "COMPLETELY ABSENT FUNCTIONAL CLUSTERS (zero genes of function):",
        f"  Strains with completely absent clusters: {len(absent_analysis['strain_analysis'])}",
        f"  BP clusters completely absent in some strains: {len(absent_analysis['absent_bp_clusters'])}",
        f"  MF clusters completely absent in some strains: {len(absent_analysis['absent_mf_clusters'])}",
        "",
        "TOP CORE FUNCTIONAL CLUSTERS:",
    ]
    
    # Add top BP clusters
    for i, cluster_info in enumerate(core_analysis['bp_cluster_summary'][:5]):
        summary_lines.append(f"  {i+1}. {cluster_info['cluster']}: {cluster_info['gene_count']} genes")
    
    summary_lines.extend([
        "",
        "TOP INCOMPLETE FUNCTIONAL CLUSTERS:",
    ])
    
    # Add top incomplete clusters
    for i, (cluster, count) in enumerate(missing_analysis['missing_bp_clusters'].most_common(5)):
        summary_lines.append(f"  {i+1}. {cluster}: incomplete in {count} strains")
    
    if absent_analysis['absent_bp_clusters']:
        summary_lines.extend([
            "",
            "TOP COMPLETELY ABSENT FUNCTIONAL CLUSTERS:",
        ])
        
        # Add top completely absent clusters
        for i, (cluster, count) in enumerate(absent_analysis['absent_bp_clusters'].most_common(5)):
            summary_lines.append(f"  {i+1}. {cluster}: completely absent in {count} strains")
    
    summary_lines.extend([
        "",
        "FILES GENERATED:",
        "  • clustered_core_*_clusters.tsv - Core genes by cluster",
        "  • clustered_incomplete_*_analysis.tsv - Incomplete core clusters analysis",
        "  • clustered_completely_absent_*_analysis.tsv - Completely absent clusters analysis",
        "  • strain_incomplete_clustered_clusters.tsv - Per-strain incomplete clusters",
        "  • strain_completely_absent_clustered_clusters.tsv - Per-strain absent clusters (TSV)",
        "  • strain_completely_absent_clustered_clusters.html - Per-strain absent clusters (Interactive HTML)",
        "  • clustered_*.png - Functional cluster visualizations"
    ])
    
    with open(output_dir / 'clustered_functional_summary.txt', 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"GOATools analysis results saved to: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='GOATools Clustering-Based Functional Core Genome Analysis')
    parser.add_argument('--gene-matrix', required=True, help='Path to gene matrix NPZ file')
    parser.add_argument('--gene-labels', required=True, help='Path to gene labels file')
    parser.add_argument('--core-genes', required=True, help='Path to core genes list file')
    parser.add_argument('--annotations', required=True, help='Path to gene annotations TSV file')
    parser.add_argument('--metadata', help='Path to proteome metadata TSV file')
    parser.add_argument('--goatools-bp', required=True, help='Path to GOATools BP clustering TSV file')
    parser.add_argument('--goatools-mf', required=True, help='Path to GOATools MF clustering TSV file')
    parser.add_argument('--allele-names', required=True, help='Path to allele names TSV file for cluster mapping')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    parser.add_argument('--debug', action='store_true', help='Enable debugging output')
    
    args = parser.parse_args()
    
    print("🧬 GOATOOLS CLUSTERING-BASED FUNCTIONAL CORE GENOME ANALYSIS")
    print("=" * 65)
    
    try:
        # Load GOATools data
        goatools_data = load_goatools_data(args.goatools_bp, args.goatools_mf)
        
        # Load other data
        df_genes, core_genes, annotations, gene_to_cluster, metadata = load_data(
            args.gene_matrix, args.gene_labels, args.core_genes, 
            args.annotations, args.allele_names, args.metadata
        )
        
        print(f"Analyzing functional clusters in {df_genes.shape[0]} genes across {df_genes.shape[1]} strains")
        
        # Map genes to GOATools clusters
        gene_to_bp_clusters, gene_to_mf_clusters = map_genes_to_goatools_clusters(annotations, goatools_data)
        
        # Analyze core genes by GOATools clusters
        core_analysis = analyze_core_goatools_clusters(core_genes, gene_to_bp_clusters, gene_to_mf_clusters, goatools_data)
        
        # Analyze missing functional clusters (can also be optimized later if needed)
        missing_analysis = analyze_missing_goatools_clusters(
            df_genes, core_genes, gene_to_bp_clusters, gene_to_mf_clusters, core_analysis, metadata
        )
        
        # Analyze completely absent functional clusters (OPTIMIZED)
        absent_analysis = analyze_completely_absent_clusters_optimized(
            df_genes, gene_to_bp_clusters, gene_to_mf_clusters, core_analysis, metadata, gene_to_cluster
        )
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate visualizations
        create_goatools_visualizations(core_analysis, missing_analysis, absent_analysis, output_dir)
        
        # Save results
        save_goatools_results(core_analysis, missing_analysis, absent_analysis, goatools_data, output_dir)
        
        print(f"\nGOATools functional analysis completed!")
        print(f"Results saved to {output_dir}")
        
        # Print key findings
        print(f"\nKey GOATools functional findings:")
        print(f"  • Core BP clusters: {len(core_analysis['bp_cluster_summary'])}")
        print(f"  • Core MF clusters: {len(core_analysis['mf_cluster_summary'])}")
        print(f"  • Annotation coverage: BP {core_analysis['stats']['core_genes_with_bp']}/{core_analysis['stats']['total_core_genes']} ({core_analysis['stats']['core_genes_with_bp']/core_analysis['stats']['total_core_genes']*100:.1f}%), MF {core_analysis['stats']['core_genes_with_mf']}/{core_analysis['stats']['total_core_genes']} ({core_analysis['stats']['core_genes_with_mf']/core_analysis['stats']['total_core_genes']*100:.1f}%)")
        print(f"  • Strains with incomplete clusters: {len(missing_analysis['strain_analysis'])}")
        print(f"  • BP clusters with incomplete core genes: {len(missing_analysis['missing_bp_clusters'])}")
        print(f"  • MF clusters with incomplete core genes: {len(missing_analysis['missing_mf_clusters'])}")
        print(f"  • Strains with completely absent clusters: {len(absent_analysis['strain_analysis'])}")
        print(f"  • BP clusters completely absent in some strains: {len(absent_analysis['absent_bp_clusters'])}")
        print(f"  • MF clusters completely absent in some strains: {len(absent_analysis['absent_mf_clusters'])}")
        
        # Top findings
        if core_analysis['bp_cluster_summary']:
            top_bp = core_analysis['bp_cluster_summary'][0]
            print(f"  • Top core BP cluster: {top_bp['cluster']} ({top_bp['gene_count']} genes)")
        
        if missing_analysis['missing_bp_clusters']:
            top_incomplete = missing_analysis['missing_bp_clusters'].most_common(1)[0]
            print(f"  • Most incomplete BP cluster: {top_incomplete[0]} (incomplete in {top_incomplete[1]} strains)")
        
        if absent_analysis['absent_bp_clusters']:
            top_absent = absent_analysis['absent_bp_clusters'].most_common(1)[0]
            print(f"  • Most absent BP cluster: {top_absent[0]} (completely absent in {top_absent[1]} strains)")
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
