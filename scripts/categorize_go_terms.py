#!/usr/bin/env python3
"""
GO term categorization using Revigo clustering with nuanced analysis:
- Core-only: GO clusters only in core genes (no expansion)
- Core-expanded: GO clusters in core genes + additional accessory copies
- Accessory-only: GO clusters in non-core genes shared by ≥2 genomes
- Unique: GO clusters found in only 1 genome
"""

import argparse
import pandas as pd
import numpy as np
from pathlib import Path
from collections import defaultdict, Counter
import sys

# Import custom sparse utilities
try:
    import sparse_utils
except ImportError:
    print("ERROR: sparse_utils module not found. Please ensure it's in your PYTHONPATH.")
    sys.exit(1)

def parse_go_terms(go_terms_str):
    """Parse GO terms from string format."""
    if pd.isna(go_terms_str) or str(go_terms_str).strip() == '' or str(go_terms_str) == 'nan':
        return []
    
    go_terms = []
    go_str = str(go_terms_str).strip()
    
    # Handle different formats
    if ';' in go_str:
        terms = go_str.split(';')
        for term in terms:
            term = term.strip()
            if '[GO:' in term:
                desc = term.split('[GO:')[0].strip()
                if desc and desc != '':
                    go_terms.append(desc)
            elif term and term != '':
                go_terms.append(term)
    else:
        if '[GO:' in go_str:
            desc = go_str.split('[GO:')[0].strip()
            if desc and desc != '':
                go_terms.append(desc)
        elif go_str and go_str != '':
            go_terms.append(go_str)
    
    return go_terms

def load_revigo_clusters(revigo_file):
    """Load and process Revigo clustering data."""
    print(f"Loading Revigo clusters from {revigo_file}")
    
    # Load Revigo data
    revigo_df = pd.read_csv(revigo_file, sep='\t')
    
    # Create cluster mapping
    go_to_cluster = {}
    cluster_to_members = defaultdict(set)
    cluster_representatives = {}
    
    for _, row in revigo_df.iterrows():
        go_id = row['TermID'].strip('"')
        go_name = row['Name'].strip('"')
        representative = row['Representative']
        
        # Determine cluster representative
        if pd.isna(representative) or representative == 'null':
            # This term is a representative
            cluster_id = go_id
            cluster_representatives[cluster_id] = go_name
        else:
            # This term belongs to another cluster
            cluster_id = f"GO:{int(representative):07d}"
        
        go_to_cluster[go_name] = cluster_id
        cluster_to_members[cluster_id].add(go_name)
    
    print(f"  Found {len(cluster_representatives)} clusters with {len(go_to_cluster)} total GO terms")
    
    return go_to_cluster, cluster_to_members, cluster_representatives

def load_proteome_metadata(metadata_file):
    """Load proteome metadata for organism names."""
    print("Loading proteome metadata...")
    metadata = pd.read_csv(metadata_file, sep='\t')
    print(f"  Loaded metadata for {len(metadata)} proteomes")
    
    # Create mapping from proteome ID to organism name
    proteome_to_organism = {}
    for _, row in metadata.iterrows():
        proteome_id = row['Proteome Id']
        organism = row['Organism']
        proteome_to_organism[proteome_id] = organism
    
    return proteome_to_organism

def load_data(all_annotations_file, gene_matrix_file, gene_labels_file, core_genes_file):
    """Load all required data."""
    print("\nLoading annotation data...")
    
    # Load all dominant alleles annotations
    all_annotations = pd.read_csv(all_annotations_file, sep='\t')
    print(f"  Loaded annotations for {len(all_annotations)} dominant alleles")
    
    # Load gene matrix to get strain information
    try:
        df_genes = sparse_utils.read_lsdf(gene_matrix_file, gene_labels_file)
        print(f"  Loaded gene matrix: {df_genes.shape}")
    except Exception as e:
        print(f"  WARNING: Could not load gene matrix: {e}")
        df_genes = None
    
    # Load core genes
    with open(core_genes_file, 'r') as f:
        core_genes = set(line.strip() for line in f if line.strip())
    print(f"  Loaded {len(core_genes)} core genes")
    
    # Debug: Count total genes in matrix
    if df_genes is not None:
        total_genes = len(df_genes.index_map)
        core_percentage = (len(core_genes) / total_genes) * 100 if total_genes > 0 else 0
        print(f"  Total genes in matrix: {total_genes}")
        print(f"  Core genes percentage: {core_percentage:.2f}%")
    
    return all_annotations, df_genes, core_genes

def analyze_go_clusters(all_annotations, df_genes, core_genes, proteome_to_organism,
                       go_to_cluster_bp, cluster_to_members_bp, cluster_representatives_bp,
                       go_to_cluster_mf, cluster_to_members_mf, cluster_representatives_mf):
    """Analyze GO clusters with nuanced categorization."""
    print("\nAnalyzing GO clusters...")
    
    # Track which strains have each cluster via core and non-core genes
    cluster_to_strains_via_core = {'BP': defaultdict(set), 'MF': defaultdict(set)}
    cluster_to_strains_via_noncore = {'BP': defaultdict(set), 'MF': defaultdict(set)}
    cluster_to_all_strains = {'BP': defaultdict(set), 'MF': defaultdict(set)}
    
    # First pass: identify which genes have which GO terms
    gene_to_go = {'BP': defaultdict(set), 'MF': defaultdict(set)}
    
    annotated_count = 0
    unmapped_go_terms = {'BP': set(), 'MF': set()}
    
    for _, row in all_annotations.iterrows():
        gene_id = row['gene']
        
        # Parse BP terms
        bp_terms = parse_go_terms(row.get('Gene Ontology (biological process)', ''))
        for term in bp_terms:
            gene_to_go['BP'][gene_id].add(term)
            if term not in go_to_cluster_bp:
                unmapped_go_terms['BP'].add(term)
        
        # Parse MF terms
        mf_terms = parse_go_terms(row.get('Gene Ontology (molecular function)', ''))
        for term in mf_terms:
            gene_to_go['MF'][gene_id].add(term)
            if term not in go_to_cluster_mf:
                unmapped_go_terms['MF'].add(term)
        
        if bp_terms or mf_terms:
            annotated_count += 1
    
    print(f"  Found {annotated_count} genes with GO annotations")
    
    if unmapped_go_terms['BP']:
        print(f"  WARNING: {len(unmapped_go_terms['BP'])} BP terms not found in Revigo clustering")
    if unmapped_go_terms['MF']:
        print(f"  WARNING: {len(unmapped_go_terms['MF'])} MF terms not found in Revigo clustering")
    
    # Get total number of strains
    total_strains = len(df_genes.columns) if df_genes is not None else 0
    print(f"\nTotal strains in analysis: {total_strains}")
    
    # Second pass: map genes to strains and track cluster presence
    if df_genes is not None:
        print("Mapping clusters to strains...")
        
        # Pre-convert to CSR for efficiency
        csr_matrix = df_genes.data.tocsr()
        
        # Process only genes that have GO terms
        genes_with_go = set()
        for go_type in ['BP', 'MF']:
            genes_with_go.update(gene_to_go[go_type].keys())
        
        processed = 0
        core_genes_with_go = 0
        
        for gene_id in genes_with_go:
            is_core = gene_id in core_genes
            if is_core:
                core_genes_with_go += 1
            
            if gene_id in df_genes.index_map:
                # Get gene position
                gene_pos = df_genes.index_map[gene_id]
                # Get row from pre-converted CSR matrix
                gene_row = csr_matrix.getrow(gene_pos)
                
                # Get the strain indices where this gene is present
                strain_indices = gene_row.nonzero()[1]
                
                # Map indices to strain names
                strains_with_gene = set()
                for idx in strain_indices:
                    if idx < len(df_genes.columns):
                        strains_with_gene.add(df_genes.columns[idx])
                
                # Map GO terms to clusters and track strains
                # BP terms
                for term in gene_to_go['BP'].get(gene_id, []):
                    if term in go_to_cluster_bp:
                        cluster_id = go_to_cluster_bp[term]
                        cluster_to_all_strains['BP'][cluster_id].update(strains_with_gene)
                        
                        if is_core:
                            cluster_to_strains_via_core['BP'][cluster_id].update(strains_with_gene)
                        else:
                            cluster_to_strains_via_noncore['BP'][cluster_id].update(strains_with_gene)
                
                # MF terms
                for term in gene_to_go['MF'].get(gene_id, []):
                    if term in go_to_cluster_mf:
                        cluster_id = go_to_cluster_mf[term]
                        cluster_to_all_strains['MF'][cluster_id].update(strains_with_gene)
                        
                        if is_core:
                            cluster_to_strains_via_core['MF'][cluster_id].update(strains_with_gene)
                        else:
                            cluster_to_strains_via_noncore['MF'][cluster_id].update(strains_with_gene)
                
                processed += 1
                if processed % 1000 == 0:
                    print(f"    Processed {processed}/{len(genes_with_go)} genes...")
        
        print(f"\n  Debug: Core gene GO annotation analysis:")
        print(f"    Core genes with GO annotations: {core_genes_with_go}/{len(core_genes)} ({core_genes_with_go/len(core_genes)*100:.1f}% of core genes)")
        print(f"    Core genes as % of annotated genes: {core_genes_with_go/len(genes_with_go)*100:.1f}%")
    else:
        print("  WARNING: Cannot map clusters to strains without gene matrix")
    
    # Categorize clusters with nuanced approach
    categories = {
        'BP': {'core_only': set(), 'core_expanded': set(), 'accessory_only': set(), 'unique': set()},
        'MF': {'core_only': set(), 'core_expanded': set(), 'accessory_only': set(), 'unique': set()}
    }
    
    unique_cluster_to_strain = {'BP': {}, 'MF': {}}
    
    # Analyze each cluster
    for go_type in ['BP', 'MF']:
        print(f"\nAnalyzing {go_type} clusters...")
        
        for cluster_id in cluster_to_all_strains[go_type]:
            total_strains_with_cluster = len(cluster_to_all_strains[go_type][cluster_id])
            strains_via_core = cluster_to_strains_via_core[go_type].get(cluster_id, set())
            strains_via_noncore = cluster_to_strains_via_noncore[go_type].get(cluster_id, set())
            
            # Check if it's in core genes
            is_in_core = len(strains_via_core) > 0
            
            # Check if it's present in all strains
            is_universal = total_strains_with_cluster == total_strains
            
            # Check for expansion beyond core
            has_expansion = len(strains_via_noncore) > 0
            
            # Categorize
            if is_in_core and is_universal and not has_expansion:
                # Present in all strains only via core genes
                categories[go_type]['core_only'].add(cluster_id)
            elif is_in_core and has_expansion:
                # Present in core genes + additional copies
                categories[go_type]['core_expanded'].add(cluster_id)
            elif not is_in_core and total_strains_with_cluster >= 2:
                # Not in core, shared by multiple strains
                categories[go_type]['accessory_only'].add(cluster_id)
            elif not is_in_core and total_strains_with_cluster == 1:
                # Not in core, unique to one strain
                categories[go_type]['unique'].add(cluster_id)
                strain_id = list(cluster_to_all_strains[go_type][cluster_id])[0]
                organism = proteome_to_organism.get(strain_id, strain_id)
                unique_cluster_to_strain[go_type][cluster_id] = {
                    'strain_id': strain_id,
                    'organism': organism
                }
    
    # Print analysis summary
    for go_type in ['BP', 'MF']:
        print(f"\n{go_type} cluster categories:")
        print(f"  Core-only: {len(categories[go_type]['core_only'])}")
        print(f"  Core-expanded: {len(categories[go_type]['core_expanded'])}")
        print(f"  Accessory-only: {len(categories[go_type]['accessory_only'])}")
        print(f"  Unique: {len(categories[go_type]['unique'])}")
    
    return (categories, cluster_to_all_strains, cluster_to_strains_via_core, 
            cluster_to_strains_via_noncore, unique_cluster_to_strain, 
            cluster_to_members_bp, cluster_representatives_bp,
            cluster_to_members_mf, cluster_representatives_mf)

def save_results(categories, cluster_to_all_strains, cluster_to_strains_via_core,
                 cluster_to_strains_via_noncore, unique_cluster_to_strain,
                 cluster_to_members_bp, cluster_representatives_bp,
                 cluster_to_members_mf, cluster_representatives_mf,
                 proteome_to_organism, output_dir):
    """Save categorization results with nuanced analysis."""
    print("\nSaving results...")
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Import matplotlib for visualizations
    try:
        import matplotlib.pyplot as plt
        import seaborn as sns
        plt.style.use('seaborn-v0_8-darkgrid')
        CAN_PLOT = True
    except ImportError:
        print("WARNING: matplotlib/seaborn not available, skipping visualizations")
        CAN_PLOT = False
    
    # Helper function to get cluster name
    def get_cluster_name(cluster_id, representatives):
        return representatives.get(cluster_id, f"Cluster {cluster_id}")
    
    # Save categorized GO clusters
    for go_type, cluster_to_members, cluster_representatives in [
        ('BP', cluster_to_members_bp, cluster_representatives_bp),
        ('MF', cluster_to_members_mf, cluster_representatives_mf)
    ]:
        go_name = 'biological_process' if go_type == 'BP' else 'molecular_function'
        
        # Core-only clusters
        core_only_data = []
        for cluster_id in sorted(categories[go_type]['core_only']):
            members = cluster_to_members.get(cluster_id, set())
            core_only_data.append({
                'cluster_id': cluster_id,
                'cluster_representative': get_cluster_name(cluster_id, cluster_representatives),
                'category': 'core_only',
                'num_go_terms': len(members),
                'go_term_members': '; '.join(sorted(members))
            })
        core_only_df = pd.DataFrame(core_only_data)
        core_only_df.to_csv(output_dir / f'go_{go_type.lower()}_core_only_clusters.tsv', 
                           sep='\t', index=False)
        
        # Core-expanded clusters
        core_expanded_data = []
        for cluster_id in sorted(categories[go_type]['core_expanded']):
            members = cluster_to_members.get(cluster_id, set())
            total_strains = len(cluster_to_all_strains[go_type][cluster_id])
            expanded_strains = len(cluster_to_strains_via_noncore[go_type].get(cluster_id, set()))
            core_expanded_data.append({
                'cluster_id': cluster_id,
                'cluster_representative': get_cluster_name(cluster_id, cluster_representatives),
                'category': 'core_expanded',
                'total_strains': total_strains,
                'strains_with_expansion': expanded_strains,
                'num_go_terms': len(members),
                'go_term_members': '; '.join(sorted(members))
            })
        core_expanded_df = pd.DataFrame(core_expanded_data)
        core_expanded_df.to_csv(output_dir / f'go_{go_type.lower()}_core_expanded_clusters.tsv', 
                               sep='\t', index=False)
        
        # Accessory-only clusters
        accessory_data = []
        for cluster_id in sorted(categories[go_type]['accessory_only']):
            members = cluster_to_members.get(cluster_id, set())
            num_strains = len(cluster_to_all_strains[go_type][cluster_id])
            accessory_data.append({
                'cluster_id': cluster_id,
                'cluster_representative': get_cluster_name(cluster_id, cluster_representatives),
                'category': 'accessory_only',
                'genome_count': num_strains,
                'num_go_terms': len(members),
                'go_term_members': '; '.join(sorted(members))
            })
        accessory_df = pd.DataFrame(accessory_data)
        if not accessory_df.empty:
            accessory_df = accessory_df.sort_values('genome_count', ascending=False)
        accessory_df.to_csv(output_dir / f'go_{go_type.lower()}_accessory_only_clusters.tsv', 
                           sep='\t', index=False)
        
        # Unique clusters
        unique_data = []
        for cluster_id in sorted(categories[go_type]['unique']):
            members = cluster_to_members.get(cluster_id, set())
            strain_info = unique_cluster_to_strain[go_type].get(cluster_id, {})
            unique_data.append({
                'cluster_id': cluster_id,
                'cluster_representative': get_cluster_name(cluster_id, cluster_representatives),
                'category': 'unique',
                'strain_id': strain_info.get('strain_id', 'Unknown'),
                'organism': strain_info.get('organism', 'Unknown'),
                'num_go_terms': len(members),
                'go_term_members': '; '.join(sorted(members))
            })
        unique_df = pd.DataFrame(unique_data)
        unique_df.to_csv(output_dir / f'go_{go_type.lower()}_unique_clusters.tsv', 
                        sep='\t', index=False)
    
    # Track strains with most unique and accessory clusters
    print("\nAnalyzing strain-specific patterns...")
    
    # Count unique clusters per strain - separate by BP and MF
    unique_clusters_per_strain_bp = Counter()
    unique_clusters_per_strain_mf = Counter()
    unique_clusters_per_strain_total = Counter()
    
    for cluster_id, strain_info in unique_cluster_to_strain['BP'].items():
        strain_id = strain_info['strain_id']
        unique_clusters_per_strain_bp[strain_id] += 1
        unique_clusters_per_strain_total[strain_id] += 1
    
    for cluster_id, strain_info in unique_cluster_to_strain['MF'].items():
        strain_id = strain_info['strain_id']
        unique_clusters_per_strain_mf[strain_id] += 1
        unique_clusters_per_strain_total[strain_id] += 1
    
    # Count non-core clusters per strain (accessory + expanded) - separate by BP and MF
    noncore_clusters_per_strain_bp = Counter()
    noncore_clusters_per_strain_mf = Counter()
    noncore_clusters_per_strain_total = Counter()
    
    # BP non-core clusters
    for cluster_id in categories['BP']['accessory_only']:
        strains = cluster_to_all_strains['BP'].get(cluster_id, set())
        for strain_id in strains:
            noncore_clusters_per_strain_bp[strain_id] += 1
            noncore_clusters_per_strain_total[strain_id] += 1
    
    for cluster_id in categories['BP']['core_expanded']:
        strains = cluster_to_strains_via_noncore['BP'].get(cluster_id, set())
        for strain_id in strains:
            noncore_clusters_per_strain_bp[strain_id] += 1
            noncore_clusters_per_strain_total[strain_id] += 1
    
    # MF non-core clusters
    for cluster_id in categories['MF']['accessory_only']:
        strains = cluster_to_all_strains['MF'].get(cluster_id, set())
        for strain_id in strains:
            noncore_clusters_per_strain_mf[strain_id] += 1
            noncore_clusters_per_strain_total[strain_id] += 1
    
    for cluster_id in categories['MF']['core_expanded']:
        strains = cluster_to_strains_via_noncore['MF'].get(cluster_id, set())
        for strain_id in strains:
            noncore_clusters_per_strain_mf[strain_id] += 1
            noncore_clusters_per_strain_total[strain_id] += 1
    
    # Create strain analysis report
    strain_analysis_data = []
    all_strains = set(unique_clusters_per_strain_total.keys()) | set(noncore_clusters_per_strain_total.keys())
    
    for strain_id in all_strains:
        organism = proteome_to_organism.get(strain_id, strain_id)
        strain_analysis_data.append({
            'strain_id': strain_id,
            'organism': organism,
            'unique_clusters_BP': unique_clusters_per_strain_bp.get(strain_id, 0),
            'unique_clusters_MF': unique_clusters_per_strain_mf.get(strain_id, 0),
            'unique_clusters_total': unique_clusters_per_strain_total.get(strain_id, 0),
            'noncore_clusters_BP': noncore_clusters_per_strain_bp.get(strain_id, 0),
            'noncore_clusters_MF': noncore_clusters_per_strain_mf.get(strain_id, 0),
            'noncore_clusters_total': noncore_clusters_per_strain_total.get(strain_id, 0),
            'total_specialized': unique_clusters_per_strain_total.get(strain_id, 0) + 
                               noncore_clusters_per_strain_total.get(strain_id, 0)
        })
    
    strain_analysis_df = pd.DataFrame(strain_analysis_data)
    strain_analysis_df = strain_analysis_df.sort_values('unique_clusters_total', ascending=False)
    strain_analysis_df.to_csv(output_dir / 'strain_cluster_analysis.tsv', sep='\t', index=False)
    
    # Create visualizations if matplotlib is available
    if CAN_PLOT:
        print("Creating visualizations...")
        
        # 1. Bar plot of strains with most unique clusters - separate BP and MF
        if len(strain_analysis_df) > 0:
            fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
            
            # Top 15 strains with most unique BP clusters
            top_unique_bp = strain_analysis_df.nlargest(15, 'unique_clusters_BP')
            top_unique_bp['short_name'] = top_unique_bp['organism'].apply(
                lambda x: x.split('(')[0].strip() if '(' in x else x
            )
            
            # Plot unique BP clusters
            ax1.barh(range(len(top_unique_bp)), top_unique_bp['unique_clusters_BP'], color='#e74c3c')
            ax1.set_yticks(range(len(top_unique_bp)))
            ax1.set_yticklabels(top_unique_bp['short_name'])
            ax1.set_xlabel('Number of Unique GO Clusters')
            ax1.set_title('Top 15 Strains - Unique Biological Process Clusters', fontsize=14, pad=10)
            ax1.invert_yaxis()
            
            # Add value labels
            for i, v in enumerate(top_unique_bp['unique_clusters_BP']):
                if v > 0:
                    ax1.text(v + 0.1, i, str(v), va='center')
            
            # Top 15 strains with most unique MF clusters
            top_unique_mf = strain_analysis_df.nlargest(15, 'unique_clusters_MF')
            top_unique_mf['short_name'] = top_unique_mf['organism'].apply(
                lambda x: x.split('(')[0].strip() if '(' in x else x
            )
            
            # Plot unique MF clusters
            ax2.barh(range(len(top_unique_mf)), top_unique_mf['unique_clusters_MF'], color='#e67e22')
            ax2.set_yticks(range(len(top_unique_mf)))
            ax2.set_yticklabels(top_unique_mf['short_name'])
            ax2.set_xlabel('Number of Unique GO Clusters')
            ax2.set_title('Top 15 Strains - Unique Molecular Function Clusters', fontsize=14, pad=10)
            ax2.invert_yaxis()
            
            # Add value labels
            for i, v in enumerate(top_unique_mf['unique_clusters_MF']):
                if v > 0:
                    ax2.text(v + 0.1, i, str(v), va='center')
            
            # Top 15 strains with most non-core BP clusters
            top_noncore_bp = strain_analysis_df.nlargest(15, 'noncore_clusters_BP')
            top_noncore_bp['short_name'] = top_noncore_bp['organism'].apply(
                lambda x: x.split('(')[0].strip() if '(' in x else x
            )
            
            # Plot non-core BP clusters
            ax3.barh(range(len(top_noncore_bp)), top_noncore_bp['noncore_clusters_BP'], color='#3498db')
            ax3.set_yticks(range(len(top_noncore_bp)))
            ax3.set_yticklabels(top_noncore_bp['short_name'])
            ax3.set_xlabel('Number of Non-Core GO Clusters')
            ax3.set_title('Top 15 Strains - Non-Core Biological Process Clusters', fontsize=14, pad=10)
            ax3.invert_yaxis()
            
            # Add value labels
            for i, v in enumerate(top_noncore_bp['noncore_clusters_BP']):
                if v > 0:
                    ax3.text(v + 1, i, str(v), va='center')
            
            # Top 15 strains with most non-core MF clusters
            top_noncore_mf = strain_analysis_df.nlargest(15, 'noncore_clusters_MF')
            top_noncore_mf['short_name'] = top_noncore_mf['organism'].apply(
                lambda x: x.split('(')[0].strip() if '(' in x else x
            )
            
            # Plot non-core MF clusters
            ax4.barh(range(len(top_noncore_mf)), top_noncore_mf['noncore_clusters_MF'], color='#2980b9')
            ax4.set_yticks(range(len(top_noncore_mf)))
            ax4.set_yticklabels(top_noncore_mf['short_name'])
            ax4.set_xlabel('Number of Non-Core GO Clusters')
            ax4.set_title('Top 15 Strains - Non-Core Molecular Function Clusters', fontsize=14, pad=10)
            ax4.invert_yaxis()
            
            # Add value labels
            for i, v in enumerate(top_noncore_mf['noncore_clusters_MF']):
                if v > 0:
                    ax4.text(v + 1, i, str(v), va='center')
            
            plt.tight_layout()
            plt.savefig(output_dir / 'strain_cluster_rankings.png', dpi=300, bbox_inches='tight')
            plt.close()
        
        # 2. Combined stacked bar chart - separate BP and MF
        if len(strain_analysis_df) > 0:
            fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(18, 8))
            
            # Add overall title
            fig.suptitle('Top 20 Strains by Specialized GO Functions (Unique + Non-Core)', fontsize=16, y=1.02)
            
            # Get top 20 strains by total specialized functions
            top_strains = strain_analysis_df.nlargest(20, 'total_specialized')
            top_strains['short_name'] = top_strains['organism'].apply(
                lambda x: x.split('(')[0].strip()[:25] if '(' in x else x[:25]
            )
            
            # BP stacked bar chart
            x = range(len(top_strains))
            width = 0.8
            
            p1 = ax1.bar(x, top_strains['unique_clusters_BP'], width, 
                         label='Unique BP', color='#e74c3c')
            p2 = ax1.bar(x, top_strains['noncore_clusters_BP'], width,
                         bottom=top_strains['unique_clusters_BP'],
                         label='Non-Core BP', color='#3498db')
            
            ax1.set_ylabel('Number of GO Clusters', fontsize=12)
            ax1.set_title('Top 20 Strains - Biological Process Functions', fontsize=14, pad=20)
            ax1.set_xticks(x)
            ax1.set_xticklabels(top_strains['short_name'], rotation=45, ha='right')
            ax1.legend()
            
            # Add total value labels on top
            for i, (idx, row) in enumerate(top_strains.iterrows()):
                total = row['unique_clusters_BP'] + row['noncore_clusters_BP']
                if total > 0:
                    ax1.text(i, total + 5, str(int(total)), ha='center', va='bottom', fontsize=8)
            
            # MF stacked bar chart
            p3 = ax2.bar(x, top_strains['unique_clusters_MF'], width, 
                         label='Unique MF', color='#e67e22')
            p4 = ax2.bar(x, top_strains['noncore_clusters_MF'], width,
                         bottom=top_strains['unique_clusters_MF'],
                         label='Non-Core MF', color='#2980b9')
            
            ax2.set_ylabel('Number of GO Clusters', fontsize=12)
            ax2.set_title('Top 20 Strains - Molecular Function', fontsize=14, pad=20)
            ax2.set_xticks(x)
            ax2.set_xticklabels(top_strains['short_name'], rotation=45, ha='right')
            ax2.legend()
            
            # Add total value labels on top
            for i, (idx, row) in enumerate(top_strains.iterrows()):
                total = row['unique_clusters_MF'] + row['noncore_clusters_MF']
                if total > 0:
                    ax2.text(i, total + 5, str(int(total)), ha='center', va='bottom', fontsize=8)
            
            plt.tight_layout()
            plt.savefig(output_dir / 'strain_noncore_clusters_stacked.png', dpi=300, bbox_inches='tight')
            plt.close()
    
    # Create enhanced summary report
    summary_lines = [
        "GO CLUSTER CATEGORIZATION SUMMARY (Nuanced Analysis with Revigo Clustering)",
        "=" * 80,
        "",
        "ANALYSIS APPROACH:",
        "- Core-only: Functions present in all strains ONLY via core genes",
        "- Core-expanded: Functions in core genes + additional accessory copies",
        "- Accessory-only: Functions not in core genes, shared by ≥2 strains",
        "- Unique: Functions not in core genes, found in only 1 strain",
        "",
        "OVERALL STATISTICS:",
        f"Total GO clusters analyzed: {sum(len(s) for cat in categories['BP'].values() for s in cat) + sum(len(s) for cat in categories['MF'].values() for s in cat)}",
        "",
        "BIOLOGICAL PROCESS (BP) CLUSTERS:",
        f"  Core-only: {len(categories['BP']['core_only'])}",
        f"  Core-expanded: {len(categories['BP']['core_expanded'])}",
        f"  Accessory-only: {len(categories['BP']['accessory_only'])}",
        f"  Unique: {len(categories['BP']['unique'])}",
        f"  Total BP clusters: {sum(len(s) for s in categories['BP'].values())}",
        "",
        "MOLECULAR FUNCTION (MF) CLUSTERS:",
        f"  Core-only: {len(categories['MF']['core_only'])}",
        f"  Core-expanded: {len(categories['MF']['core_expanded'])}",
        f"  Accessory-only: {len(categories['MF']['accessory_only'])}",
        f"  Unique: {len(categories['MF']['unique'])}",
        f"  Total MF clusters: {sum(len(s) for s in categories['MF'].values())}",
        "",
        "KEY INSIGHTS:",
        f"  Functions with expansion beyond core: {len(categories['BP']['core_expanded']) + len(categories['MF']['core_expanded'])}",
        f"  Functions restricted to core genes: {len(categories['BP']['core_only']) + len(categories['MF']['core_only'])}",
        "",
        "STRAIN ANALYSIS:",
        ""
    ]
    
    # Add top strains with unique clusters
    if unique_clusters_per_strain_total:
        top_unique_strains = unique_clusters_per_strain_total.most_common(5)
        summary_lines.append("Top 5 strains with most unique GO clusters (BP+MF combined):")
        for strain_id, count in top_unique_strains:
            organism = proteome_to_organism.get(strain_id, strain_id)
            bp_count = unique_clusters_per_strain_bp.get(strain_id, 0)
            mf_count = unique_clusters_per_strain_mf.get(strain_id, 0)
            summary_lines.append(f"  • {organism}: {count} total ({bp_count} BP, {mf_count} MF)")
    
    summary_lines.append("")
    
    # Add top strains with non-core clusters
    if noncore_clusters_per_strain_total:
        top_noncore_strains = noncore_clusters_per_strain_total.most_common(5)
        summary_lines.append("Top 5 strains with most non-core GO clusters (BP+MF combined):")
        for strain_id, count in top_noncore_strains:
            organism = proteome_to_organism.get(strain_id, strain_id)
            bp_count = noncore_clusters_per_strain_bp.get(strain_id, 0)
            mf_count = noncore_clusters_per_strain_mf.get(strain_id, 0)
            summary_lines.append(f"  • {organism}: {count} total ({bp_count} BP, {mf_count} MF)")
    
    summary_lines.extend([
        "",
        "TOP CORE-EXPANDED FUNCTIONS:",
        "",
        "Biological Process (most expanded):"
    ])
    
    # Add top core-expanded BP functions
    if categories['BP']['core_expanded']:
        bp_expanded_sorted = sorted(
            [(cluster_id, len(cluster_to_strains_via_noncore['BP'].get(cluster_id, set()))) 
             for cluster_id in categories['BP']['core_expanded']],
            key=lambda x: x[1], reverse=True
        )[:5]
        
        for i, (cluster_id, expanded_count) in enumerate(bp_expanded_sorted, 1):
            cluster_name = get_cluster_name(cluster_id, cluster_representatives_bp)
            num_terms = len(cluster_to_members_bp.get(cluster_id, set()))
            summary_lines.append(f"  {i}. {cluster_name} (+{expanded_count} strains with extra copies, {num_terms} GO terms)")
    
    summary_lines.extend(["", "Molecular Function (most expanded):"])
    
    # Add top core-expanded MF functions
    if categories['MF']['core_expanded']:
        mf_expanded_sorted = sorted(
            [(cluster_id, len(cluster_to_strains_via_noncore['MF'].get(cluster_id, set()))) 
             for cluster_id in categories['MF']['core_expanded']],
            key=lambda x: x[1], reverse=True
        )[:5]
        
        for i, (cluster_id, expanded_count) in enumerate(mf_expanded_sorted, 1):
            cluster_name = get_cluster_name(cluster_id, cluster_representatives_mf)
            num_terms = len(cluster_to_members_mf.get(cluster_id, set()))
            summary_lines.append(f"  {i}. {cluster_name} (+{expanded_count} strains with extra copies, {num_terms} GO terms)")
    
    # Save summary
    with open(output_dir / 'go_cluster_categorization_summary.txt', 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"Results saved to {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Categorize GO clusters with nuanced analysis')
    parser.add_argument('--all-annotations', required=True, 
                       help='All dominant alleles annotations TSV file')
    parser.add_argument('--gene-matrix', required=True, 
                       help='Gene matrix NPZ file')
    parser.add_argument('--gene-labels', required=True, 
                       help='Gene labels file')
    parser.add_argument('--core-genes', required=True, 
                       help='Core genes list file')
    parser.add_argument('--revigo-all-bp', required=True,
                       help='Revigo BP table for all GO terms')
    parser.add_argument('--revigo-all-mf', required=True,
                       help='Revigo MF table for all GO terms')
    parser.add_argument('--proteome-metadata', required=True,
                       help='Proteome metadata TSV file')
    parser.add_argument('--output-dir', required=True, 
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Load Revigo clustering data
    (go_to_cluster_bp, cluster_to_members_bp, 
     cluster_representatives_bp) = load_revigo_clusters(args.revigo_all_bp)
    
    (go_to_cluster_mf, cluster_to_members_mf, 
     cluster_representatives_mf) = load_revigo_clusters(args.revigo_all_mf)
    
    # Load proteome metadata
    proteome_to_organism = load_proteome_metadata(args.proteome_metadata)
    
    # Load annotation data
    all_annotations, df_genes, core_genes = load_data(
        args.all_annotations, args.gene_matrix, 
        args.gene_labels, args.core_genes
    )
    
    # Analyze GO clusters with nuanced approach
    (categories, cluster_to_all_strains, cluster_to_strains_via_core,
     cluster_to_strains_via_noncore, unique_cluster_to_strain,
     cluster_to_members_bp, cluster_representatives_bp,
     cluster_to_members_mf, cluster_representatives_mf) = analyze_go_clusters(
        all_annotations, df_genes, core_genes, proteome_to_organism,
        go_to_cluster_bp, cluster_to_members_bp, cluster_representatives_bp,
        go_to_cluster_mf, cluster_to_members_mf, cluster_representatives_mf
    )
    
    # Save results
    save_results(categories, cluster_to_all_strains, cluster_to_strains_via_core,
                 cluster_to_strains_via_noncore, unique_cluster_to_strain,
                 cluster_to_members_bp, cluster_representatives_bp,
                 cluster_to_members_mf, cluster_representatives_mf,
                 proteome_to_organism, args.output_dir)
    
    # Print summary
    print("\nNuanced Analysis Summary:")
    for go_type in ['BP', 'MF']:
        go_name = 'Biological Process' if go_type == 'BP' else 'Molecular Function'
        print(f"\n{go_name} Clusters:")
        print(f"  Core-only: {len(categories[go_type]['core_only'])}")
        print(f"  Core-expanded: {len(categories[go_type]['core_expanded'])}")
        print(f"  Accessory-only: {len(categories[go_type]['accessory_only'])}")
        print(f"  Unique: {len(categories[go_type]['unique'])}")

if __name__ == '__main__':
    main()
