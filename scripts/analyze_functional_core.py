#!/usr/bin/env python3
"""
Functional Core Genome Analysis
Focuses on functional annotation analysis using GO terms, protein families, and pathways.
Handles: GO enrichment, functional missing genes, pathway analysis, protein family clustering.
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
from scipy import stats
from itertools import combinations
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

def load_data(gene_matrix, gene_labels, core_genes_file, annotations_file, metadata_file=None):
    """Load all required data files for functional analysis."""
    print("Loading data files for functional analysis...")
    
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
    
    # Load metadata if provided
    metadata = None
    if metadata_file and Path(metadata_file).exists():
        try:
            metadata = pd.read_csv(metadata_file, sep='\t')
            print(f"Loaded metadata for {len(metadata)} proteomes")
        except Exception as e:
            print(f"Warning: Could not load metadata: {e}")
    
    return df_genes, core_genes, annotations, metadata

def get_full_organism_name(proteome_id, metadata):
    """Get complete organism name from metadata - preserves full names."""
    if metadata is None:
        return proteome_id
    
    matching_row = metadata[metadata['Proteome Id'] == proteome_id]
    if not matching_row.empty:
        organism_name = matching_row.iloc[0]['Organism']
        return organism_name  # Return complete name as-is
    
    return proteome_id

def parse_go_terms(go_terms_str):
    """Parse GO terms from string format."""
    if pd.isna(go_terms_str) or str(go_terms_str).strip() == '' or str(go_terms_str) == 'nan':
        return []
    
    go_terms = []
    go_str = str(go_terms_str)
    
    if ';' in go_str:
        terms = go_str.split(';')
        for term in terms:
            term = term.strip()
            if '[GO:' in term:
                desc = term.split('[GO:')[0].strip()
                if desc and desc != '':
                    go_terms.append(desc)
    else:
        if '[GO:' in go_str:
            desc = go_str.split('[GO:')[0].strip()
            if desc and desc != '':
                go_terms.append(desc)
    
    return go_terms

def create_functional_mappings(annotations):
    """Create comprehensive functional mappings from annotations."""
    print("Creating functional mappings...")
    
    gene_mappings = {
        'gene_names': {},
        'protein_names': {},
        'ec_numbers': {},
        'go_bp': {},       # Biological Process
        'go_mf': {},       # Molecular Function
        'go_cc': {},       # Cellular Component
        'protein_families': {},
        'organisms': {}
    }
    
    for _, row in annotations.iterrows():
        gene_id = row['gene']
        
        # Gene names
        if pd.notna(row.get('Gene Names (primary)', '')):
            gene_mappings['gene_names'][gene_id] = str(row['Gene Names (primary)']).strip()
        elif pd.notna(row.get('Gene Names', '')):
            gene_names_str = str(row['Gene Names']).strip()
            if gene_names_str and gene_names_str != 'nan':
                gene_mappings['gene_names'][gene_id] = gene_names_str.split()[0]
        
        # Protein names
        if pd.notna(row.get('Protein names', '')):
            protein_name = str(row['Protein names']).strip()
            if protein_name != 'nan':
                gene_mappings['protein_names'][gene_id] = protein_name
        
        # EC numbers
        if pd.notna(row.get('EC number', '')):
            ec_str = str(row['EC number']).strip()
            if ec_str != 'nan':
                gene_mappings['ec_numbers'][gene_id] = ec_str
        
        # GO terms
        for go_type, column in [('go_bp', 'Gene Ontology (biological process)'), 
                               ('go_mf', 'Gene Ontology (molecular function)'),
                               ('go_cc', 'Gene Ontology (cellular component)')]:
            go_terms = parse_go_terms(row.get(column, ''))
            if go_terms:
                gene_mappings[go_type][gene_id] = go_terms
        
        # Organisms
        if pd.notna(row.get('Organism', '')):
            organism = str(row['Organism']).strip()
            if organism != 'nan':
                gene_mappings['organisms'][gene_id] = organism
    
    print(f"Functional mappings created:")
    for mapping_type, mapping_dict in gene_mappings.items():
        print(f"  {mapping_type}: {len(mapping_dict)} genes")
    
    return gene_mappings

def analyze_core_functional_composition(df_genes, core_genes, gene_mappings):
    """Analyze functional composition of core genes."""
    print("\n=== Core Functional Composition Analysis ===")
    
    core_genes_in_matrix = [g for g in core_genes if g in df_genes.index]
    
    # Analyze GO term distribution in core genes
    go_analysis = {}
    for go_type in ['go_bp', 'go_mf', 'go_cc']:
        go_terms_counter = Counter()
        
        for gene in core_genes_in_matrix:
            if gene in gene_mappings[go_type]:
                go_terms_counter.update(gene_mappings[go_type][gene])
        
        go_analysis[go_type] = {
            'total_terms': len(go_terms_counter),
            'total_annotations': sum(go_terms_counter.values()),
            'most_common': go_terms_counter.most_common(20)
        }
    
    # Analyze protein families
    protein_families = Counter()
    ec_numbers = Counter()
    
    for gene in core_genes_in_matrix:
        if gene in gene_mappings['protein_names']:
            # Extract protein family keywords
            protein_name = gene_mappings['protein_names'][gene].lower()
            
            # Common protein family keywords
            family_keywords = ['kinase', 'phosphatase', 'transferase', 'hydrolase', 'oxidoreductase',
                              'ligase', 'lyase', 'isomerase', 'transporter', 'ribosomal', 'polymerase',
                              'protease', 'synthase', 'dehydrogenase', 'reductase', 'transcriptional',
                              'regulatory', 'binding', 'chaperone', 'peptidase']
            
            for keyword in family_keywords:
                if keyword in protein_name:
                    protein_families[keyword] += 1
        
        if gene in gene_mappings['ec_numbers']:
            ec_number = gene_mappings['ec_numbers'][gene]
            # Get EC class (first digit)
            ec_class = ec_number.split('.')[0] if '.' in ec_number else ec_number
            ec_numbers[ec_class] += 1
    
    functional_composition = {
        'go_analysis': go_analysis,
        'protein_families': protein_families.most_common(15),
        'ec_classes': ec_numbers.most_common()
    }
    
    print(f"Core functional composition:")
    print(f"  GO Biological Process terms: {go_analysis['go_bp']['total_terms']}")
    print(f"  GO Molecular Function terms: {go_analysis['go_mf']['total_terms']}")
    print(f"  GO Cellular Component terms: {go_analysis['go_cc']['total_terms']}")
    print(f"  Protein family keywords: {len(protein_families)}")
    print(f"  EC number classes: {len(ec_numbers)}")
    
    return functional_composition

def analyze_missing_genes_functional(df_genes, core_genes, gene_mappings, metadata):
    """Analyze functional patterns in missing core genes."""
    print("\n=== Missing Genes Functional Analysis ===")
    
    core_genes_in_matrix = [g for g in core_genes if g in df_genes.index]
    core_df = df_genes.labelslice(indices=core_genes_in_matrix)
    
    # Identify missing genes per strain
    missing_functional_analysis = {}
    strain_functional_profiles = {}
    
    for strain_idx, strain_id in enumerate(df_genes.columns):
        strain_column = core_df.data.tocsc().getcol(strain_idx)
        present_gene_indices = strain_column.nonzero()[0]
        present_genes = set(core_genes_in_matrix[i] for i in present_gene_indices)
        missing_genes = set(core_genes_in_matrix) - present_genes
        
        if len(missing_genes) == 0:
            continue
        
        # Analyze GO terms in missing genes
        missing_go_analysis = {'go_bp': Counter(), 'go_mf': Counter(), 'go_cc': Counter()}
        missing_protein_families = Counter()
        missing_ec_classes = Counter()
        
        for gene in missing_genes:
            # GO terms
            for go_type in ['go_bp', 'go_mf', 'go_cc']:
                if gene in gene_mappings[go_type]:
                    missing_go_analysis[go_type].update(gene_mappings[go_type][gene])
            
            # Protein families
            if gene in gene_mappings['protein_names']:
                protein_name = gene_mappings['protein_names'][gene].lower()
                family_keywords = ['kinase', 'phosphatase', 'transferase', 'hydrolase', 'oxidoreductase',
                                  'ligase', 'lyase', 'isomerase', 'transporter', 'ribosomal', 'polymerase',
                                  'protease', 'synthase', 'dehydrogenase', 'reductase', 'transcriptional']
                
                for keyword in family_keywords:
                    if keyword in protein_name:
                        missing_protein_families[keyword] += 1
            
            # EC classes
            if gene in gene_mappings['ec_numbers']:
                ec_number = gene_mappings['ec_numbers'][gene]
                ec_class = ec_number.split('.')[0] if '.' in ec_number else ec_number
                missing_ec_classes[ec_class] += 1
        
        organism = get_full_organism_name(strain_id, metadata)
        
        strain_functional_profiles[strain_id] = {
            'organism': organism,
            'missing_count': len(missing_genes),
            'missing_go_bp': missing_go_analysis['go_bp'],
            'missing_go_mf': missing_go_analysis['go_mf'],
            'missing_go_cc': missing_go_analysis['go_cc'],
            'missing_protein_families': missing_protein_families,
            'missing_ec_classes': missing_ec_classes
        }
    
    # Aggregate functional analysis across all strains
    aggregated_missing_functions = {
        'go_bp': Counter(),
        'go_mf': Counter(),
        'go_cc': Counter(),
        'protein_families': Counter(),
        'ec_classes': Counter()
    }
    
    for strain_data in strain_functional_profiles.values():
        aggregated_missing_functions['go_bp'].update(strain_data['missing_go_bp'])
        aggregated_missing_functions['go_mf'].update(strain_data['missing_go_mf'])
        aggregated_missing_functions['go_cc'].update(strain_data['missing_go_cc'])
        aggregated_missing_functions['protein_families'].update(strain_data['missing_protein_families'])
        aggregated_missing_functions['ec_classes'].update(strain_data['missing_ec_classes'])
    
    missing_functional_analysis = {
        'strain_profiles': strain_functional_profiles,
        'aggregated_functions': aggregated_missing_functions
    }
    
    print(f"Missing genes functional analysis:")
    print(f"  Strains with functional missing analysis: {len(strain_functional_profiles)}")
    print(f"  Most common missing BP terms: {len(aggregated_missing_functions['go_bp'])}")
    print(f"  Most common missing MF terms: {len(aggregated_missing_functions['go_mf'])}")
    print(f"  Most common missing protein families: {len(aggregated_missing_functions['protein_families'])}")
    
    return missing_functional_analysis

def analyze_functional_enrichment(df_genes, core_genes, gene_mappings):
    """Analyze functional enrichment in different gene categories."""
    print("\n=== Functional Enrichment Analysis ===")
    
    n_strains = len(df_genes.columns)
    
    # Calculate gene frequencies
    gene_counts = df_genes.sum(axis='index')
    
    # Define gene categories
    core_genes_set = set(core_genes)
    unique_genes = set(df_genes.index[gene_counts == 1])
    rare_genes = set(df_genes.index[(gene_counts > 1) & (gene_counts <= 0.05 * n_strains)])
    accessory_genes = set(df_genes.index[(gene_counts > 0.05 * n_strains) & (gene_counts < 0.95 * n_strains)])
    
    # Remove core genes from accessory for cleaner analysis
    accessory_genes = accessory_genes - core_genes_set
    
    categories = {
        'Core': core_genes_set & set(df_genes.index),
        'Accessory': accessory_genes,
        'Rare': rare_genes,
        'Unique': unique_genes
    }
    
    # Analyze GO term enrichment for each category
    enrichment_analysis = {}
    
    for category_name, gene_set in categories.items():
        category_go_analysis = {'go_bp': Counter(), 'go_mf': Counter(), 'go_cc': Counter()}
        category_protein_families = Counter()
        
        for gene in gene_set:
            # GO terms
            for go_type in ['go_bp', 'go_mf', 'go_cc']:
                if gene in gene_mappings[go_type]:
                    category_go_analysis[go_type].update(gene_mappings[go_type][gene])
            
            # Protein families
            if gene in gene_mappings['protein_names']:
                protein_name = gene_mappings['protein_names'][gene].lower()
                family_keywords = ['kinase', 'phosphatase', 'transferase', 'hydrolase', 'oxidoreductase',
                                  'ligase', 'lyase', 'isomerase', 'transporter', 'ribosomal', 'polymerase']
                
                for keyword in family_keywords:
                    if keyword in protein_name:
                        category_protein_families[keyword] += 1
        
        enrichment_analysis[category_name] = {
            'gene_count': len(gene_set),
            'go_bp': category_go_analysis['go_bp'].most_common(10),
            'go_mf': category_go_analysis['go_mf'].most_common(10),
            'go_cc': category_go_analysis['go_cc'].most_common(10),
            'protein_families': category_protein_families.most_common(10)
        }
    
    print(f"Functional enrichment analysis:")
    for category, data in enrichment_analysis.items():
        print(f"  {category}: {data['gene_count']} genes")
    
    return enrichment_analysis

def analyze_functional_similarity_patterns(missing_functional_analysis, metadata):
    """Analyze functional similarity patterns between strains with missing genes."""
    print("\n=== Functional Similarity Pattern Analysis ===")
    
    strain_profiles = missing_functional_analysis['strain_profiles']
    
    # Filter to strains with significant missing genes (≥10)
    significant_strains = {k: v for k, v in strain_profiles.items() if v['missing_count'] >= 10}
    
    if len(significant_strains) < 2:
        print("Not enough strains with significant missing genes for similarity analysis")
        return {}
    
    print(f"Analyzing functional similarity for {len(significant_strains)} strains with ≥10 missing genes")
    
    # Calculate functional similarity between strains
    strain_list = list(significant_strains.keys())
    similarity_data = []
    
    for i, strain1 in enumerate(strain_list):
        for strain2 in strain_list[i+1:]:
            profile1 = significant_strains[strain1]
            profile2 = significant_strains[strain2]
            
            # Calculate GO BP similarity (Jaccard index)
            bp_terms1 = set(profile1['missing_go_bp'].keys())
            bp_terms2 = set(profile2['missing_go_bp'].keys())
            bp_jaccard = len(bp_terms1 & bp_terms2) / len(bp_terms1 | bp_terms2) if len(bp_terms1 | bp_terms2) > 0 else 0
            
            # Calculate protein family similarity
            pf_terms1 = set(profile1['missing_protein_families'].keys())
            pf_terms2 = set(profile2['missing_protein_families'].keys())
            pf_jaccard = len(pf_terms1 & pf_terms2) / len(pf_terms1 | pf_terms2) if len(pf_terms1 | pf_terms2) > 0 else 0
            
            organism1 = get_full_organism_name(strain1, metadata)
            organism2 = get_full_organism_name(strain2, metadata)
            
            similarity_data.append({
                'strain1': strain1,
                'strain2': strain2,
                'organism1': organism1,
                'organism2': organism2,
                'missing_count1': profile1['missing_count'],
                'missing_count2': profile2['missing_count'],
                'bp_similarity': bp_jaccard,
                'protein_family_similarity': pf_jaccard,
                'combined_similarity': (bp_jaccard + pf_jaccard) / 2
            })
    
    similarity_df = pd.DataFrame(similarity_data)
    
    # Statistical analysis
    if len(similarity_df) > 0:
        mean_bp_similarity = similarity_df['bp_similarity'].mean()
        mean_pf_similarity = similarity_df['protein_family_similarity'].mean()
        
        print(f"Functional similarity analysis:")
        print(f"  Mean GO BP similarity: {mean_bp_similarity:.3f}")
        print(f"  Mean protein family similarity: {mean_pf_similarity:.3f}")
        print(f"  Strain pairs analyzed: {len(similarity_df)}")
    
    return {
        'similarity_data': similarity_df,
        'significant_strains': significant_strains
    }

def create_functional_visualizations(functional_composition, missing_functional_analysis, 
                                   enrichment_analysis, similarity_analysis, 
                                   gene_mappings, metadata, output_dir):
    """Create comprehensive functional visualizations."""
    print("\n=== Creating Functional Visualizations ===")
    
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # 1. Core functional composition
    create_core_functional_plots(functional_composition, output_dir)
    
    # 2. Missing genes functional analysis
    create_missing_functional_plots(missing_functional_analysis, metadata, output_dir)
    
    # 3. Functional enrichment analysis
    create_enrichment_plots(enrichment_analysis, output_dir)
    
    # 4. Functional similarity analysis
    if similarity_analysis:
        create_similarity_plots(similarity_analysis, output_dir)
    
    # 5. GO term network analysis
    create_go_network_analysis(functional_composition, missing_functional_analysis, output_dir)

def create_core_functional_plots(functional_composition, output_dir):
    """Create core functional composition plots."""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    # 1. GO Biological Process in core genes
    go_bp_data = functional_composition['go_analysis']['go_bp']['most_common']
    if go_bp_data:
        terms = [item[0][:50] + '...' if len(item[0]) > 50 else item[0] for item in go_bp_data[:15]]
        counts = [item[1] for item in go_bp_data[:15]]
        
        bars1 = ax1.barh(range(len(terms)), counts, color='lightblue', alpha=0.8)
        ax1.set_yticks(range(len(terms)))
        ax1.set_yticklabels(terms, fontsize=9)
        ax1.set_xlabel('Number of Core Genes')
        ax1.set_title('Top GO Biological Processes in Core Genes', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars1, counts):
            width = bar.get_width()
            ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    # 2. GO Molecular Function in core genes
    go_mf_data = functional_composition['go_analysis']['go_mf']['most_common']
    if go_mf_data:
        terms = [item[0][:50] + '...' if len(item[0]) > 50 else item[0] for item in go_mf_data[:15]]
        counts = [item[1] for item in go_mf_data[:15]]
        
        bars2 = ax2.barh(range(len(terms)), counts, color='lightgreen', alpha=0.8)
        ax2.set_yticks(range(len(terms)))
        ax2.set_yticklabels(terms, fontsize=9)
        ax2.set_xlabel('Number of Core Genes')
        ax2.set_title('Top GO Molecular Functions in Core Genes', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars2, counts):
            width = bar.get_width()
            ax2.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    # 3. Protein families in core genes
    protein_families = functional_composition['protein_families']
    if protein_families:
        families = [item[0] for item in protein_families]
        counts = [item[1] for item in protein_families]
        
        bars3 = ax3.bar(range(len(families)), counts, color='lightsalmon', alpha=0.8)
        ax3.set_xticks(range(len(families)))
        ax3.set_xticklabels(families, rotation=45, ha='right', fontsize=9)
        ax3.set_ylabel('Number of Core Genes')
        ax3.set_title('Core Gene Protein Families', fontsize=12, fontweight='bold')
        ax3.grid(True, alpha=0.3, axis='y')
        
        for bar, count in zip(bars3, counts):
            height = bar.get_height()
            ax3.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{count}', ha='center', va='bottom', fontsize=8)
    
    # 4. EC number classes in core genes
    ec_classes = functional_composition['ec_classes']
    if ec_classes:
        # EC class descriptions
        ec_descriptions = {
            '1': 'Oxidoreductases',
            '2': 'Transferases', 
            '3': 'Hydrolases',
            '4': 'Lyases',
            '5': 'Isomerases',
            '6': 'Ligases'
        }
        
        classes = [f"EC {item[0]}: {ec_descriptions.get(item[0], 'Unknown')}" for item in ec_classes]
        counts = [item[1] for item in ec_classes]
        
        colors = plt.cm.Set3(np.linspace(0, 1, len(classes)))
        wedges, texts, autotexts = ax4.pie(counts, labels=classes, autopct='%1.1f%%', 
                                          colors=colors, startangle=90)
        ax4.set_title('EC Number Classes in Core Genes', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'core_functional_composition.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_missing_functional_plots(missing_functional_analysis, metadata, output_dir):
    """Create missing genes functional analysis plots."""
    
    aggregated_functions = missing_functional_analysis['aggregated_functions']
    strain_profiles = missing_functional_analysis['strain_profiles']
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    # 1. Most frequently missing GO BP terms
    missing_bp = aggregated_functions['go_bp'].most_common(15)
    if missing_bp:
        terms = [item[0][:45] + '...' if len(item[0]) > 45 else item[0] for item in missing_bp]
        counts = [item[1] for item in missing_bp]
        
        bars1 = ax1.barh(range(len(terms)), counts, color='lightcoral', alpha=0.8)
        ax1.set_yticks(range(len(terms)))
        ax1.set_yticklabels(terms, fontsize=9)
        ax1.set_xlabel('Number of Strains Missing This Function')
        ax1.set_title('Most Frequently Missing GO Biological Processes', fontsize=12, fontweight='bold')
        ax1.grid(True, alpha=0.3, axis='x')
        
        for bar, count in zip(bars1, counts):
            width = bar.get_width()
            ax1.text(width + 0.1, bar.get_y() + bar.get_height()/2,
                    f'{count}', ha='left', va='center', fontsize=8)
    
    # 2. Most frequently missing protein families
    missing_families = aggregated_functions['protein_families'].most_common(12)
    if missing_families:
        families = [item[0] for item in missing_families]
        counts = [item[1] for item in missing_families]
        
        bars2 = ax2.bar(range(len(families)), counts, color='lightsteelblue', alpha=0.8)
        ax2.set_xticks(range(len(families)))
        ax2.set_xticklabels(families, rotation=45, ha='right', fontsize=9)
        ax2.set_ylabel('Number of Strains Missing This Family')
        ax2.set_title('Most Frequently Missing Protein Families', fontsize=12, fontweight='bold')
        ax2.grid(True, alpha=0.3, axis='y')
        
        for bar, count in zip(bars2, counts):
            height = bar.get_height()
            ax2.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                    f'{count}', ha='center', va='bottom', fontsize=8)
    
    # 3. Strains with most missing functional categories
    strain_functional_counts = []
    for strain_id, profile in strain_profiles.items():
        total_missing_functions = (len(profile['missing_go_bp']) + 
                                 len(profile['missing_go_mf']) + 
                                 len(profile['missing_protein_families']))
        strain_functional_counts.append((strain_id, profile['organism'], 
                                       profile['missing_count'], total_missing_functions))
    
    # Sort by total missing functions and take top 15
    strain_functional_counts.sort(key=lambda x: x[3], reverse=True)
    top_strains = strain_functional_counts[:15]
    
    if top_strains:
        organism_names = []
        functional_counts = []
        missing_gene_counts = []
        
        for strain_id, organism, missing_genes, missing_functions in top_strains:
            if len(organism) > 45:
                organism_names.append(organism[:42] + '...')
            else:
                organism_names.append(organism)
            functional_counts.append(missing_functions)
            missing_gene_counts.append(missing_genes)
        
        x_pos = np.arange(len(organism_names))
        width = 0.35
        
        bars3a = ax3.bar(x_pos - width/2, missing_gene_counts, width, label='Missing Genes', 
                        color='orange', alpha=0.7)
        bars3b = ax3.bar(x_pos + width/2, functional_counts, width, label='Missing Functions', 
                        color='purple', alpha=0.7)
        
        ax3.set_xticks(x_pos)
        ax3.set_xticklabels(organism_names, rotation=45, ha='right', fontsize=8)
        ax3.set_ylabel('Count')
        ax3.set_title('Strains with Most Missing Functional Categories', fontsize=12, fontweight='bold')
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Missing EC classes distribution
    missing_ec = aggregated_functions['ec_classes'].most_common()
    if missing_ec:
        ec_descriptions = {
            '1': 'Oxidoreductases',
            '2': 'Transferases', 
            '3': 'Hydrolases',
            '4': 'Lyases',
            '5': 'Isomerases',
            '6': 'Ligases'
        }
        
        classes = [f"EC {item[0]}: {ec_descriptions.get(item[0], 'Unknown')}" for item in missing_ec]
        counts = [item[1] for item in missing_ec]
        
        colors = plt.cm.Reds(np.linspace(0.3, 0.9, len(classes)))
        wedges, texts, autotexts = ax4.pie(counts, labels=classes, autopct='%1.1f%%', 
                                          colors=colors, startangle=90)
        ax4.set_title('Missing EC Number Classes Distribution', fontsize=12, fontweight='bold')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'missing_functional_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_enrichment_plots(enrichment_analysis, output_dir):
    """Create functional enrichment comparison plots."""
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(18, 14))
    
    categories = ['Core', 'Accessory', 'Rare', 'Unique']
    colors = ['#66b3ff', '#99ff99', '#ffcc99', '#ff9999']
    
    # 1. Gene count comparison
    gene_counts = [enrichment_analysis[cat]['gene_count'] for cat in categories]
    bars1 = ax1.bar(categories, gene_counts, color=colors, alpha=0.8)
    ax1.set_ylabel('Number of Genes')
    ax1.set_title('Gene Count by Category', fontsize=12, fontweight='bold')
    ax1.grid(True, alpha=0.3, axis='y')
    
    for bar, count in zip(bars1, gene_counts):
        height = bar.get_height()
        ax1.text(bar.get_x() + bar.get_width()/2., height + 0.01*max(gene_counts),
                f'{count:,}', ha='center', va='bottom', fontsize=10, fontweight='bold')
    
    # 2. Top protein families comparison
    family_comparison = {}
    for category in categories:
        families = enrichment_analysis[category]['protein_families']
        for family, count in families:
            if family not in family_comparison:
                family_comparison[family] = {}
            family_comparison[family][category] = count
    
    # Get top families across all categories
    top_families = sorted(family_comparison.keys(), 
                         key=lambda x: sum(family_comparison[x].values()), reverse=True)[:8]
    
    x = np.arange(len(top_families))
    width = 0.2
    
    for i, category in enumerate(categories):
        counts = [family_comparison[family].get(category, 0) for family in top_families]
        ax2.bar(x + i*width, counts, width, label=category, color=colors[i], alpha=0.8)
    
    ax2.set_xlabel('Protein Families')
    ax2.set_ylabel('Number of Genes')
    ax2.set_title('Protein Family Distribution by Gene Category', fontsize=12, fontweight='bold')
    ax2.set_xticks(x + width * 1.5)
    ax2.set_xticklabels(top_families, rotation=45, ha='right')
    ax2.legend()
    ax2.grid(True, alpha=0.3, axis='y')
    
    # 3. GO BP term comparison (top terms from core)
    core_bp_terms = [item[0] for item in enrichment_analysis['Core']['go_bp'][:6]]
    
    if core_bp_terms:
        bp_comparison = {}
        for term in core_bp_terms:
            bp_comparison[term] = {}
            for category in categories:
                bp_data = dict(enrichment_analysis[category]['go_bp'])
                bp_comparison[term][category] = bp_data.get(term, 0)
        
        x = np.arange(len(core_bp_terms))
        width = 0.2
        
        for i, category in enumerate(categories):
            counts = [bp_comparison[term].get(category, 0) for term in core_bp_terms]
            ax3.bar(x + i*width, counts, width, label=category, color=colors[i], alpha=0.8)
        
        ax3.set_xlabel('GO Biological Process Terms')
        ax3.set_ylabel('Number of Genes')
        ax3.set_title('Key GO BP Terms Distribution by Category', fontsize=12, fontweight='bold')
        ax3.set_xticks(x + width * 1.5)
        
        # Truncate labels for readability
        truncated_terms = [term[:20] + '...' if len(term) > 20 else term for term in core_bp_terms]
        ax3.set_xticklabels(truncated_terms, rotation=45, ha='right', fontsize=9)
        ax3.legend()
        ax3.grid(True, alpha=0.3, axis='y')
    
    # 4. Functional diversity comparison
    diversity_metrics = []
    for category in categories:
        bp_diversity = len(enrichment_analysis[category]['go_bp'])
        mf_diversity = len(enrichment_analysis[category]['go_mf'])
        pf_diversity = len(enrichment_analysis[category]['protein_families'])
        
        diversity_metrics.append({
            'category': category,
            'bp_diversity': bp_diversity,
            'mf_diversity': mf_diversity,
            'pf_diversity': pf_diversity
        })
    
    df_diversity = pd.DataFrame(diversity_metrics)
    
    x = np.arange(len(categories))
    width = 0.25
    
    ax4.bar(x - width, df_diversity['bp_diversity'], width, label='GO BP Terms', 
           color='lightblue', alpha=0.8)
    ax4.bar(x, df_diversity['mf_diversity'], width, label='GO MF Terms', 
           color='lightgreen', alpha=0.8)
    ax4.bar(x + width, df_diversity['pf_diversity'], width, label='Protein Families', 
           color='lightsalmon', alpha=0.8)
    
    ax4.set_xlabel('Gene Categories')
    ax4.set_ylabel('Number of Unique Functional Terms')
    ax4.set_title('Functional Diversity by Gene Category', fontsize=12, fontweight='bold')
    ax4.set_xticks(x)
    ax4.set_xticklabels(categories)
    ax4.legend()
    ax4.grid(True, alpha=0.3, axis='y')
    
    plt.tight_layout()
    plt.savefig(output_dir / 'functional_enrichment_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_similarity_plots(similarity_analysis, output_dir):
    """Create functional similarity analysis plots."""
    
    similarity_df = similarity_analysis['similarity_data']
    significant_strains = similarity_analysis['significant_strains']
    
    if similarity_df.empty:
        return
    
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(16, 12))
    
    # 1. Functional similarity scatter plot
    ax1.scatter(similarity_df['bp_similarity'], similarity_df['protein_family_similarity'], 
               alpha=0.6, s=60, color='purple')
    ax1.set_xlabel('GO Biological Process Similarity')
    ax1.set_ylabel('Protein Family Similarity')
    ax1.set_title('Functional Similarity Between Strain Pairs')
    ax1.grid(True, alpha=0.3)
    
    # Add correlation
    corr = similarity_df['bp_similarity'].corr(similarity_df['protein_family_similarity'])
    ax1.text(0.05, 0.95, f'r = {corr:.3f}', transform=ax1.transAxes,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightblue', alpha=0.8))
    
    # 2. Distribution of similarities
    ax2.hist(similarity_df['combined_similarity'], bins=20, alpha=0.7, 
            color='lightcoral', edgecolor='black')
    ax2.axvline(similarity_df['combined_similarity'].mean(), color='red', 
               linestyle='--', label=f'Mean: {similarity_df["combined_similarity"].mean():.3f}')
    ax2.set_xlabel('Combined Functional Similarity')
    ax2.set_ylabel('Number of Strain Pairs')
    ax2.set_title('Distribution of Functional Similarities')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    
    # 3. Missing gene count vs similarity
    missing_diff = abs(similarity_df['missing_count1'] - similarity_df['missing_count2'])
    ax3.scatter(missing_diff, similarity_df['combined_similarity'], 
               alpha=0.6, s=60, color='green')
    ax3.set_xlabel('Difference in Missing Gene Counts')
    ax3.set_ylabel('Combined Functional Similarity')
    ax3.set_title('Missing Gene Count Difference vs Functional Similarity')
    ax3.grid(True, alpha=0.3)
    
    # Add correlation
    corr2 = missing_diff.corr(similarity_df['combined_similarity'])
    ax3.text(0.05, 0.95, f'r = {corr2:.3f}', transform=ax3.transAxes,
            bbox=dict(boxstyle='round,pad=0.3', facecolor='lightgreen', alpha=0.8))
    
    # 4. Top similar strain pairs
    top_similar = similarity_df.nlargest(10, 'combined_similarity')
    
    pair_labels = []
    similarities = []
    for _, row in top_similar.iterrows():
        # Truncate organism names for display
        org1 = row['organism1'][:20] + '...' if len(row['organism1']) > 20 else row['organism1']
        org2 = row['organism2'][:20] + '...' if len(row['organism2']) > 20 else row['organism2']
        pair_labels.append(f"{org1}\nvs\n{org2}")
        similarities.append(row['combined_similarity'])
    
    bars = ax4.barh(range(len(pair_labels)), similarities, color='skyblue', alpha=0.8)
    ax4.set_yticks(range(len(pair_labels)))
    ax4.set_yticklabels(pair_labels, fontsize=8)
    ax4.set_xlabel('Combined Functional Similarity')
    ax4.set_title('Top 10 Most Functionally Similar Strain Pairs')
    ax4.grid(True, alpha=0.3, axis='x')
    
    for bar, sim in zip(bars, similarities):
        width = bar.get_width()
        ax4.text(width + 0.01, bar.get_y() + bar.get_height()/2,
                f'{sim:.3f}', ha='left', va='center', fontsize=8)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'functional_similarity_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def create_go_network_analysis(functional_composition, missing_functional_analysis, output_dir):
    """Create GO term network analysis visualization."""
    
    # Extract core GO terms and missing GO terms
    core_go_terms = functional_composition['go_analysis']['go_bp']['most_common'][:20]
    missing_go_terms = missing_functional_analysis['aggregated_functions']['go_bp'].most_common(15)
    
    if not core_go_terms and not missing_go_terms:
        return
    
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))
    
    # 1. Core vs Missing GO terms comparison
    core_terms_dict = dict(core_go_terms)
    missing_terms_dict = dict(missing_go_terms)
    
    # Find overlapping terms
    all_terms = set(core_terms_dict.keys()) | set(missing_terms_dict.keys())
    comparison_data = []
    
    for term in all_terms:
        core_count = core_terms_dict.get(term, 0)
        missing_count = missing_terms_dict.get(term, 0)
        comparison_data.append({
            'term': term[:40] + '...' if len(term) > 40 else term,
            'core_count': core_count,
            'missing_count': missing_count,
            'total': core_count + missing_count
        })
    
    # Sort by total and take top 15
    comparison_data.sort(key=lambda x: x['total'], reverse=True)
    top_comparison = comparison_data[:15]
    
    if top_comparison:
        terms = [item['term'] for item in top_comparison]
        core_counts = [item['core_count'] for item in top_comparison]
        missing_counts = [item['missing_count'] for item in top_comparison]
        
        y_pos = np.arange(len(terms))
        width = 0.35
        
        bars1 = ax1.barh(y_pos - width/2, core_counts, width, label='Core Genes', 
                        color='lightblue', alpha=0.8)
        bars2 = ax1.barh(y_pos + width/2, missing_counts, width, label='Missing Genes', 
                        color='lightcoral', alpha=0.8)
        
        ax1.set_yticks(y_pos)
        ax1.set_yticklabels(terms, fontsize=9)
        ax1.set_xlabel('Number of Genes')
        ax1.set_title('GO Biological Process: Core vs Missing Genes', fontsize=12, fontweight='bold')
        ax1.legend()
        ax1.grid(True, alpha=0.3, axis='x')
    
    # 2. GO term frequency heatmap (simplified)
    # Create a matrix of GO terms vs functional categories
    strain_profiles = missing_functional_analysis['strain_profiles']
    
    # Get top GO terms from missing analysis
    top_missing_terms = [item[0] for item in missing_go_terms[:10]]
    
    if top_missing_terms and strain_profiles:
        # Create matrix: strains vs GO terms
        strain_go_matrix = []
        strain_names = []
        
        # Take a subset of strains for visualization
        strain_subset = list(strain_profiles.keys())[:15]  # Top 15 strains
        
        for strain_id in strain_subset:
            profile = strain_profiles[strain_id]
            strain_names.append(profile['organism'][:25] + '...' if len(profile['organism']) > 25 else profile['organism'])
            
            strain_go_counts = []
            for term in top_missing_terms:
                count = profile['missing_go_bp'].get(term, 0)
                strain_go_counts.append(count)
            strain_go_matrix.append(strain_go_counts)
        
        strain_go_matrix = np.array(strain_go_matrix)
        
        if strain_go_matrix.size > 0:
            # Truncate GO terms for labels
            term_labels = [term[:25] + '...' if len(term) > 25 else term for term in top_missing_terms]
            
            im = ax2.imshow(strain_go_matrix, cmap='Reds', aspect='auto')
            ax2.set_xticks(range(len(term_labels)))
            ax2.set_xticklabels(term_labels, rotation=45, ha='right', fontsize=8)
            ax2.set_yticks(range(len(strain_names)))
            ax2.set_yticklabels(strain_names, fontsize=8)
            ax2.set_title('Missing GO Terms Heatmap\n(Strains vs GO Biological Processes)', 
                         fontsize=12, fontweight='bold')
            
            # Add colorbar
            cbar = plt.colorbar(im, ax=ax2, shrink=0.6)
            cbar.set_label('Number of Missing Genes', rotation=270, labelpad=15)
    
    plt.tight_layout()
    plt.savefig(output_dir / 'go_network_analysis.png', dpi=300, bbox_inches='tight')
    plt.close()

def save_functional_results(functional_composition, missing_functional_analysis, 
                          enrichment_analysis, similarity_analysis, gene_mappings, 
                          metadata, output_dir):
    """Save all functional analysis results."""
    print("\n=== Saving Functional Analysis Results ===")
    
    output_dir = Path(output_dir)
    
    # 1. Core functional composition
    core_go_bp = pd.DataFrame(functional_composition['go_analysis']['go_bp']['most_common'],
                             columns=['GO_BP_Term', 'Gene_Count'])
    core_go_bp.to_csv(output_dir / 'core_go_biological_process.tsv', sep='\t', index=False)
    
    core_go_mf = pd.DataFrame(functional_composition['go_analysis']['go_mf']['most_common'],
                             columns=['GO_MF_Term', 'Gene_Count'])
    core_go_mf.to_csv(output_dir / 'core_go_molecular_function.tsv', sep='\t', index=False)
    
    core_families = pd.DataFrame(functional_composition['protein_families'],
                                columns=['Protein_Family', 'Gene_Count'])
    core_families.to_csv(output_dir / 'core_protein_families.tsv', sep='\t', index=False)
    
    # 2. Missing functional analysis
    strain_missing_data = []
    for strain_id, profile in missing_functional_analysis['strain_profiles'].items():
        strain_missing_data.append({
            'strain_id': strain_id,
            'organism': profile['organism'],
            'missing_gene_count': profile['missing_count'],
            'missing_bp_terms': len(profile['missing_go_bp']),
            'missing_mf_terms': len(profile['missing_go_mf']),
            'missing_protein_families': len(profile['missing_protein_families']),
            'top_missing_bp': '; '.join([f"{term}({count})" for term, count in profile['missing_go_bp'].most_common(5)]),
            'top_missing_families': '; '.join([f"{family}({count})" for family, count in profile['missing_protein_families'].most_common(3)])
        })
    
    df_strain_missing = pd.DataFrame(strain_missing_data)
    df_strain_missing = df_strain_missing.sort_values('missing_gene_count', ascending=False)
    df_strain_missing.to_csv(output_dir / 'strain_missing_functional_analysis.tsv', sep='\t', index=False)
    
    # 3. Aggregated missing functions
    missing_bp = pd.DataFrame(missing_functional_analysis['aggregated_functions']['go_bp'].most_common(),
                             columns=['Missing_GO_BP_Term', 'Strain_Count'])
    missing_bp.to_csv(output_dir / 'frequently_missing_go_bp.tsv', sep='\t', index=False)
    
    missing_families = pd.DataFrame(missing_functional_analysis['aggregated_functions']['protein_families'].most_common(),
                                   columns=['Missing_Protein_Family', 'Strain_Count'])
    missing_families.to_csv(output_dir / 'frequently_missing_protein_families.tsv', sep='\t', index=False)
    
    # 4. Enrichment analysis by category
    for category, data in enrichment_analysis.items():
        category_go_bp = pd.DataFrame(data['go_bp'], columns=['GO_BP_Term', 'Gene_Count'])
        category_go_bp.to_csv(output_dir / f'{category.lower()}_genes_go_bp.tsv', sep='\t', index=False)
        
        category_families = pd.DataFrame(data['protein_families'], columns=['Protein_Family', 'Gene_Count'])
        category_families.to_csv(output_dir / f'{category.lower()}_genes_protein_families.tsv', sep='\t', index=False)
    
    # 5. Functional similarity analysis
    if similarity_analysis and 'similarity_data' in similarity_analysis:
        similarity_df = similarity_analysis['similarity_data']
        similarity_df.to_csv(output_dir / 'functional_similarity_analysis.tsv', sep='\t', index=False)
    
    # 6. Comprehensive summary
    summary_lines = [
        "FUNCTIONAL CORE GENOME ANALYSIS SUMMARY",
        "=" * 50,
        f"Analysis Date: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}",
        "",
        "CORE FUNCTIONAL COMPOSITION:",
        f"  GO Biological Process terms: {functional_composition['go_analysis']['go_bp']['total_terms']}",
        f"  GO Molecular Function terms: {functional_composition['go_analysis']['go_mf']['total_terms']}",
        f"  GO Cellular Component terms: {functional_composition['go_analysis']['go_cc']['total_terms']}",
        f"  Protein family keywords identified: {len(functional_composition['protein_families'])}",
        f"  EC number classes found: {len(functional_composition['ec_classes'])}",
        "",
        "MISSING FUNCTIONAL ANALYSIS:",
        f"  Strains with missing functional annotations: {len(missing_functional_analysis['strain_profiles'])}",
        f"  Unique missing GO BP terms: {len(missing_functional_analysis['aggregated_functions']['go_bp'])}",
        f"  Unique missing protein families: {len(missing_functional_analysis['aggregated_functions']['protein_families'])}",
        "",
        "FUNCTIONAL ENRICHMENT:",
    ]
    
    for category, data in enrichment_analysis.items():
        summary_lines.append(f"  {category} genes: {data['gene_count']:,} genes")
        summary_lines.append(f"    GO BP terms: {len(data['go_bp'])}")
        summary_lines.append(f"    Protein families: {len(data['protein_families'])}")
    
    if similarity_analysis and 'similarity_data' in similarity_analysis:
        similarity_df = similarity_analysis['similarity_data']
        summary_lines.extend([
            "",
            "FUNCTIONAL SIMILARITY ANALYSIS:",
            f"  Strain pairs analyzed: {len(similarity_df)}",
            f"  Mean GO BP similarity: {similarity_df['bp_similarity'].mean():.3f}",
            f"  Mean protein family similarity: {similarity_df['protein_family_similarity'].mean():.3f}",
        ])
    
    summary_lines.extend([
        "",
        "TOP CORE FUNCTIONAL CATEGORIES:",
    ])
    
    # Add top categories
    for term, count in functional_composition['go_analysis']['go_bp']['most_common'][:5]:
        summary_lines.append(f"  {term}: {count} genes")
    
    summary_lines.extend([
        "",
        "TOP MISSING FUNCTIONAL CATEGORIES:",
    ])
    
    for term, count in missing_functional_analysis['aggregated_functions']['go_bp'].most_common(5):
        summary_lines.append(f"  {term}: missing in {count} strains")
    
    summary_lines.extend([
        "",
        "FILES GENERATED:",
        "  • core_go_*.tsv - Core gene GO term analysis",
        "  • *_genes_*.tsv - Enrichment analysis by gene category",
        "  • strain_missing_functional_analysis.tsv - Per-strain missing functions",
        "  • frequently_missing_*.tsv - Most commonly missing functions",
        "  • functional_similarity_analysis.tsv - Strain functional similarity",
        "  • Multiple functional visualization PNG files"
    ])
    
    with open(output_dir / 'functional_analysis_summary.txt', 'w') as f:
        f.write('\n'.join(summary_lines))
    
    print(f"Functional analysis results saved to: {output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Functional Core Genome Analysis')
    parser.add_argument('--gene-matrix', required=True, help='Path to gene matrix NPZ file')
    parser.add_argument('--gene-labels', required=True, help='Path to gene labels file')
    parser.add_argument('--core-genes', required=True, help='Path to core genes list file')
    parser.add_argument('--annotations', required=True, help='Path to gene annotations TSV file')
    parser.add_argument('--metadata', help='Path to proteome metadata TSV file')
    parser.add_argument('--output-dir', required=True, help='Output directory')
    
    args = parser.parse_args()
    
    print("🧬 FUNCTIONAL CORE GENOME ANALYSIS")
    print("=" * 45)
    
    try:
        # Load data
        df_genes, core_genes, annotations, metadata = load_data(
            args.gene_matrix, args.gene_labels, args.core_genes, 
            args.annotations, args.metadata
        )
        
        print(f"Analyzing functional patterns in {df_genes.shape[0]} genes across {df_genes.shape[1]} strains")
        
        # Create functional mappings
        gene_mappings = create_functional_mappings(annotations)
        
        # Run functional analyses
        functional_composition = analyze_core_functional_composition(df_genes, core_genes, gene_mappings)
        missing_functional_analysis = analyze_missing_genes_functional(df_genes, core_genes, gene_mappings, metadata)
        enrichment_analysis = analyze_functional_enrichment(df_genes, core_genes, gene_mappings)
        similarity_analysis = analyze_functional_similarity_patterns(missing_functional_analysis, metadata)
        
        # Create output directory
        output_dir = Path(args.output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate visualizations
        create_functional_visualizations(
            functional_composition, missing_functional_analysis, enrichment_analysis, 
            similarity_analysis, gene_mappings, metadata, output_dir
        )
        
        # Save results
        save_functional_results(
            functional_composition, missing_functional_analysis, enrichment_analysis, 
            similarity_analysis, gene_mappings, metadata, output_dir
        )
        
        print(f"\nFunctional analysis completed!")
        print(f"Results saved to {output_dir}")
        
        # Print key findings
        print(f"\nKey functional findings:")
        print(f"  • Core GO BP terms: {functional_composition['go_analysis']['go_bp']['total_terms']}")
        print(f"  • Core protein families: {len(functional_composition['protein_families'])}")
        print(f"  • Strains with missing functions: {len(missing_functional_analysis['strain_profiles'])}")
        
        if similarity_analysis and 'similarity_data' in similarity_analysis:
            similarity_df = similarity_analysis['similarity_data']
            print(f"  • Functional similarity pairs: {len(similarity_df)}")
            print(f"  • Mean functional similarity: {similarity_df['combined_similarity'].mean():.3f}")
        
        # Top functional categories
        top_core_function = functional_composition['go_analysis']['go_bp']['most_common'][0]
        print(f"  • Top core function: {top_core_function[0]} ({top_core_function[1]} genes)")
        
        if missing_functional_analysis['aggregated_functions']['go_bp']:
            top_missing_function = missing_functional_analysis['aggregated_functions']['go_bp'].most_common(1)[0]
            print(f"  • Top missing function: {top_missing_function[0]} (missing in {top_missing_function[1]} strains)")
        
    except Exception as e:
        print(f"ERROR: {str(e)}")
        import traceback
        traceback.print_exc()
        sys.exit(1)

if __name__ == '__main__':
    main()
