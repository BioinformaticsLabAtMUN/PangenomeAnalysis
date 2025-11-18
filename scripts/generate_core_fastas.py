#!/usr/bin/env python3
import os
import glob
import argparse
import pandas as pd
from collections import defaultdict, Counter
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

def parse_allele_mappings(allele_names_file):
    """
    Parse the enhanced allele names TSV file to create bidirectional mappings.
    
    Returns:
        original_to_cluster: dict mapping original ID -> cluster ID
        original_to_allele: dict mapping original ID -> allele ID
        cluster_to_originals: dict mapping cluster ID -> list of original IDs
    """
    print(f"Parsing allele mappings from {allele_names_file}...")
    
    # Check if the file has a header
    with open(allele_names_file, 'r') as f:
        first_line = f.readline().strip()
        has_header = 'Original_ID' in first_line or 'New_ID' in first_line
    
    # Read the TSV file
    if has_header:
        df = pd.read_csv(allele_names_file, sep='\t')
        # Handle both old and new format
        if 'Original_ID' in df.columns and 'Cluster_ID' in df.columns:
            # New format
            original_to_cluster = dict(zip(df['Original_ID'], df['Cluster_ID']))
            original_to_allele = dict(zip(df['Original_ID'], df['New_ID']))
        else:
            # Old format - first column is allele ID, others are originals
            original_to_cluster = {}
            original_to_allele = {}
            for _, row in df.iterrows():
                allele_id = row.iloc[0]
                cluster_id = allele_id.rsplit('A', 1)[0]  # Remove allele suffix
                for original_id in row.iloc[1:]:
                    if pd.notna(original_id):
                        original_to_cluster[original_id] = cluster_id
                        original_to_allele[original_id] = allele_id
    else:
        # Handle the old format without header
        original_to_cluster = {}
        original_to_allele = {}
        with open(allele_names_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) < 2:
                    continue
                allele_id = fields[0]
                cluster_id = allele_id.rsplit('A', 1)[0]  # Remove allele suffix
                for original_id in fields[1:]:
                    if original_id:
                        original_to_cluster[original_id] = cluster_id
                        original_to_allele[original_id] = allele_id
    
    # Create reverse mapping
    cluster_to_originals = defaultdict(list)
    for original_id, cluster_id in original_to_cluster.items():
        cluster_to_originals[cluster_id].append(original_id)
    
    print(f"Found {len(original_to_cluster)} original IDs mapping to {len(cluster_to_originals)} clusters")
    return original_to_cluster, original_to_allele, cluster_to_originals

def parse_core_genes(core_genes_file):
    """Read the list of core gene cluster IDs."""
    print(f"Reading core genes from {core_genes_file}...")
    with open(core_genes_file) as f:
        core_genes = {line.strip() for line in f if line.strip()}
    print(f"Found {len(core_genes)} core genes")
    return core_genes

def extract_genome_id(filename):
    """Extract genome ID from filename, handling various formats."""
    # Extract just the filename without path
    basename = os.path.basename(filename)
    # Remove .faa extension
    basename = basename.replace('.faa', '')
    
    # Common genome ID patterns
    if basename.startswith('GCF_') or basename.startswith('GCA_'):
        return basename
    
    # For other formats, use the entire basename
    return basename

def collect_genome_sequences(input_dir, original_to_cluster, original_to_allele, core_genes):
    """Scan all FAA files recursively and collect sequences for core genes."""
    cluster_sequences = defaultdict(list)
    genome_stats = defaultdict(int)
    print(f"\nScanning for FAA files recursively in {input_dir}...")
    
    # Recursively find all .faa files
    faa_files = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if file.endswith('.faa') and not any(x in file for x in ['renamed', 'cdhit', 'nr', 'swift']):
                faa_files.append(os.path.join(root, file))
    
    print(f"Found {len(faa_files)} FAA files for genome analysis")
    
    # Process each FAA file
    for faa in faa_files:
        genome_id = extract_genome_id(faa)
        print(f"Processing {os.path.basename(faa)} (Genome: {genome_id})")
        
        genes_found = set()
        for record in SeqIO.parse(faa, "fasta"):
            protein_id = record.id
            if protein_id in original_to_cluster:
                cluster_id = original_to_cluster[protein_id]
                allele_id = original_to_allele.get(protein_id, f"{cluster_id}A?")
                
                # Only process core genes
                if cluster_id in core_genes:
                    # Create a better FASTA header
                    display_id = f"{genome_id}|{allele_id}"
                    description = f"original_id={protein_id} cluster={cluster_id}"
                    
                    cluster_sequences[cluster_id].append((
                        genome_id, 
                        allele_id, 
                        protein_id, 
                        record.seq,
                        display_id,
                        description
                    ))
                    genes_found.add(cluster_id)
        
        # Update stats for this genome
        genome_stats[genome_id] = len(genes_found)
        
    print(f"Collected sequences from {len(genome_stats)} genomes")
    return cluster_sequences, genome_stats

def generate_core_fastas(core_genes, output_dir, cluster_sequences, genome_stats, prefix="", method=""):
    """Generate multi-FASTA files for each core gene cluster and output summary statistics."""
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"\nGenerating core gene FASTA files using {method} clustering...")
    summary_path = os.path.join(output_dir, f"{prefix}core_genes_summary.txt")
    
    # Prepare detailed CSV with better formatting
    with open(summary_path, 'w') as summary_f:
        summary_f.write("Cluster_ID\tNum_Genomes\tTotal_Sequences\tGenome_Coverage\tGenomes\n")
        
        unique_genome_counts = []
        total_genomes = len(genome_stats)
        
        for cluster_id in sorted(core_genes):
            sequences = cluster_sequences.get(cluster_id, [])
            if sequences:
                # Create a FASTA file for each core gene
                output_file = os.path.join(output_dir, f"{prefix}{cluster_id}.faa")
                records = []
                
                # Track unique genomes for this cluster
                genomes = set()
                
                for genome_id, allele_id, protein_id, sequence, display_id, description in sequences:
                    record = SeqRecord(
                        sequence,
                        id=display_id,
                        description=description
                    )
                    records.append(record)
                    genomes.add(genome_id)
                
                # Write the FASTA file
                SeqIO.write(records, output_file, "fasta")
                
                # Update statistics
                num_genomes = len(genomes)
                unique_genome_counts.append(num_genomes)
                genome_coverage = round(num_genomes / total_genomes * 100, 2) if total_genomes > 0 else 0
                
                # Write to summary
                summary_f.write(f"{cluster_id}\t{num_genomes}\t{len(sequences)}\t{genome_coverage}%\t{','.join(sorted(genomes))}\n")
            else:
                print(f"Warning: No sequences found for core gene cluster {cluster_id}")
    
    # Generate genome coverage report
    coverage_path = os.path.join(output_dir, f"{prefix}genome_coverage.tsv")
    with open(coverage_path, 'w') as f:
        f.write("Genome_ID\tCore_Genes_Found\tCoverage_Percent\n")
        for genome_id, count in sorted(genome_stats.items()):
            coverage = round(count / len(core_genes) * 100, 2) if core_genes else 0
            f.write(f"{genome_id}\t{count}\t{coverage}%\n")
    
    # Generate statistics and plots
    if unique_genome_counts:
        generate_statistics(unique_genome_counts, output_dir, prefix, total_genomes)
    else:
        print("No core gene clusters processed.")

def generate_statistics(unique_genome_counts, output_dir, prefix, total_genomes):
    """Generate statistical analysis and plots for core gene distribution."""
    import numpy as np
    
    # Calculate basic statistics
    total_clusters = len(unique_genome_counts)
    mean_genomes = np.mean(unique_genome_counts)
    median_genomes = np.median(unique_genome_counts)
    min_genomes = np.min(unique_genome_counts)
    max_genomes = np.max(unique_genome_counts)
    
    # Print summary statistics
    print("\n=== Core Gene Cluster Statistics ===")
    print(f"Total clusters: {total_clusters}")
    print(f"Total genomes analyzed: {total_genomes}")
    print(f"Average unique genomes per cluster: {mean_genomes:.2f}")
    print(f"Median unique genomes per cluster: {median_genomes}")
    print(f"Minimum unique genomes in a cluster: {min_genomes}")
    print(f"Maximum unique genomes in a cluster: {max_genomes}")
    
    # Save statistics to file
    stats_path = os.path.join(output_dir, f"{prefix}statistics.txt")
    with open(stats_path, 'w') as f:
        f.write("=== Core Gene Cluster Statistics ===\n")
        f.write(f"Total clusters: {total_clusters}\n")
        f.write(f"Total genomes analyzed: {total_genomes}\n")
        f.write(f"Average unique genomes per cluster: {mean_genomes:.2f}\n")
        f.write(f"Median unique genomes per cluster: {median_genomes}\n")
        f.write(f"Minimum unique genomes in a cluster: {min_genomes}\n")
        f.write(f"Maximum unique genomes in a cluster: {max_genomes}\n\n")
        
        # Distribution details
        f.write("Distribution (Unique Genome Count: Frequency):\n")
        distribution = Counter(unique_genome_counts)
        for genome_count, freq in sorted(distribution.items()):
            percent = round(freq / total_clusters * 100, 2)
            f.write(f"{genome_count}: {freq} ({percent}%)\n")
    
    # Generate histogram
    try:
        import matplotlib.pyplot as plt
        plt.figure(figsize=(10, 6))
        
        # Create histogram
        plt.hist(unique_genome_counts, bins=range(min_genomes, max_genomes + 1), 
                 edgecolor='black', alpha=0.7, color='steelblue')
        
        # Add title and labels
        plt.xlabel("Number of Unique Genomes", fontsize=12)
        plt.ylabel("Number of Core Gene Clusters", fontsize=12)
        plt.title(f"Distribution of Unique Genomes per Core Gene Cluster\n(Total: {total_clusters} clusters across {total_genomes} genomes)", 
                 fontsize=14)
        
        # Add grid for better readability
        plt.grid(axis='y', alpha=0.3)
        
        # Add average line
        plt.axvline(x=mean_genomes, color='red', linestyle='--', 
                   label=f'Mean: {mean_genomes:.2f}')
        plt.axvline(x=median_genomes, color='green', linestyle='-', 
                   label=f'Median: {median_genomes}')
        
        plt.legend()
        plt.tight_layout()
        
        # Save the plot
        hist_path = os.path.join(output_dir, f"{prefix}unique_genomes_distribution.png")
        plt.savefig(hist_path, dpi=300)
        plt.close()
        print(f"Histogram saved to {hist_path}")
        
    except ImportError:
        print("matplotlib or numpy is not installed; skipping histogram plot.")

def main():
    parser = argparse.ArgumentParser(description='Generate core gene FASTA files with enhanced features')
    parser.add_argument('--input-dir', required=True, help='Input directory containing genome FAA files')
    parser.add_argument('--core-genes', required=True, help='File containing core gene cluster IDs')
    parser.add_argument('--allele-names', required=True, help='TSV file mapping original IDs to allele IDs')
    parser.add_argument('--output-dir', required=True, help='Output directory for FASTA files')
    parser.add_argument('--prefix', default="", help='Prefix for output files')
    parser.add_argument('--method', default="", help='Clustering method (cdhit or swift)')
    
    args = parser.parse_args()
    
    # Parse input files
    original_to_cluster, original_to_allele, cluster_to_originals = parse_allele_mappings(args.allele_names)
    core_genes = parse_core_genes(args.core_genes)
    
    # Collect sequences from all genomes
    cluster_sequences, genome_stats = collect_genome_sequences(
        args.input_dir, original_to_cluster, original_to_allele, core_genes
    )
    
    # Generate FASTA files and statistics
    generate_core_fastas(
        core_genes, args.output_dir, cluster_sequences, genome_stats, 
        prefix=args.prefix, method=args.method
    )
    
    print("\nProcess completed successfully!")

if __name__ == "__main__":
    main()
