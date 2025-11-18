#!/usr/bin/env python3
import argparse
import csv
import re
import sys
from pathlib import Path
import logging

# Set up logging
logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')
logger = logging.getLogger(__name__)

def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description='Check essential genes in core genome.')
    parser.add_argument('--core-genes', required=True, help='Path to core genes list file')
    parser.add_argument('--allele-names', required=True, help='Path to allele names mapping file')
    parser.add_argument('--sco-mapping', required=True, help='Path to SCO-NP mapping file')
    parser.add_argument('--output', required=True, help='Path to output TSV file')
    parser.add_argument('--method', required=True, help='Clustering method (cdhit or swift)')
    parser.add_argument('--verbose', '-v', action='store_true', help='Enable verbose output')
    return parser.parse_args()

def load_core_genes(core_genes_file):
    """Load core genes from file into a set."""
    core_genes = set()
    with open(core_genes_file, 'r') as f:
        for line in f:
            gene_id = line.strip()
            if gene_id:  # Skip empty lines
                core_genes.add(gene_id)
    return core_genes

def build_mappings_enhanced(allele_names_file):
    """
    Build various mappings from the enhanced TSV format allele names file:
    1. protein_to_cluster: Maps protein IDs to cluster IDs
    2. np_to_cluster: Maps only NP_* IDs to cluster IDs (with and without version numbers)
    """
    protein_to_cluster = {}
    np_to_cluster = {}
    cluster_to_proteins = {}
    
    with open(allele_names_file, 'r') as f:
        csvreader = csv.reader(f, delimiter='\t')
        header = next(csvreader)  # Skip header
        
        for row in csvreader:
            if len(row) < 4:  # Need at least Original_ID, New_ID, Cluster_ID, Allele_Number
                continue
                
            protein_id = row[0].strip()  # Original_ID
            cluster_id = row[2].strip()  # Cluster_ID
            
            # Map all original protein IDs to this cluster
            protein_to_cluster[protein_id] = cluster_id
            
            # For NP IDs, also store without version
            if protein_id.startswith('NP_'):
                base_id = protein_id.split('.')[0]
                np_to_cluster[base_id] = cluster_id
                np_to_cluster[protein_id] = cluster_id
                
            # Store mapping from cluster to proteins for reverse lookup
            if cluster_id not in cluster_to_proteins:
                cluster_to_proteins[cluster_id] = []
            cluster_to_proteins[cluster_id].append(protein_id)
    
    return protein_to_cluster, np_to_cluster, cluster_to_proteins

def check_essential_genes(sco_mapping_file, protein_to_cluster, np_to_cluster, cluster_to_proteins, core_genes, output_file, method):
    """
    Check which SCO genes are in the core genome.
    """
    not_found_count = 0
    core_count = 0
    not_core_count = 0
    
    with open(sco_mapping_file, 'r') as infile, open(output_file, 'w') as outfile:
        # Write header
        outfile.write("SCO_ID\tNP_ID\tCluster_ID\tGene_Function\tStatus\n")
        
        # Parse the CSV file
        csv_reader = csv.reader(infile)
        header = next(csv_reader)  # Skip header
        
        for row in csv_reader:
            if len(row) < 3:
                continue  # Skip rows with insufficient columns
            
            sco_id = row[0].strip()
            np_id = row[1].strip()
            gene_function = row[2].strip() if len(row) > 2 else ""
            
            # Look up cluster ID directly from NP ID mapping
            cluster_id = np_to_cluster.get(np_id, "")
            
            # If not found, try adding version numbers
            if not cluster_id:
                for version in ['.1', '.2', '.3']:
                    versioned_np_id = f"{np_id}{version}"
                    if versioned_np_id in protein_to_cluster:
                        cluster_id = protein_to_cluster[versioned_np_id]
                        logger.info(f"Found match using versioned ID: {versioned_np_id} -> {cluster_id}")
                        break
            
            # Determine status
            if not cluster_id:
                status = "NOT_FOUND"
                not_found_count += 1
                logger.debug(f"NOT_FOUND: {sco_id} - {np_id}")
            elif cluster_id in core_genes:
                status = "CORE"
                core_count += 1
            else:
                status = "NOT_CORE"
                not_core_count += 1
            
            # Write output row
            outfile.write(f"{sco_id}\t{np_id}\t{cluster_id}\t{gene_function}\t{status}\n")
    
    # Print summary statistics
    logger.info(f"Method: {method} - Summary: {core_count} CORE, {not_core_count} NOT_CORE, {not_found_count} NOT_FOUND")
    
    return core_count, not_core_count, not_found_count

def main():
    """Main function."""
    args = parse_args()
    
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    
    # Check if input files exist
    for input_file in [args.core_genes, args.allele_names, args.sco_mapping]:
        if not Path(input_file).is_file():
            logger.error(f"Input file {input_file} does not exist.")
            sys.exit(1)
    
    # Load core genes
    logger.info(f"Loading core genes from {args.core_genes}...")
    core_genes = load_core_genes(args.core_genes)
    logger.info(f"Loaded {len(core_genes)} core genes.")
    
    # Build protein-to-cluster mapping
    logger.info(f"Building protein-to-cluster mapping from {args.allele_names}...")
    protein_to_cluster, np_to_cluster, cluster_to_proteins = build_mappings_enhanced(args.allele_names)
    logger.info(f"Mapped {len(protein_to_cluster)} total protein IDs to cluster IDs.")
    logger.info(f"Mapped {len(np_to_cluster)} NP protein IDs to cluster IDs.")
    
    # Count how many NP_ proteins we have
    np_ids = [pid for pid in protein_to_cluster if pid.startswith('NP_')]
    logger.info(f"Found {len(np_ids)} NP IDs in the allele names file")
    if np_ids:
        logger.info(f"Sample NP IDs: {np_ids[:5]}")
    
    # Check essential genes
    logger.info(f"Checking essential genes from {args.sco_mapping}...")
    core_count, not_core_count, not_found_count = check_essential_genes(
        args.sco_mapping, protein_to_cluster, np_to_cluster, cluster_to_proteins, 
        core_genes, args.output, args.method
    )
    
    # If all genes are NOT_FOUND, suggest potential solutions
    if not_found_count > 0 and core_count == 0 and not_core_count == 0:
        logger.warning("All genes were marked as NOT_FOUND. This suggests a serious mapping issue.")
        logger.warning("Possible solutions:")
        logger.warning("1. Check if the NP IDs in SCO_NP_mapping.csv correspond to the same organism as your strain files")
        logger.warning("2. You may need to perform a BLAST search to map the SCO genes to your strain's genome")
        logger.warning("3. Check if there's a different version of the SCO mapping file that matches your strains")
    
    logger.info(f"Results written to {args.output}.")

if __name__ == "__main__":
    main()
