#!/usr/bin/env python3
"""
validate_pangenome.py

This script validates two aspects of your pan-genome pipeline:
  1. Gene Table Validation: For each genome, it converts the raw FASTA headers 
     (from the original files in --input-dir) into renamed headers using the allele
     names mapping (from --allele-names), then collapses them to gene names (by 
     removing the allele suffix) and compares that set to the gene x genome table.
     
  2. Allele Table Validation: For each genome, it converts raw FASTA headers into 
     renamed headers and compares that set to the allele x genome table.

Both gene and allele tables are stored in LSDF (NPZ plus labels) format and are loaded 
via your existing sparse_utils module.
"""

import os
import sys
import argparse
import hashlib
import traceback
from datetime import datetime
import numpy as np
import scipy.sparse
from multiprocessing import Pool, cpu_count
import time

# Import your existing sparse_utils module (as in your original code)
import sparse_utils

# --- Helper Functions ---

def log_message(message, f=None):
    """Prints a timestamped message to stdout and, if provided, to file f."""
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    formatted = f"[{timestamp}] {message}"
    print(formatted)
    if f:
        f.write(formatted + "\n")
        f.flush()

def __get_gene_from_allele__(allele):
    """
    Converts a renamed allele (e.g. 'Strep_C13798A0') to its gene (e.g. 'Strep_C13798')
    by removing the trailing allele part. Assumes the allele name splits on 'A'.
    Optimized version using rsplit for efficiency.
    """
    if 'A' in allele:
        return allele.rsplit('A', 1)[0]
    return allele
    
def batch_get_genes_from_alleles(alleles):
    """
    Efficiently converts a set of alleles to their gene names.
    """
    return {__get_gene_from_allele__(allele) for allele in alleles}

def __hash_sequence__(seq):
    """Returns the SHA-256 digest (bytes) for a given sequence string."""
    return hashlib.sha256(seq.encode('utf-8')).digest()

def __get_genome_from_filename__(filepath):
    """Extracts the genome name from a filepath by stripping directory and extension."""
    filename = os.path.basename(filepath)  # More efficient than os.path.split
    return os.path.splitext(filename)[0]

def parse_fasta_header(raw):
    """
    Given the first token of a FASTA header,
    extract and return the protein ID.
    
    Handles:
    - taxid|protein_id format (e.g., "73044|A0A4P6TQ31")
    - UniProt format (e.g., "tr|D7BUT1|D7BUT1_STRBB")
    - Plain IDs without pipes
    """
    parts = raw.split('|')
    if len(parts) == 3:  # UniProt format
        return parts[1]
    elif len(parts) == 2:  # taxid|protein_id format
        return parts[1]
    return raw  # Return raw for any other format

def load_fasta_headers(fasta_file):
    """
    Loads headers (first token after '>') from a FASTA file.
    Returns a set of raw headers. Optimized to split only once.
    """
    headers = set()
    with open(fasta_file, 'r') as f:
        for line in f:
            if line.startswith('>'):
                # Split only once to get the first token (more efficient)
                header = line[1:].strip().split(None, 1)[0]
                headers.add(header)
    return headers

def load_inverse_allele_mapping(allele_names_file):
    """
    Loads the allele names mapping from a TSV file.
    Supports both old format and new format with header row.
    Returns a dict mapping each original header to its renamed allele.
    Optimized for larger files.
    """
    # Check if file has a header row by examining first line
    with open(allele_names_file, 'r') as f:
        first_line = f.readline().strip().split('\t')
        has_header = len(first_line) > 1 and first_line[0] == "Original_ID" and first_line[1] == "New_ID"
    
    mapping = {}
    # Reopen file since we read the first line
    with open(allele_names_file, 'r') as f:
        # Skip header if present
        if has_header:
            next(f)
            
        for line in f:
            tokens = line.strip().split('\t')
            if len(tokens) < 2:
                continue
                
            if has_header:
                # New format: Original_ID, New_ID, etc.
                # For taxid|protein_id format, we need to map both the full ID and just the protein_id
                orig_id = tokens[0]
                new_id = tokens[1]
                mapping[orig_id] = new_id
                
                # If the original_id contains a pipe, also map the protein_id alone
                if '|' in orig_id:
                    protein_id = parse_fasta_header(orig_id)
                    mapping[protein_id] = new_id
            else:
                # Old format: Renamed, Original_header
                renamed = tokens[0]
                # Add each original header to mapping
                for orig in tokens[1:]:
                    mapping[orig] = renamed
                    # If the original contains a pipe, also map the protein_id alone
                    if '|' in orig:
                        protein_id = parse_fasta_header(orig)
                        mapping[protein_id] = renamed
                    
    return mapping

def get_lsdf_column_as_set(df, genome_name):
    """
    Efficiently extracts a column from a LSDF DataFrame as a set of index values.
    """
    if genome_name not in df.columns:
        return set()
        
    col_idx = df.columns.get_loc(genome_name) if hasattr(df.columns, 'get_loc') else list(df.columns).index(genome_name)
    # Get column using CSC format (column-oriented, efficient for column access)
    col = df.data.tocsc()[:, col_idx]
    # Get non-zero indices
    nonzero_indices = col.nonzero()[0]
    # Map to actual index values
    return {df.index[i] for i in nonzero_indices}

# --- Validation Functions (Parallelized) ---

def validate_genome_combined(args):
    """
    Worker function for parallel validation of both gene and allele tables in one pass.
    Returns ((genome_name, gene_result, gene_diff), (genome_name, allele_result, allele_diff))
    Where result is: True for match, False for mismatch, None for error
    """
    fasta_path, df_genes_obj, df_alleles_obj, allele_mapping = args
    genome_name = "unknown"  # Default in case extraction fails
    
    try:
        # Step 1: Extract genome name from filename
        try:
            genome_name = __get_genome_from_filename__(fasta_path)
            # Debug line to show processing started for this genome
            # print(f"DEBUG: Started processing {genome_name} from {fasta_path}")
        except Exception as e:
            return ((genome_name, None, f"Error extracting genome name: {str(e)}"), 
                    (genome_name, None, f"Error extracting genome name: {str(e)}"))
        
        # Step 2: Load FASTA headers
        try:
            raw_headers = load_fasta_headers(fasta_path)
            # Debug: print first few headers
            # sample_headers = list(raw_headers)[:3] if raw_headers else []
            # print(f"DEBUG: {genome_name} - {len(raw_headers)} raw headers. Sample: {sample_headers}")
        except Exception as e:
            return ((genome_name, None, f"Error loading FASTA headers: {str(e)}"), 
                    (genome_name, None, f"Error loading FASTA headers: {str(e)}"))
            
        # Step 3: Convert raw headers to renamed headers
        try:
            renamed_headers = set()
            for h in raw_headers:
                # First try the full header
                if h in allele_mapping:
                    renamed_headers.add(allele_mapping[h])
                else:
                    # Try to extract and use the protein_id part
                    protein_id = parse_fasta_header(h)
                    if protein_id in allele_mapping:
                        renamed_headers.add(allele_mapping[protein_id])
                    else:
                        renamed_headers.add(h)  # Keep original if not in mapping
            # Debug: print count and any issues
            # print(f"DEBUG: {genome_name} - {len(renamed_headers)} renamed headers")
        except Exception as e:
            return ((genome_name, None, f"Error converting headers using mapping: {str(e)}"), 
                    (genome_name, None, f"Error converting headers using mapping: {str(e)}"))
        
        # Step 4: Generate gene headers by collapsing allele identifiers
        try:
            gene_headers = batch_get_genes_from_alleles(renamed_headers)
            # Debug: print counts
            # print(f"DEBUG: {genome_name} - {len(gene_headers)} gene headers")
        except Exception as e:
            return ((genome_name, None, f"Error generating gene headers: {str(e)}"), 
                    (genome_name, None, f"Error generating gene headers: {str(e)}"))
        
        # Step 5: Get gene table data for this genome
        try:
            genes_lsdf = get_lsdf_column_as_set(df_genes_obj, genome_name)
            # Debug: print count
            # print(f"DEBUG: {genome_name} - {len(genes_lsdf)} genes in LSDF")
        except Exception as e:
            return ((genome_name, None, f"Error retrieving gene data from matrix: {str(e)}"), 
                    (genome_name, None, f"Error retrieving gene data from matrix: {str(e)}"))
        
        # Step 6: Get allele table data for this genome
        try:
            alleles_lsdf = get_lsdf_column_as_set(df_alleles_obj, genome_name)
            # Debug: print count
            # print(f"DEBUG: {genome_name} - {len(alleles_lsdf)} alleles in LSDF")
        except Exception as e:
            return ((genome_name, None, f"Error retrieving allele data from matrix: {str(e)}"), 
                    (genome_name, None, f"Error retrieving allele data from matrix: {str(e)}"))
        
        # Step 7: Check genes
        try:
            if not genes_lsdf:
                gene_result = (genome_name, None, "Not found in gene table")
            elif genes_lsdf != gene_headers:
                diff = genes_lsdf.symmetric_difference(gene_headers)
                # Limit diff size for reporting to avoid huge output
                if len(diff) > 10:
                    gene_result = (genome_name, False, f"{len(diff)} differences (showing first 10): {list(diff)[:10]}")
                else:
                    gene_result = (genome_name, False, diff)
            else:
                gene_result = (genome_name, True, None)
        except Exception as e:
            gene_result = (genome_name, None, f"Error comparing gene sets: {str(e)}")
            
        # Step 8: Check alleles
        try:
            if not alleles_lsdf:
                allele_result = (genome_name, None, "Not found in allele table")
            elif alleles_lsdf != renamed_headers:
                diff = alleles_lsdf.symmetric_difference(renamed_headers)
                # Limit diff size for reporting to avoid huge output
                if len(diff) > 10:
                    allele_result = (genome_name, False, f"{len(diff)} differences (showing first 10): {list(diff)[:10]}")
                else:
                    allele_result = (genome_name, False, diff)
            else:
                allele_result = (genome_name, True, None)
        except Exception as e:
            allele_result = (genome_name, None, f"Error comparing allele sets: {str(e)}")
            
        # Success - return results
        return (gene_result, allele_result)
    
    except Exception as e:
        # Catch-all for any unexpected errors
        import traceback
        tb = traceback.format_exc()
        error_msg = f"Unexpected error processing {genome_name}: {str(e)}\n{tb}"
        return ((genome_name, None, error_msg), (genome_name, None, error_msg))

def validate_tables_from_files(df_genes, df_alleles, input_dir, allele_mapping, workers=None, log_group=10, output_file=None, batch_size=100):
    """
    Parallelized version of combined gene and allele table validation.
    Uses multiple workers to process genomes simultaneously, checking both tables in one pass.
    Processes genomes in batches to reduce memory pressure.
    Returns (gene_inconsistencies, allele_inconsistencies)
    """
    start_time = time.time()
    log_message("Validating gene and allele tables from raw FASTA files (combined pass)...", output_file)
    
    # Gather genome FASTA files
    fasta_paths = []
    for root, _, files in os.walk(input_dir):
        for file in files:
            if (file.endswith('.faa') or file.endswith('.fa') or file.endswith('.fasta')) and not any(x in file for x in ['renamed', 'cdhit', 'nr']):
                fasta_paths.append(os.path.join(root, file))
    
    sorted_paths = sorted(fasta_paths)
    log_message(f"Found {len(sorted_paths)} genome FASTA files to process", output_file)
    
    # Sample test: Try to process one genome synchronously to detect any issues
    try:
        log_message("Testing single genome processing before launching parallel jobs...", output_file)
        test_path = sorted_paths[0]
        test_genome = __get_genome_from_filename__(test_path)
        log_message(f"Testing with genome: {test_genome} from {test_path}", output_file)
        
        # Load headers
        raw_headers = load_fasta_headers(test_path)
        log_message(f"Test genome has {len(raw_headers)} raw headers", output_file)
        
        # Try mapping headers
        sample_headers = list(raw_headers)[:3]
        log_message(f"Sample raw headers: {sample_headers}", output_file)
        
        # Test mapping to renamed headers
        renamed_headers = set()
        for h in raw_headers:
            # First try the full header
            if h in allele_mapping:
                renamed_headers.add(allele_mapping[h])
            else:
                # Try to extract and use the protein_id part
                protein_id = parse_fasta_header(h)
                if protein_id in allele_mapping:
                    renamed_headers.add(allele_mapping[protein_id])
                else:
                    renamed_headers.add(h)  # Keep original if not in mapping
        
        log_message(f"Test genome has {len(renamed_headers)} renamed headers", output_file)
        
        # Test gene extraction
        gene_headers = batch_get_genes_from_alleles(renamed_headers)
        log_message(f"Test genome has {len(gene_headers)} gene headers", output_file)
        
        # Test matrix access
        genes_lsdf = get_lsdf_column_as_set(df_genes, test_genome)
        log_message(f"Test genome has {len(genes_lsdf)} genes in LSDF", output_file)
        
        alleles_lsdf = get_lsdf_column_as_set(df_alleles, test_genome)
        log_message(f"Test genome has {len(alleles_lsdf)} alleles in LSDF", output_file)
        
        log_message("Test processing successful, proceeding with parallel processing", output_file)
    except Exception as e:
        import traceback
        tb = traceback.format_exc()
        log_message(f"ERROR in test processing: {str(e)}\n{tb}", output_file)
        log_message("Will proceed with sequential processing due to test failure", output_file)
        workers = 1  # Fall back to single process mode
    
    # Use multiprocessing if workers > 1
    num_workers = min(workers or max(1, min(cpu_count(), 16)), 16)  # Use at most 16 cores by default
    log_message(f"Using {num_workers} worker processes", output_file)
    
    # Process in batches to reduce memory pressure
    gene_inconsistencies = 0
    allele_inconsistencies = 0
    processed = 0
    total_genomes = len(sorted_paths)
    
    # Process genomes in batches to reduce memory usage
    for batch_start in range(0, total_genomes, batch_size):
        batch_end = min(batch_start + batch_size, total_genomes)
        current_batch = sorted_paths[batch_start:batch_end]
        log_message(f"Processing batch {batch_start//batch_size + 1}/{(total_genomes + batch_size - 1)//batch_size}: "
                  f"genomes {batch_start+1}-{batch_end}", output_file)
        
        # Prepare task arguments for this batch
        task_args = [(path, df_genes, df_alleles, allele_mapping) for path in current_batch]
        
        # Process this batch
        if num_workers > 1:
            # Parallel processing for this batch
            with Pool(processes=num_workers) as pool:
                try:
                    # Use imap with a small chunksize to get more frequent output
                    results_iter = pool.imap(validate_genome_combined, task_args, chunksize=1)
                    
                    # Process results as they arrive
                    for result in results_iter:
                        try:
                            # Unpack result
                            (genome_name, gene_valid, gene_msg), (_, allele_valid, allele_msg) = result
                            processed += 1
                            
                            # Process and log results
                            _process_validation_result(genome_name, gene_valid, gene_msg, allele_valid, allele_msg, 
                                                     output_file, processed, len(sorted_paths), start_time, log_group,
                                                     gene_inconsistencies, allele_inconsistencies)
                            
                            # Update inconsistency counts
                            if gene_valid is False and isinstance(gene_msg, set):
                                gene_inconsistencies += len(gene_msg)
                            elif gene_valid is False:
                                gene_inconsistencies += 1
                                
                            if allele_valid is False and isinstance(allele_msg, set):
                                allele_inconsistencies += len(allele_msg)
                            elif allele_valid is False:
                                allele_inconsistencies += 1
                                
                        except Exception as e:
                            log_message(f"ERROR processing result: {str(e)}", output_file)
                            
                except Exception as e:
                    import traceback
                    tb = traceback.format_exc()
                    log_message(f"ERROR in parallel processing batch: {str(e)}\n{tb}", output_file)
                    log_message("Switching to sequential processing for remaining genomes", output_file)
                    num_workers = 1  # Fall back to sequential mode for future batches
        
        else:
            # Sequential processing for this batch
            for args in task_args:
                try:
                    result = validate_genome_combined(args)
                    (genome_name, gene_valid, gene_msg), (_, allele_valid, allele_msg) = result
                    processed += 1
                    
                    # Process and log results
                    _process_validation_result(genome_name, gene_valid, gene_msg, allele_valid, allele_msg, 
                                             output_file, processed, len(sorted_paths), start_time, log_group,
                                             gene_inconsistencies, allele_inconsistencies)
                    
                    # Update inconsistency counts
                    if gene_valid is False and isinstance(gene_msg, set):
                        gene_inconsistencies += len(gene_msg)
                    elif gene_valid is False:
                        gene_inconsistencies += 1
                        
                    if allele_valid is False and isinstance(allele_msg, set):
                        allele_inconsistencies += len(allele_msg)
                    elif allele_valid is False:
                        allele_inconsistencies += 1
                        
                except Exception as e:
                    import traceback
                    tb = traceback.format_exc()
                    genome_name = __get_genome_from_filename__(args[0]) if args else "unknown"
                    log_message(f"ERROR processing genome {genome_name}: {str(e)}\n{tb}", output_file)
    
    elapsed = time.time() - start_time
    log_message(f"Combined validation complete. Gene inconsistencies: {gene_inconsistencies}, Allele inconsistencies: {allele_inconsistencies}", output_file)
    log_message(f"Total time: {elapsed:.1f}s, Average: {elapsed/total_genomes:.2f}s per genome", output_file)
    
    return (gene_inconsistencies, allele_inconsistencies)

def _process_validation_result(genome_name, gene_valid, gene_msg, allele_valid, allele_msg, 
                              output_file, processed, total, start_time, log_group,
                              gene_inconsistencies, allele_inconsistencies):
    """Helper function to process and log validation results"""
    # Process gene result
    if gene_valid is None:
        log_message(f"WARNING: Genome {genome_name}: {gene_msg} (Gene Table)", output_file)
    elif not gene_valid:
        if isinstance(gene_msg, set):
            log_message(f"{genome_name} gene table: {len(gene_msg)} inconsistencies", output_file)
        else:
            log_message(f"{genome_name} gene table: 1 inconsistency", output_file)
    else:
        log_message(f"{genome_name} gene table: 0 inconsistencies", output_file)
        
    # Process allele result
    if allele_valid is None:
        log_message(f"WARNING: Genome {genome_name}: {allele_msg} (Allele Table)", output_file)
    elif not allele_valid:
        if isinstance(allele_msg, set):
            log_message(f"{genome_name} allele table: {len(allele_msg)} inconsistencies", output_file)
        else:
            log_message(f"{genome_name} allele table: 1 inconsistency", output_file)
    else:
        log_message(f"{genome_name} allele table: 0 inconsistencies", output_file)
        
    # Log progress periodically
    if processed % log_group == 0:
        elapsed = time.time() - start_time
        log_message(f"Processed {processed}/{total} genomes ({processed/total*100:.1f}%). "
                  f"Time elapsed: {elapsed:.1f}s, Avg: {elapsed/processed:.2f}s per genome", output_file)

# --- Main Function ---

def main():
    try:
        parser = argparse.ArgumentParser(description="Validate pan-genome tables (LSDF format) using allele mapping")
        parser.add_argument('--gene-matrix', required=True,
                          help="Path to gene x genome NPZ file (LSDF format)")
        parser.add_argument('--allele-matrix', required=True,
                          help="Path to allele x genome NPZ file (LSDF format)")
        parser.add_argument('--input-dir', required=True,
                          help="Directory containing original genome FASTA files (e.g., .faa)")
        parser.add_argument('--allele-names', required=True,
                          help="Path to allele names mapping TSV file (e.g., Strep_allele_names.tsv)")
        parser.add_argument('--output-dir', required=False, default=".",
                          help="Output directory (unused, for pipeline compatibility)")
        parser.add_argument('--workers', type=int, default=None,
                          help="Number of worker processes (default: auto-detect, up to 16)")
        parser.add_argument('--log-group', type=int, default=10,
                          help="Log progress after processing this many genomes (default: 10)")
        parser.add_argument('--batch-size', type=int, default=100,
                          help="Process genomes in batches of this size to reduce memory usage (default: 100)")
        args = parser.parse_args()

        # Create explicit output file handler
        log_file_path = os.path.join(args.output_dir, 'validation_summary.txt')
        with open(log_file_path, 'w') as log_file:
            overall_start = time.time()
            log_message("Starting pan-genome validation", log_file)
            log_message(f"Gene matrix: {args.gene_matrix}", log_file)
            log_message(f"Allele matrix: {args.allele_matrix}", log_file)
            log_message(f"Allele names: {args.allele_names}", log_file)
            log_message(f"Workers: {args.workers or 'auto'}", log_file)
            log_message(f"Batch size: {args.batch_size}", log_file)
            
            # Validate file existence
            if not os.path.exists(args.gene_matrix):
                raise FileNotFoundError(f"Gene matrix not found: {args.gene_matrix}")
            if not os.path.exists(args.allele_matrix):
                raise FileNotFoundError(f"Allele matrix not found: {args.allele_matrix}")
            if not os.path.exists(args.allele_names):
                raise FileNotFoundError(f"Allele names file not found: {args.allele_names}")

            # Load data
            loading_start = time.time()
            log_message("Loading sparse matrices and allele mapping...", log_file)
            
            df_genes = sparse_utils.read_lsdf(args.gene_matrix)
            df_alleles = sparse_utils.read_lsdf(args.allele_matrix)
            allele_mapping = load_inverse_allele_mapping(args.allele_names)
            
            loading_time = time.time() - loading_start
            log_message(f"Data loaded in {loading_time:.2f}s", log_file)
            log_message(f"Gene matrix shape: {df_genes.shape}", log_file)
            log_message(f"Allele matrix shape: {df_alleles.shape}", log_file)
            log_message(f"Allele mapping entries: {len(allele_mapping)}", log_file)

            log_message("\n--- Combined Gene and Allele Table Validation ---", log_file)
            gene_incon, allele_incon = validate_tables_from_files(
                df_genes, df_alleles, args.input_dir, allele_mapping, 
                workers=args.workers, log_group=args.log_group, 
                output_file=log_file, batch_size=args.batch_size
            )
            
            overall_time = time.time() - overall_start
            log_message(f"\nValidation complete in {overall_time:.2f}s", log_file)
            log_message(f"Total gene table inconsistencies: {gene_incon}", log_file)
            log_message(f"Total allele table inconsistencies: {allele_incon}", log_file)
            log_message(f"Results written to {log_file_path}", log_file)
            
            # Final summary - writing this directly, not just to log file
            print(f"\nValidation complete. Results:")
            print(f"- Gene table inconsistencies: {gene_incon}")
            print(f"- Allele table inconsistencies: {allele_incon}")
            print(f"- Detailed results saved to: {log_file_path}")

    except Exception as e:
        error_msg = f"Validation failed: {str(e)}\nTraceback:\n{traceback.format_exc()}"
        print(error_msg)
        try:
            if 'log_file' in locals() and log_file:
                log_file.write(error_msg)
        except:
            pass
        sys.exit(1)

if __name__ == '__main__':
    main()
