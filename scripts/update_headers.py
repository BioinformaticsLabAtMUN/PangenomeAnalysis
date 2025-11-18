#!/usr/bin/env python3
import os
import sys
import json
import re
from pathlib import Path

def convert_header_format(header, taxid_dict, debug=True):
    """
    Convert a FASTA header to include taxid for SwiftOrtho.
    Format: >taxid|protein_id
    """
    if debug:
        print(f"\nDEBUG: Converting header: {header}")
    
    # Remove leading '>' if present
    header = header.lstrip('>')
    
    # Find text within square brackets
    matches = re.findall(r"\[(.+?)\]", header)
    
    if debug:
        print(f"DEBUG: Found matches: {matches}")
    
    if not matches:
        if debug:
            print(f"WARNING: No organism found in brackets for header: {header}")
        return f">{header}"
        
    # Extract the organism from the last match
    organism = matches[-1].strip()
    
    if debug:
        print(f"DEBUG: Extracted organism: {organism}")
    
    taxid = taxid_dict.get(organism)
    
    if taxid:
        if debug:
            print(f"DEBUG: Found taxid: {taxid}")
        # Get the sequence ID (everything before the first space)
        seq_id = header.split()[0]
        new_header = f">{taxid}|{seq_id}"
        if debug:
            print(f"DEBUG: Created new header: {new_header}")
        return new_header
    else:
        if debug:
            print(f"WARNING: No taxid found for organism: {organism}")
        return f">{header}"

def verify_file_contents(filepath, description):
    """Verify file exists and has content."""
    path = Path(filepath)
    if not path.exists():
        raise FileNotFoundError(f"{description} not found: {filepath}")
    
    if path.stat().st_size == 0:
        raise ValueError(f"{description} is empty: {filepath}")
    
    with open(filepath, 'r') as f:
        first_line = f.readline().strip()
        if not first_line:
            raise ValueError(f"{description} exists but appears to be empty (no first line)")
        print(f"DEBUG: First line of {description}: {first_line}")

def update_fasta_headers(nr_faa_input, taxid_json, nr_faa_output):
    """
    Update FASTA headers to include taxid for SwiftOrtho.
    No shared headers processing.
    """
    # Verify input files exist and have content
    print("\n=== Verifying Input Files ===")
    verify_file_contents(nr_faa_input, "Non-redundant FAA input")
    verify_file_contents(taxid_json, "Taxonomy JSON")
    
    print(f"\n=== Loading Taxonomy ===")
    try:
        with open(taxid_json, 'r') as json_file:
            taxid_dict = json.load(json_file)
        print(f"Successfully loaded {len(taxid_dict)} taxonomy entries")
        print("Sample entries:")
        for k, v in list(taxid_dict.items())[:3]:
            print(f"  {k}: {v}")
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON file: {e}")
        raise
    except Exception as e:
        print(f"Error loading taxonomy file: {e}")
        raise

    print(f"\n=== Processing Non-redundant Sequences ===")
    headers_processed = 0
    sequences_processed = 0
    try:
        with open(nr_faa_input, 'r') as in_faa, open(nr_faa_output, 'w') as out_faa:
            for line in in_faa:
                if line.startswith('>'):
                    headers_processed += 1
                    new_header = convert_header_format(line.strip(), taxid_dict)
                    out_faa.write(f"{new_header}\n")
                else:
                    sequences_processed += 1
                    out_faa.write(line)
        print(f"Processed {headers_processed} headers and {sequences_processed} sequence lines")
    except Exception as e:
        print(f"Error processing non-redundant sequences: {e}")
        raise

    # Verify output file was created and has content
    print("\n=== Verifying Output File ===")
    verify_file_contents(nr_faa_output, "Non-redundant FAA output")

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python update_headers.py <nr_faa_input> <taxid_json> <nr_faa_output>")
        sys.exit(1)

    nr_faa_input = sys.argv[1]
    taxid_json = sys.argv[2]
    nr_faa_output = sys.argv[3]

    print("\n=== Starting Header Update Process ===")
    print(f"Input FAA: {nr_faa_input}")
    print(f"Taxonomy JSON: {taxid_json}")
    print(f"Output FAA: {nr_faa_output}")

    try:
        update_fasta_headers(nr_faa_input, taxid_json, nr_faa_output)
        print("\n=== Header Update Process Completed Successfully ===")
    except Exception as e:
        print(f"\nERROR: Header update process failed: {e}")
        sys.exit(1)
