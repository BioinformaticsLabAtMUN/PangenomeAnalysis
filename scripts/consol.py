#!/usr/bin/env python3
import os
import sys
import hashlib
import re

def extract_header_info(header_line):
    """
    Extracts protein ID, description, and organism (if present) from a header line.
    The header_line should not include the initial ">".
    Returns a tuple of (protein_id, description, organism).
    """
    if not header_line:
        return None, None, None
    parts = header_line.split(maxsplit=1)
    protein_id = parts[0] if parts else None
    description = parts[1] if len(parts) > 1 else ""
    # Extract organism from square brackets if present
    organism = None
    organism_match = re.search(r"\[(.+?)\]", header_line)
    if organism_match:
        organism = organism_match.group(1)
    return protein_id, description, organism

def consolidate_seqs(input_directory, output_nr_faa, output_swift_faa, output_shared_headers, 
                     output_missing_headers, output_descriptions, output_organism_map):
    """
    Merges all .faa files under input_directory into a single non-redundant set.
    
    Produces two FASTA outputs:
      - output_nr_faa: FASTA with simplified headers (only protein IDs)
      - output_swift_faa: FASTA with complete original headers
     
    Also produces:
      - output_shared_headers: if multiple protein IDs share an identical sequence, they are listed (tab-separated)
      - output_missing_headers: any protein IDs missing a sequence
      - output_descriptions: mapping of protein ID to its description
      - output_organism_map: mapping of protein ID to its organism
    """
    # Dictionary keyed by sequence hash. For each unique sequence, we store:
    #   "ids": list of protein IDs (for simplified FASTA)
    #   "full_headers": list of complete original headers (for SwiftOrtho)
    #   "sequence": the protein sequence
    #   "organisms": list of organism strings (if found)
    non_redundant = {}
    encounter_order = []
    missing_headers = []
    protein_descriptions = {}  # protein_id -> description
    protein_organisms = {}     # protein_id -> organism

    def process_header_seq(protein_id, full_header, seq, organism, description):
        if not protein_id:
            return
        if not seq:
            missing_headers.append(protein_id)
            return

        # Store description and organism if not already recorded
        if protein_id not in protein_descriptions:
            protein_descriptions[protein_id] = description
        if organism and protein_id not in protein_organisms:
            protein_organisms[protein_id] = organism

        # Compute hash of the sequence
        seqhash = hashlib.sha256(seq.encode('utf-8')).digest()
        if seqhash not in non_redundant:
            non_redundant[seqhash] = {
                "ids": [protein_id],
                "full_headers": [full_header],
                "sequence": seq,
                "organisms": [organism] if organism else []
            }
            encounter_order.append(seqhash)
        else:
            # Append new protein ID and full header if not already present
            if protein_id not in non_redundant[seqhash]["ids"]:
                non_redundant[seqhash]["ids"].append(protein_id)
            if full_header not in non_redundant[seqhash]["full_headers"]:
                non_redundant[seqhash]["full_headers"].append(full_header)
            if organism and organism not in non_redundant[seqhash]["organisms"]:
                non_redundant[seqhash]["organisms"].append(organism)

    # Traverse the input directory to process all .faa files
    for root, _, files in os.walk(input_directory):
        for fname in files:
            if fname.endswith('.faa'):
                faa_path = os.path.join(root, fname)
                print(f"Processing file: {faa_path}")
                with open(faa_path, 'r') as f:
                    current_full_header = None  # stores the complete header (with '>')
                    current_seq = []
                    current_description = ""
                    current_organism = None
                    for line in f:
                        line_strip = line.rstrip("\n")
                        if line_strip.startswith(">"):
                            # Process previous record if present
                            if current_full_header is not None:
                                pid, _, _ = extract_header_info(current_full_header[1:].strip())
                                if pid:
                                    process_header_seq(pid, current_full_header, "".join(current_seq), current_organism, current_description)
                                else:
                                    print(f"Warning: Malformed header encountered: '{current_full_header}'", file=sys.stderr)
                            # Start new record
                            current_full_header = line_strip  # keep full header unmodified
                            header_content = line_strip[1:].strip()  # header without '>'
                            pid, description, organism = extract_header_info(header_content)
                            current_description = description
                            current_organism = organism
                            current_seq = []
                        else:
                            current_seq.append(line_strip)
                    # Process last record in the file
                    if current_full_header is not None:
                        pid, _, _ = extract_header_info(current_full_header[1:].strip())
                        if pid:
                            process_header_seq(pid, current_full_header, "".join(current_seq), current_organism, current_description)
                        else:
                            print(f"Warning: Malformed header encountered at end of file: '{current_full_header}'", file=sys.stderr)

    # Write output files
    print("Writing output files...")
    with open(output_nr_faa, 'w') as nr_out, \
         open(output_swift_faa, 'w') as swift_out, \
         open(output_shared_headers, 'w') as shared_out, \
         open(output_missing_headers, 'w') as missing_out, \
         open(output_descriptions, 'w') as desc_out, \
         open(output_organism_map, 'w') as org_out:

        # Write protein descriptions mapping
        for protein_id, description in protein_descriptions.items():
            desc_out.write(f"{protein_id}\t{description}\n")
            
        # Write organism mapping
        for protein_id, organism in protein_organisms.items():
            org_out.write(f"{protein_id}\t{organism}\n")

        # For each unique sequence, write to the two FASTA outputs
        for seqhash in encounter_order:
            info = non_redundant[seqhash]
            ids = info["ids"]                # list of protein IDs (simplified headers)
            full_headers = info["full_headers"]  # list of complete original headers
            seq = info["sequence"]

            # Use the first encountered values as the representative headers
            rep_id = ids[0]
            rep_full_header = full_headers[0]

            # Write to simplified FASTA (only protein ID in header)
            nr_out.write(f">{rep_id}\n")
            nr_out.write(seq + "\n")

            # Write to SwiftOrtho FASTA (complete original header)
            swift_out.write(f"{rep_full_header}\n")
            swift_out.write(seq + "\n")

            # Write shared headers only for the simplified FASTA version if there are duplicates
            if len(ids) > 1:
                unique_ids = list(dict.fromkeys(ids))
                shared_out.write("\t".join(unique_ids) + "\n")

        # Write any missing headers
        for mh in missing_headers:
            missing_out.write(f"{mh}\n")

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print("Usage: python consolidate_sequences.py <input_directory> <output_nr_faa> <output_swift_faa> "
              "<output_shared_headers> <output_missing_headers> <output_descriptions> <output_organism_map>")
        sys.exit(1)

    input_dir = sys.argv[1]
    output_nr_faa = sys.argv[2]
    output_swift_faa = sys.argv[3]
    output_shared_headers = sys.argv[4]
    output_missing_headers = sys.argv[5]
    output_descriptions = sys.argv[6]
    output_organism_map = sys.argv[7]

    if not os.path.isdir(input_dir):
        print(f"Error: Input directory '{input_dir}' does not exist.")
        sys.exit(1)

    consolidate_seqs(input_dir, output_nr_faa, output_swift_faa, output_shared_headers, 
                     output_missing_headers, output_descriptions, output_organism_map)

