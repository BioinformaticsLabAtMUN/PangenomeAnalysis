#!/usr/bin/env python3
import os
import sys
import argparse
import scipy.sparse
from pathlib import Path
import sparse_utils

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

def load_header_to_allele(allele_names_file, shared_headers_file=None):
    """Load mapping from original headers to allele names (handles synonyms)."""
    # load shared synonyms if provided
    shared = {}
    if shared_headers_file and os.path.exists(shared_headers_file):
        with open(shared_headers_file) as f:
            for line in f:
                headers = line.rstrip('\n').split('\t')
                if headers:
                    rep = headers[0]
                    for h in headers:
                        shared[h] = rep

    # detect old vs new format
    with open(allele_names_file) as f:
        first = f.readline().rstrip('\n').split('\t')
    has_header = len(first) >= 2 and first[0] == "Original_ID" and first[1] == "New_ID"

    mapping = {}
    with open(allele_names_file) as f:
        if has_header:
            next(f)
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            if has_header:
                orig, allele = parts[0], parts[1]
            else:
                allele, orig = parts[0], parts[1]
            mapping[orig] = allele
            # also map any shared synonyms
            if orig in shared:
                mapping[shared[orig]] = allele

    print(f"Loaded {len(mapping)} header→allele mappings (with {len(shared)} shared synonyms)")
    return mapping

def build_genetic_feature_tables(strain_faa_paths, header_to_allele):
    """Build presence/absence matrices for alleles and genes."""
    # extract strain names
    strains = [Path(p).stem for p in strain_faa_paths]

    # all alleles
    allele_list = sorted(set(header_to_allele.values()))
    allele_idx = {a:i for i,a in enumerate(allele_list)}

    # derive genes by dropping trailing "A<number>"
    gene_list = []
    last = None
    for a in allele_list:
        g = a.rsplit('A',1)[0]
        if g != last:
            gene_list.append(g)
            last = g
    gene_idx = {g:i for i,g in enumerate(gene_list)}

    # prepare DOK matrices
    A = scipy.sparse.dok_matrix((len(allele_list), len(strains)), dtype=int)
    G = scipy.sparse.dok_matrix((len(gene_list),   len(strains)), dtype=int)

    for col, fasta in enumerate(strain_faa_paths):
        name = strains[col]
        print(f"→ processing {name} ({col+1}/{len(strains)})")
        with open(fasta) as f:
            for line in f:
                if not line.startswith('>'): 
                    continue
                raw_id = line[1:].split(None,1)[0]
                hdr    = parse_fasta_header(raw_id)
                if hdr not in header_to_allele:
                    continue
                allele = header_to_allele[hdr]
                ai = allele_idx[allele]
                gi = gene_idx[allele.rsplit('A',1)[0]]
                A[ai, col] = 1
                G[gi, col] = 1

    print("Allele nonzeros:", A.count_nonzero())
    print("Gene   nonzeros:", G.count_nonzero())

    dfA = sparse_utils.LightSparseDataFrame(
        index=allele_list, columns=strains, data=A.tocoo())
    dfG = sparse_utils.LightSparseDataFrame(
        index=gene_list,   columns=strains, data=G.tocoo())
    return dfA, dfG

def main():
    p = argparse.ArgumentParser(description="Build pangenome tables")
    p.add_argument("--input-dir",    required=True)
    p.add_argument("--output-dir",   required=True)
    p.add_argument("--name",         required=True)
    p.add_argument("--allele-names", required=True)
    p.add_argument("--shared-headers", required=False)
    args = p.parse_args()

    indir  = Path(args.input_dir)
    outdir = Path(args.output_dir)
    outdir.mkdir(parents=True, exist_ok=True)

    # find all FASTA files
    paths = sorted(indir.rglob("*.fasta")) + sorted(indir.rglob("*.faa")) + sorted(indir.rglob("*.fa"))
    # drop anything with these substrings
    paths = [p for p in paths if not any(x in p.name for x in ("renamed","cdhit","nr","swift"))]
    if not paths:
        sys.exit(" No FASTA files found in " + str(indir))

    print(f"Found {len(paths)} FASTAs → building tables")
    mapping = load_header_to_allele(args.allele_names, args.shared_headers)
    dfA, dfG  = build_genetic_feature_tables(paths, mapping)

    # save
    prefix = outdir/args.name
    print("Saving allele matrix…")
    dfA.to_npz(str(prefix) + "_strain_by_allele.npz")
    print("Saving gene   matrix…")
    dfG.to_npz(str(prefix) + "_strain_by_gene.npz")
    print("Done.")

if __name__ == "__main__":
    main()
