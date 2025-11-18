#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Extract dominant alleles per gene using sparse_utils LightSparseDataFrame.

Inputs (same interface you’re already using):
  --allele-matrix           .../Strep_cdhit_strain_by_allele.npz
  --allele-labels           .../Strep_cdhit_strain_by_allele.npz.labels.txt
  --allele-names            .../Strep_cdhit_allele_names.tsv
  --allele-faa              .../Strep_cdhit_renamed.fasta
  --core-genes              .../Strep_core_genes.txt
Outputs:
  --all-dominant-faa        .../Strep_cdhit_all_dominant_alleles.faa
  --all-dominant-summary    .../Strep_cdhit_all_dominant_summary.tsv
  --core-dominant-faa       .../Strep_cdhit_core_dominant_alleles.faa
  --core-dominant-summary   .../Strep_cdhit_core_dominant_summary.tsv

Definition of dominance:
  For each gene (Cluster_ID), pick the allele (New_ID) present in the
  largest number of strains (presence/absence from the sparse matrix).
  Ties break by allele ID (lexicographic).
  dominant_frequency = dominant_count / N_strains_overall  (fraction, not %)
"""

import sys
import os
import re
import argparse
from collections import defaultdict, Counter
from typing import Dict, List, Tuple, Iterable, Set

import numpy as np
import pandas as pd
import sparse_utils
from scipy.sparse import csr_matrix, csc_matrix

ALLELE_RE = re.compile(r'^[A-Za-z]+_C\d+A\d+$')  # e.g., Strep_C123A4
GENE_RE   = re.compile(r'^[A-Za-z]+_C\d+$')      # e.g., Strep_C123

def is_allele_label(s: str) -> bool:
    return bool(ALLELE_RE.match(s or ""))

def is_gene_label(s: str) -> bool:
    return bool(GENE_RE.match(s or ""))

def read_core_genes(path: str) -> Set[str]:
    core = set()
    with open(path, 'r', encoding='utf-8', errors='ignore') as fh:
        for line in fh:
            s = line.strip()
            if s and not s.startswith("#"):
                core.add(s.split()[0])
    return core

def stream_write_selected_fasta(in_faa: str, out_faa: str, wanted: Set[str]):
    """Write only sequences whose header first token matches wanted allele IDs."""
    if not wanted:
        open(out_faa, "w").close()
        return
    n = 0
    write = False
    with open(in_faa, 'r', encoding='utf-8', errors='ignore') as fin, \
         open(out_faa, 'w', encoding='utf-8') as fout:
        for line in fin:
            if line.startswith('>'):
                first = line[1:].strip().split()[0]
                write = (first in wanted)
                if write:
                    n += 1
            if write:
                fout.write(line)
    print(f"  FASTA: wrote {n} sequences to {os.path.basename(out_faa)}", file=sys.stderr)

def main():
    ap = argparse.ArgumentParser(description="Dominant allele extraction using sparse_utils LSDF.")
    ap.add_argument('--allele-matrix', required=True)
    ap.add_argument('--allele-labels', required=True)
    ap.add_argument('--allele-names', required=True)
    ap.add_argument('--allele-faa', required=True)
    ap.add_argument('--core-genes', required=True)

    ap.add_argument('--all-dominant-faa', required=True)
    ap.add_argument('--all-dominant-summary', required=True)
    ap.add_argument('--core-dominant-faa', required=True)
    ap.add_argument('--core-dominant-summary', required=True)

    ap.add_argument('--max-alternatives', type=int, default=3)
    args = ap.parse_args()

    # -------------------------------------------------------------------------
    # 1) Load sparse matrix via sparse_utils (keeps correct labels/orientation)
    # -------------------------------------------------------------------------
    print("Loading LSDF ...", file=sys.stderr)
    lsdf = sparse_utils.read_lsdf(args.allele_matrix, args.allele_labels)
    n_rows, n_cols = lsdf.shape
    print(f"  LSDF shape: {n_rows} x {n_cols}", file=sys.stderr)

    # Detect which axis are alleles vs strains
    # Prefer to identify by label format
    rows_look_like_alleles = any(is_allele_label(x) for x in lsdf.index[: min(1000, len(lsdf.index))])
    cols_look_like_alleles = any(is_allele_label(x) for x in lsdf.columns[: min(1000, len(lsdf.columns))])

    if rows_look_like_alleles and not cols_look_like_alleles:
        orientation = "rows_are_alleles"
    elif cols_look_like_alleles and not rows_look_like_alleles:
        orientation = "cols_are_alleles"
    else:
        # Fallback: if neither clearly matches, assume rows are alleles (historical format)
        orientation = "rows_are_alleles"
        print("  WARNING: Could not clearly detect orientation from labels; assuming rows=alleles.", file=sys.stderr)

    print(f"  Orientation: {orientation}", file=sys.stderr)

    # Binarize presence (treat any nonzero as 1)
    if orientation == "rows_are_alleles":
        M = lsdf.data.tocsr()
        M.data[:] = 1
        allele_labels = lsdf.index
        n_strains = M.shape[1]
        # counts per allele (row nnz)
        allele_counts = M.getnnz(axis=1).astype(np.int64)  # length = #rows
    else:
        M = lsdf.data.tocsc()
        M.data[:] = 1
        allele_labels = lsdf.columns
        n_strains = M.shape[0]
        # counts per allele (col nnz)
        allele_counts = M.getnnz(axis=0).astype(np.int64)  # length = #cols

    print(f"  Total strains: {n_strains}", file=sys.stderr)
    print(f"  Alleles seen on matrix axis: {len(allele_labels)}", file=sys.stderr)

    # Map: allele -> (gene, [WPs...]); cluster_id -> representative WP
    print("Reading allele-names TSV (chunked) ...", file=sys.stderr)
    allele_to_gene: Dict[str, str] = {}
    allele_to_wps: Dict[str, List[str]] = defaultdict(list)
    cluster_to_rep: Dict[str, str] = {}

    # Only keep mapping for alleles that actually appear on the matrix axis
    matrix_alleles_set = set(allele_labels)

    # Some files mark 'representative' as 'True'/'False' or 1/0
    usecols = ['Original_ID', 'New_ID', 'Cluster_ID', 'representative']
    # Allow missing 'representative' column
    try:
        first_chunk = pd.read_csv(args.allele_names, sep='\t', nrows=1)
        if 'representative' not in first_chunk.columns:
            usecols = ['Original_ID', 'New_ID', 'Cluster_ID']
    except Exception:
        pass

    chunksize = 200000
    seen_reps = set()

    for chunk in pd.read_csv(args.allele_names, sep='\t', usecols=usecols, chunksize=chunksize, dtype=str):
        chunk = chunk.fillna('')
        # Normalize representative flag
        rep_mask = False
        if 'representative' in chunk.columns:
            rep_mask = chunk['representative'].astype(str).str.lower().isin(['true', '1', 'yes', 'y'])

        # Restrict to alleles in matrix to keep dicts small
        if 'New_ID' in chunk.columns:
            chunk = chunk[chunk['New_ID'].isin(matrix_alleles_set)]

        for row in chunk.itertuples(index=False):
            if 'representative' in usecols and len(usecols) == 4:
                Original_ID, New_ID, Cluster_ID, representative = row
            else:
                Original_ID, New_ID, Cluster_ID = row
                representative = ''

            if not New_ID:
                continue

            # Map allele -> gene (first one wins; all rows for same New_ID share same Cluster_ID)
            if New_ID not in allele_to_gene:
                allele_to_gene[New_ID] = Cluster_ID

            # Collect WPs per allele (unique)
            if Original_ID and (not allele_to_wps[New_ID] or allele_to_wps[New_ID][-1] != Original_ID):
                # Ensure uniqueness without O(n^2)
                if Original_ID not in allele_to_wps[New_ID]:
                    allele_to_wps[New_ID].append(Original_ID)

            # Cluster representative
            if representative and (str(representative).lower() in ('true', '1', 'yes', 'y')):
                if Cluster_ID and Cluster_ID not in seen_reps:
                    cluster_to_rep[Cluster_ID] = Original_ID
                    seen_reps.add(Cluster_ID)

    print(f"  Mapped alleles -> genes: {len(allele_to_gene)}", file=sys.stderr)
    print(f"  Clusters with explicit representative: {len(cluster_to_rep)}", file=sys.stderr)

    # Build per-gene buckets of (allele_index, allele_id, count)
    print("Grouping alleles by gene ...", file=sys.stderr)
    gene_to_items: Dict[str, List[Tuple[int, str, int]]] = defaultdict(list)

    # Fast map from allele label to its position on the matrix axis
    pos_map = {lab: i for i, lab in enumerate(allele_labels)}

    # Iterate only alleles that have a gene mapping
    matched = 0
    for lab in matrix_alleles_set:
        g = allele_to_gene.get(lab)
        if not g:
            continue
        i = pos_map[lab]
        c = int(allele_counts[i] if orientation == "rows_are_alleles" else allele_counts[i])
        gene_to_items[g].append((i, lab, c))
        matched += 1

    print(f"  Alleles matched to genes: {matched}", file=sys.stderr)
    print(f"  Total genes with ≥1 allele on matrix: {len(gene_to_items)}", file=sys.stderr)

    # Helper to compute union of presence across alleles within a gene
    def strains_with_gene_count(indices: List[int]) -> int:
        if not indices:
            return 0
        if orientation == "rows_are_alleles":
            sub = M[indices, :]               # rows slice, shape (#alleles_in_gene, #strains)
            row_sum = np.asarray(sub.sum(axis=0)).ravel()
            return int((row_sum > 0).sum())
        else:
            sub = M[:, indices]               # columns slice, shape (#strains, #alleles_in_gene)
            col_sum = np.asarray(sub.sum(axis=1)).ravel()
            return int((col_sum > 0).sum())

    # Select dominant per gene and assemble summary rows
    print("Selecting dominant allele per gene ...", file=sys.stderr)
    rows = []
    dominant_ids_all = set()

    for gene, items in gene_to_items.items():
        # items: list of (pos, allele_id, count)
        items.sort(key=lambda t: (-t[2], t[1]))  # by count desc, then allele id asc
        dom_pos, dom_allele, dom_count = items[0]
        indices = [t[0] for t in items]

        total_strains_gene = strains_with_gene_count(indices)
        dom_freq = (dom_count / n_strains) if n_strains > 0 else 0.0

        # dominant proteins (choose 1st WP mapped to that allele if available)
        dom_wps = allele_to_wps.get(dom_allele, [])
        dom_wp_str = dom_wps[0] if dom_wps else ""

        # cluster representative
        rep_wp = cluster_to_rep.get(gene, "")
        if not rep_wp:
            rep_wp = dom_wp_str

        # alternatives: collect non-dominant WPs across other alleles (limited)
        alt_list: List[str] = []
        seen_alt = set([dom_wp_str, rep_wp]) if rep_wp else set([dom_wp_str])
        for _, alt_aid, _ in items[1:]:
            for wp in allele_to_wps.get(alt_aid, []):
                if wp and wp not in seen_alt:
                    alt_list.append(wp)
                    seen_alt.add(wp)
                    if len(alt_list) >= args.max_alternatives:
                        break
            if len(alt_list) >= args.max_alternatives:
                break

        rows.append({
            'gene': gene,
            'dominant_allele': dom_allele,
            'dominant_count': int(dom_count),
            'dominant_frequency': float(dom_freq),          # fraction (not percent)
            'dominant_proteins': dom_wp_str,
            'cluster_representative': rep_wp,
            'alternative_proteins': ",".join(alt_list),
            'total_alleles': int(len(items)),
            'total_strains_with_gene': int(total_strains_gene),
        })
        dominant_ids_all.add(dom_allele)

    # Write ALL summary
    all_df = pd.DataFrame(rows).sort_values('gene').reset_index(drop=True)
    os.makedirs(os.path.dirname(os.path.abspath(args.all_dominant_summary)) or ".", exist_ok=True)
    all_df.to_csv(args.all_dominant_summary, sep='\t', index=False)
    print(f"Wrote ALL summary: {args.all_dominant_summary} ({len(all_df)} rows)", file=sys.stderr)

    # Write ALL FASTA (stream)
    stream_write_selected_fasta(args.allele_faa, args.all_dominant_faa, dominant_ids_all)

    # -------------------------------------------------------------------------
    # 2) Subset to CORE
    # -------------------------------------------------------------------------
    core_set = read_core_genes(args.core_genes)
    print(f"Core gene list: {len(core_set)} IDs", file=sys.stderr)

    if core_set:
        core_df = all_df[all_df['gene'].isin(core_set)].copy()
        core_df.to_csv(args.core_dominant_summary, sep='\t', index=False)
        print(f"Wrote CORE summary: {args.core_dominant_summary} ({len(core_df)} rows)", file=sys.stderr)

        core_dom_ids = set(core_df['dominant_allele'].tolist())
        stream_write_selected_fasta(args.allele_faa, args.core_dominant_faa, core_dom_ids)
    else:
        # Still create empty artifacts with header
        pd.DataFrame(columns=[
            'gene','dominant_allele','dominant_count','dominant_frequency',
            'dominant_proteins','cluster_representative','alternative_proteins',
            'total_alleles','total_strains_with_gene'
        ]).to_csv(args.core_dominant_summary, sep='\t', index=False)
        open(args.core_dominant_faa, 'w').close()
        print("Core list empty → wrote empty CORE outputs.", file=sys.stderr)

    # Quick sanity echo
    if len(core_set) and core_df.empty:
        # help the user see why there’s a mismatch
        head_core = ", ".join(list(sorted(core_set))[:5])
        print("WARNING: 0 core genes matched in ALL summary.", file=sys.stderr)
        print(f"  core head: {head_core}", file=sys.stderr)
        example_gene = next(iter(core_set))
        print(f"  Try grep for that gene in ALL summary: {example_gene}", file=sys.stderr)

    print("Done.", file=sys.stderr)

if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        print(f"FATAL: {e}", file=sys.stderr)
        raise

