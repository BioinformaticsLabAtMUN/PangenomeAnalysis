#!/usr/bin/env python3
"""
Annotate dominant alleles using local NCBI GTF files.

Inputs:
  --input-summary  : <prefix>_all_dominant_summary.tsv  (from extract_dominant_alleles_unified.py)
  --core-summary   : <prefix>_core_dominant_summary.tsv
  --gtf-directory  : directory with *.gtf files (RefSeq)
  --annotate-scope : "core" or "all" (if "all", writes both core & accessory)

Outputs (changed):
  --output-core-annotations            (written)
  --output-core-merged                 (IGNORED, not written)
  --output-accessory-annotations       (written)
  --output-accessory-merged            (IGNORED, not written)
"""

import argparse
import os
import re
import sys
import glob
import pandas as pd
from typing import Dict, List

# ---------- constants ----------
ANNOTATION_COLS = [
    "protein_id", "locus_tag", "gene_id", "product",
    "go_ids", "go_process", "go_function", "go_component"
]

# These columns exist in the dominant summary produced upstream
MERGED_PREFIX_COLS = [
    "gene", "dominant_allele", "dominant_count", "dominant_frequency",
    "dominant_proteins", "cluster_representative", "alternative_proteins",
    "total_alleles", "total_strains_with_gene"
]

WP_RE = re.compile(r'WP_\d+(?:\.\d+)?')


def parse_attributes(attr_field: str) -> Dict[str, List[str]]:
    """Parse GTF attributes (column 9)."""
    out: Dict[str, List[str]] = {}
    for m in re.finditer(r'(\S+)\s+"([^"]*)"', attr_field):
        k = m.group(1)
        v = m.group(2)
        out.setdefault(k, []).append(v)
    return out


def score_record(rec: Dict[str, str]) -> int:
    """Choose 'best' annotation when multiple genomes provide same protein_id."""
    s = 0
    if rec.get("go_ids"):
        s += 2
    if rec.get("go_process") or rec.get("go_function") or rec.get("go_component"):
        s += 2
    prod = rec.get("product", "")
    if prod and "hypothetical" not in prod.lower():
        s += 1
    return s


def build_gtf_index(gtf_dir: str) -> pd.DataFrame:
    """
    Scan *.gtf, collect annotations for protein_id (WP_...), preferring richer records.
    Returns DataFrame with unique protein_id rows.
    """
    mapping: Dict[str, Dict[str, str]] = {}
    gtf_paths = sorted(glob.glob(os.path.join(gtf_dir, "*.gtf")))
    if not gtf_paths:
        print(f"ERROR: No .gtf files found in {gtf_dir}", file=sys.stderr)
        sys.exit(1)

    for path in gtf_paths:
        print(f"Processing {os.path.basename(path)}...")
        try:
            with open(path, "r", encoding="utf-8", errors="ignore") as fh:
                for line in fh:
                    if not line or line.startswith("#"):
                        continue
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 9:
                        continue
                    if cols[2] != "CDS":
                        continue
                    attrs = parse_attributes(cols[8])

                    # protein_id usually present; if not, try db_xref GenBank:WP_...
                    prot = None
                    if "protein_id" in attrs:
                        prot = attrs["protein_id"][0]
                    if not prot and "db_xref" in attrs:
                        for v in attrs["db_xref"]:
                            m = WP_RE.search(v.replace("GenBank:", ""))
                            if m:
                                prot = m.group(0)
                                break
                    if not prot or not WP_RE.match(prot):
                        continue

                    rec = {
                        "protein_id": prot,
                        "locus_tag": attrs.get("locus_tag", [""])[0],
                        "gene_id": attrs.get("gene_id", [""])[0],
                        "product": attrs.get("product", [""])[0],
                        "go_process": "|".join(attrs.get("go_process", [])),
                        "go_function": "|".join(attrs.get("go_function", [])),
                        "go_component": "|".join(attrs.get("go_component", [])),
                        "go_ids": ";".join(sorted({v for v in attrs.get("Ontology_term", []) if v.startswith("GO:")}))
                    }

                    if prot in mapping:
                        if score_record(rec) > score_record(mapping[prot]):
                            mapping[prot] = rec
                    else:
                        mapping[prot] = rec
        except Exception as e:
            print(f"Warning: failed reading {path}: {e}", file=sys.stderr)
            continue

    print(f"Loaded annotations for {len(mapping)} proteins")
    if not mapping:
        return pd.DataFrame(columns=ANNOTATION_COLS)

    return pd.DataFrame(mapping.values(), columns=ANNOTATION_COLS)


def pick_protein_id(row: pd.Series) -> str:
    """
    Decide which protein_id to annotate:
      1) first WP_ in dominant_proteins (comma/space separated)
      2) cluster_representative if WP_
      3) first WP_ in alternative_proteins
    """
    for col in ("dominant_proteins", "cluster_representative", "alternative_proteins"):
        raw = str(row.get(col, "") or "")
        for tok in re.split(r"[,\s]+", raw):
            if WP_RE.match(tok):
                return tok
    return ""


def ensure_tsv(path: str, df: pd.DataFrame):
    """Write a TSV with headers; if df is empty, still write headers."""
    os.makedirs(os.path.dirname(os.path.abspath(path)) or ".", exist_ok=True)
    df.to_csv(path, sep="\t", index=False)


def annotate_block(block_df: pd.DataFrame, ann: pd.DataFrame) -> pd.DataFrame:
    """
    Return the annotations-only table:
      columns = ['gene','dominant_allele'] + ANNOTATION_COLS
    """
    if block_df.empty:
        return pd.DataFrame(columns=["gene", "dominant_allele"] + ANNOTATION_COLS)

    # guard missing expected columns
    for c in MERGED_PREFIX_COLS:
        if c not in block_df.columns:
            block_df[c] = ""

    work = block_df.copy()
    work["protein_id"] = work.apply(pick_protein_id, axis=1)

    merged = work.merge(ann, on="protein_id", how="left")

    # annotations-only view (one row per gene), now includes dominant_allele
    ann_cols = ["gene", "dominant_allele"] + ANNOTATION_COLS
    ann_only = merged[ann_cols].copy()
    return ann_only


def main():
    p = argparse.ArgumentParser(description="Annotate dominant alleles from NCBI GTFs")
    p.add_argument("--input-summary", required=True)
    p.add_argument("--core-summary", required=True)
    p.add_argument("--gtf-directory", required=True)
    p.add_argument("--annotate-scope", choices=["core", "all"], default="all")
    p.add_argument("--output-core-annotations", required=True)
    p.add_argument("--output-core-merged", required=True)           # kept for API compatibility (ignored)
    p.add_argument("--output-accessory-annotations", required=True)
    p.add_argument("--output-accessory-merged", required=True)      # kept for API compatibility (ignored)
    args = p.parse_args()

    # load summaries
    all_df = pd.read_csv(args.input_summary, sep="\t")
    core_df = pd.read_csv(args.core_summary, sep="\t") if os.path.exists(args.core_summary) else pd.DataFrame(columns=["gene"])
    core_genes = set(core_df["gene"].tolist()) if "gene" in core_df.columns else set()

    print(f"Input genes: {len(all_df)}")
    print(f"Core genes:  {len(core_genes)}")

    # split into core/accessory if needed
    if args.annotate_scope == "core":
        accessory_df = pd.DataFrame(columns=all_df.columns)
    else:
        accessory_df = all_df[~all_df["gene"].isin(core_genes)].copy()

    # index annotations from GTFs
    ann = build_gtf_index(args.gtf_directory)

    # annotate core
    core_ann_only = annotate_block(core_df, ann)
    ensure_tsv(args.output_core_annotations, core_ann_only)

    # annotate accessory
    acc_ann_only = annotate_block(accessory_df, ann)
    ensure_tsv(args.output_accessory_annotations, acc_ann_only)

    # do NOT write merged files (explicitly requested)
    # touch nothing for args.output_core_merged / args.output_accessory_merged

    # tiny summary
    core_hits = int(core_ann_only["protein_id"].notna().sum()) if not core_ann_only.empty else 0
    acc_hits = int(acc_ann_only["protein_id"].notna().sum()) if not acc_ann_only.empty else 0
    print("\nResults:")
    print(f"  Core genes annotated:      {core_hits}/{len(core_ann_only)}")
    print(f"  Accessory genes annotated: {acc_hits}/{len(acc_ann_only)}")
    print("\nNote: *_genes_with_annotations.tsv outputs are intentionally suppressed.")
    

if __name__ == "__main__":
    main()

