#!/usr/bin/env python3
"""
Cascading UniProt Annotator (Batched, Bulletproof)
- Input IDs are UniProt accessions (no mapping step)
- Uses /uniprotkb/search with fields=... (TSV) + cursor pagination (max 500/page)
- Validates responses before parsing; dumps non-TSV pages for inspection
- Exponential backoff on 429/5xx; gzip supported; output schema unchanged
"""

import argparse
import os
import re
import json
import time
import random
import threading
import gzip
from pathlib import Path
from io import StringIO, BytesIO

import pandas as pd
import requests

UNI_BASE = "https://rest.uniprot.org/uniprotkb/search"
UNI_FIELDS = "accession,id,protein_name,go_p,go_f,go_c,organism_name,gene_names"
UNI_PAGE_MAX = 500  # UniProt hard limit per page

class CascadingUniProtAnnotator:
    def __init__(self, cache_dir=None, delay=0.1, max_workers=6):
        self.delay = max(0.0, delay)
        self.max_workers = max_workers  # kept for CLI parity
        self.session = requests.Session()
        self.session.headers.update({
            "User-Agent": "PangenomeAnnotator/3.1 (research; contact: you@example.org)",
            "Accept": "text/tab-separated-values, text/plain; q=0.9, */*; q=0.1",
            "Accept-Encoding": "gzip",
        })
        self.lock = threading.Lock()
        self.total_api_calls = 0

        self.cache_dir = Path(cache_dir) if cache_dir else None
        if self.cache_dir:
            self.cache_dir.mkdir(parents=True, exist_ok=True)

        self._memo_annotations = {}

    # -------- cache helpers --------
    def _cache_path(self, acc): return self.cache_dir / f"{acc}.json" if self.cache_dir else None

    def _cache_get(self, acc):
        if not self.cache_dir: return None
        p = self._cache_path(acc)
        if not p or not p.exists(): return None
        try:
            return json.loads(p.read_text())
        except Exception:
            return None

    def _cache_put(self, acc, payload):
        if not self.cache_dir: return
        try:
            self._cache_path(acc).write_text(json.dumps(payload))
        except Exception:
            pass

    # -------- HTTP with retry/backoff --------
    def _retry_get(self, url, params, max_tries=6, timeout=60):
        last = None
        for attempt in range(1, max_tries + 1):
            with self.lock:
                self.total_api_calls += 1
            last = self.session.get(url, params=params, timeout=timeout)
            if last.status_code < 400:
                return last
            if last.status_code in (429, 500, 502, 503, 504):
                sleep_for = min(60, (2 ** attempt) * (0.5 + random.random()))
                print(f"[WARN] HTTP {last.status_code}; backoff {sleep_for:.1f}s")
                time.sleep(sleep_for)
                continue
            return last
        return last

    # -------- Batched UniProt fetch (cursor pagination) --------
    def fetch_batch(self, protein_ids, desc=None):
        if not protein_ids:
            return {}

        # sanitize + de-dup + honor memo/cache
        ids = []
        for pid in map(str, protein_ids):
            pid = pid.strip()
            if not pid:
                continue
            if pid in self._memo_annotations:
                continue
            cached = self._cache_get(pid)
            if cached is not None:
                self._memo_annotations[pid] = cached
                continue
            # Skip only truly invalid entries (empty or whitespace)
            # All non-empty strings are treated as potential UniProt IDs
            ids.append(pid)

        if not ids:
            return {k: v for k, v in self._memo_annotations.items() if k in protein_ids}

        out = {}
        chunk_size = 400  # fits query length; UniProt pages will still return 500 max/page

        for i in range(0, len(ids), chunk_size):
            chunk = ids[i:i + chunk_size]
            if not chunk:
                continue

            query = " OR ".join([f"accession:{acc}" for acc in chunk])
            params = {
                "query": f"({query})",
                "format": "tsv",
                "fields": UNI_FIELDS,
                "size": str(UNI_PAGE_MAX),
                "compressed": "true",
            }

            cursor = None
            page_no = 0
            while True:
                if cursor:
                    params["cursor"] = cursor

                r = self._retry_get(UNI_BASE, params)

                # Decompress response if compressed parameter was used
                response_text = r.text
                if params.get("compressed") == "true" and r.status_code == 200:
                    try:
                        # The response content is gzipped, decompress it
                        decompressed_bytes = gzip.decompress(r.content)
                        response_text = decompressed_bytes.decode('utf-8')
                        print(f"[DEBUG] Successfully decompressed response (from {len(r.content)} bytes to {len(response_text)} chars)")
                    except Exception as e:
                        print(f"[WARN] Failed to decompress response: {e}. Trying as plain text.")
                        response_text = r.text

                # debug headers (one-liner)
                ct = (r.headers.get("Content-Type", "") or "").lower()
                link_hdr = r.headers.get("Link", "") or ""
                print(f"[DEBUG] {desc or 'batch'} chunk={i//chunk_size+1} page={page_no} "
                      f"HTTP={r.status_code} CT={ct} len={len(response_text)} next={('rel=\"next\"' in link_hdr)}")

                if r.status_code != 200 or not response_text:
                    print(f"[WARN] Bad page (HTTP {r.status_code}); skipping.")
                    break

                # fast smoke-test: TSV header must contain 'Entry' and tabs
                lines = response_text.splitlines()
                if not lines or ("\t" not in lines[0]) or ("Entry" not in lines[0]):
                    # dump sample for inspection
                    import tempfile
                    dbg = os.path.join(tempfile.gettempdir(), f"uniprot_bad_page_{int(time.time())}.txt")
                    with open(dbg, "w", encoding="utf-8") as fh:
                        fh.write(response_text[:5000])
                    print(f"[WARN] Non-TSV or unexpected header: {lines[0][:50] if lines else '<empty>'}. "
                          f"Saved sample: {dbg}. Skipping page.")
                    break

                # parse TSV safely
                try:
                    import csv
                    df = pd.read_csv(StringIO(response_text), sep="\t", dtype=str,
                                     engine="python", quoting=csv.QUOTE_MINIMAL)
                except Exception as e:
                    print(f"[WARN] pandas failed to parse TSV ({e}). Skipping page.")
                    break

                for _, row in df.iterrows():
                    acc = row.get("Entry")
                    if not isinstance(acc, str) or not acc:
                        continue
                    ann = {
                        "Entry": row.get("Entry", acc),
                        "Entry Name": row.get("Entry Name", ""),
                        "Protein names": row.get("Protein names", ""),
                        "Gene Ontology (biological process)": row.get("Gene Ontology (biological process)", ""),
                        "Gene Ontology (molecular function)": row.get("Gene Ontology (molecular function)", ""),
                        "Gene Ontology (cellular component)": row.get("Gene Ontology (cellular component)", ""),
                        "Organism": row.get("Organism", ""),
                        "Gene Names": row.get("Gene Names", ""),
                    }
                    self._memo_annotations[acc] = ann
                    self._cache_put(acc, ann)
                    out[acc] = ann

                # cursor pagination via Link header
                if 'rel="next"' in link_hdr:
                    m = re.search(r'[?&]cursor=([^&>]+)', link_hdr)
                    cursor = m.group(1) if m else None
                    page_no += 1
                    if not cursor:
                        break
                    continue
                break

        # return only what the caller asked for
        return {k: self._memo_annotations.get(k) for k in protein_ids if k in self._memo_annotations}

    # -------- scoring & checks --------
    @staticmethod
    def has_go_annotations(ann):
        if not ann: return False
        for k in ("Gene Ontology (biological process)",
                  "Gene Ontology (molecular function)",
                  "Gene Ontology (cellular component)"):
            v = ann.get(k, "")
            if isinstance(v, str) and v.strip():
                return True
        return False

    def score(self, ann):
        if not ann: return 0
        s = 1
        if str(ann.get("Protein names", "")).strip(): s += 2
        # Give GO annotations much higher weight
        if self.has_go_annotations(ann): s += 100  # Changed from 10 to 100
        return s

    # -------- cascade orchestrators --------
    def annotate_genes_cascading(self, df, max_alts=4, is_accessory=False, chunk_size=5000,
                                 sample_size=None, min_go_rate=0.7):
        if sample_size and len(df) > sample_size:
            print(f"Sampling {sample_size} genes from {len(df)} total")
            df = df.sample(n=sample_size, random_state=42)

        print("\nCascading UniProt Annotation")
        print(f"Total genes: {len(df)}")
        print(f"Max alternatives: {max_alts}")
        if is_accessory and len(df) > 10000:
            print(f"Chunk size: {chunk_size}")
            print(f"Early stopping GO rate: {min_go_rate*100:.0f}%")
        print("="*60)

        if is_accessory and len(df) > chunk_size:
            return self._annotate_chunked(df, max_alts, chunk_size, min_go_rate)
        return self._annotate_standard(df, max_alts)

    def _annotate_standard(self, df, max_alts):
        start = time.time()
        out = {}

        # Round 1: dominant
        dom_map, dom_ids = {}, []
        for idx, row in df.iterrows():
            gene = row["gene"]
            dominant = str(row.get("dominant_proteins", "")).strip()
            if dominant and dominant.lower() != "nan":
                pid = dominant.split(",")[0].strip()
                if pid:
                    dom_map[gene] = (pid, idx)
                    dom_ids.append(pid)

        dom_ids = list(set(dom_ids))
        if not dom_ids:
            print("Round 1: no dominant accessions; skipping")
            dom_res = {}
        else:
            dom_res = self.fetch_batch(dom_ids, desc="Dominant")

        annotated, go_hits = set(), 0
        for gene, (pid, idx) in dom_map.items():
            ann = dom_res.get(pid)
            if ann:
                out[gene] = {"annotation": ann, "protein_id": pid, "source": "dominant_protein", "idx": idx}
                annotated.add(gene)
                if self.has_go_annotations(ann): go_hits += 1

        total = len(df)
        print(f"Round 1: {len(annotated)}/{total} ({(len(annotated)/max(1,total))*100:.1f}%) "
              f"with GO: {go_hits} ({(go_hits/max(1,len(annotated)))*100:.1f}%)")

        # Round 2: representatives for those without GO
        need_go = [g for g in out if not self.has_go_annotations(out[g]["annotation"])]
        if need_go:
            rep_ids, rep_for = [], {}
            for g in need_go:
                idx = out[g]["idx"]
                row = df.loc[idx]
                rep = str(row.get("cluster_representative", "")).strip()
                if rep and rep.lower() != "nan":
                    rep_ids.append(rep)
                    rep_for[g] = rep
            if rep_ids:
                rep_res = self.fetch_batch(list(set(rep_ids)), desc="Representatives")
                improved = 0
                go_found = 0
                for g, rep in rep_for.items():
                    ann = rep_res.get(rep)
                    if ann:
                        # Prioritize GO annotations heavily
                        if self.has_go_annotations(ann):
                            out[g] = {"annotation": ann, "protein_id": rep, "source": "cluster_representative",
                                      "idx": out[g]["idx"]}
                            improved += 1
                            go_found += 1
                        elif self.score(ann) > self.score(out[g]["annotation"]):
                            out[g] = {"annotation": ann, "protein_id": rep, "source": "cluster_representative",
                                      "idx": out[g]["idx"]}
                            improved += 1
                print(f"Round 2: improved {improved} genes ({go_found} with GO)")

        # Rounds 3+: alternatives (position r)
        for r in range(max_alts):
            still = [g for g in out if not self.has_go_annotations(out[g]["annotation"])]
            if not still: break
            alt_ids, alt_for = [], {}
            for g in still:
                idx = out[g]["idx"]
                row = df.loc[idx]
                alts = str(row.get("alternative_proteins", "")).strip()
                if not alts or alts.lower() == "nan": continue
                arr = [x.strip() for x in alts.split(",") if x.strip()]
                if len(arr) > r:
                    pid = arr[r]
                    alt_ids.append(pid)
                    alt_for[g] = pid
            if not alt_ids:
                print(f"Round {3+r}: no alternatives at position {r+1}")
                continue
            alt_res = self.fetch_batch(list(set(alt_ids)), desc=f"Alt#{r+1}")
            improved = 0
            go_found = 0
            for g, pid in alt_for.items():
                ann = alt_res.get(pid)
                if ann:
                    # Always accept if it has GO annotations
                    if self.has_go_annotations(ann):
                        out[g] = {"annotation": ann, "protein_id": pid, "source": f"alternative_protein_{r+1}",
                                  "idx": out[g]["idx"]}
                        improved += 1
                        go_found += 1
                    # Otherwise only accept if it improves the score (and we still don't have GO)
                    elif not self.has_go_annotations(out[g]["annotation"]) and self.score(ann) > self.score(out[g]["annotation"]):
                        out[g] = {"annotation": ann, "protein_id": pid, "source": f"alternative_protein_{r+1}",
                                  "idx": out[g]["idx"]}
                        improved += 1
            print(f"Round {3+r}: improved {improved} genes ({go_found} with GO)")

        elapsed = time.time() - start
        with_go = sum(1 for v in out.values() if self.has_go_annotations(v["annotation"]))
        print("="*60)
        print(f"Cascading Annotation Completed in {elapsed:.1f}s; API calls: {self.total_api_calls}")
        print(f"Annotated: {len(out)}/{len(df)} ({(len(out)/max(1,len(df)))*100:.1f}%)")
        print(f"With GO:  {with_go}/{len(df)} ({(with_go/max(1,len(df)))*100:.1f}%)")
        return out

    def _annotate_chunked(self, df, max_alts, chunk_size, min_go_rate):
        start = time.time()
        out = {}
        n = len(df)
        chunks = (n + chunk_size - 1) // chunk_size
        print(f"\nProcessing {n} genes in {chunks} chunks of {chunk_size}")

        for c in range(chunks):
            s = c * chunk_size
            e = min((c + 1) * chunk_size, n)
            cdf = df.iloc[s:e]
            print(f"\nChunk {c+1}/{chunks}: genes {s+1}-{e}")

            # Round 1: dominant
            dom_map, dom_ids = {}, []
            for idx, row in cdf.iterrows():
                gene = row["gene"]
                dominant = str(row.get("dominant_proteins", "")).strip()
                if dominant and dominant.lower() != "nan":
                    pid = dominant.split(",")[0].strip()
                    if pid:
                        dom_map[gene] = (pid, idx)
                        dom_ids.append(pid)

            dom_ids = list(set(dom_ids))
            dom_res = {}
            if not dom_ids:
                print(f"  Chunk {c+1}: no dominant accessions; skipping reps/alts")
            else:
                dom_res = self.fetch_batch(dom_ids, desc=f"Chunk {c+1} Dominant")

            ann_ct, go_ct = 0, 0
            for gene, (pid, idx) in dom_map.items():
                ann = dom_res.get(pid)
                if ann:
                    out[gene] = {"annotation": ann, "protein_id": pid, "source": "dominant_protein", "idx": idx}
                    ann_ct += 1
                    if self.has_go_annotations(ann): go_ct += 1

            go_rate = (go_ct / len(cdf)) if len(cdf) else 0
            print(f"  Annotated {ann_ct}/{len(cdf)} | GO rate: {go_rate*100:.1f}%")

            if go_rate >= min_go_rate:
                print(f"  GO good ({go_rate*100:.1f}%) -> Skip reps/alts for this chunk")
                continue

            # Round 2: representatives (only for genes lacking GO)
            genes_no_go = [g for g in dom_map if g in out and not self.has_go_annotations(out[g]["annotation"])]
            rep_ids, rep_for = [], {}
            for g in genes_no_go:
                idx = out[g]["idx"]
                row = df.loc[idx]
                rep = str(row.get("cluster_representative", "")).strip()
                if rep and rep.lower() != "nan":
                    rep_ids.append(rep)
                    rep_for[g] = rep

            if rep_ids:
                rep_res = self.fetch_batch(list(set(rep_ids)), desc=f"Chunk {c+1} Reps")
                improved = 0
                go_found = 0
                for g, rep in rep_for.items():
                    ann = rep_res.get(rep)
                    if ann:
                        # Prioritize GO annotations heavily
                        if self.has_go_annotations(ann):
                            out[g] = {"annotation": ann, "protein_id": rep, "source": "cluster_representative",
                                      "idx": out[g]["idx"]}
                            improved += 1
                            go_found += 1
                        elif self.score(ann) > self.score(out[g]["annotation"]):
                            out[g] = {"annotation": ann, "protein_id": rep, "source": "cluster_representative",
                                      "idx": out[g]["idx"]}
                            improved += 1
                print(f"  Representatives improved {improved} ({go_found} with GO)")

            # Alternatives
            for r in range(max_alts):
                unresolved = [g for g in dom_map if g in out and not self.has_go_annotations(out[g]["annotation"])]
                if not unresolved:
                    break
                alt_ids, alt_for = [], {}
                for g in unresolved:
                    idx = out[g]["idx"]
                    row = df.loc[idx]
                    alts = str(row.get("alternative_proteins", "")).strip()
                    if not alts or alts.lower() == "nan": continue
                    arr = [x.strip() for x in alts.split(",") if x.strip()]
                    if len(arr) > r:
                        pid = arr[r]
                        alt_ids.append(pid)
                        alt_for[g] = pid
                if not alt_ids:
                    print(f"  Alt #{r+1}: none available")
                    continue
                alt_res = self.fetch_batch(list(set(alt_ids)), desc=f"Chunk {c+1} Alt#{r+1}")
                improved = 0
                go_found = 0
                for g, pid in alt_for.items():
                    ann = alt_res.get(pid)
                    if ann:
                        # Always accept if it has GO annotations
                        if self.has_go_annotations(ann):
                            out[g] = {"annotation": ann, "protein_id": pid, "source": f"alternative_protein_{r+1}",
                                      "idx": out[g]["idx"]}
                            improved += 1
                            go_found += 1
                        # Otherwise only accept if it improves the score (and we still don't have GO)
                        elif not self.has_go_annotations(out[g]["annotation"]) and self.score(ann) > self.score(out[g]["annotation"]):
                            out[g] = {"annotation": ann, "protein_id": pid, "source": f"alternative_protein_{r+1}",
                                      "idx": out[g]["idx"]}
                            improved += 1
                print(f"  Alt #{r+1}: improved {improved} ({go_found} with GO)")

        elapsed = time.time() - start
        with_go = sum(1 for v in out.values() if self.has_go_annotations(v["annotation"]))
        print("="*60)
        print(f"Chunked Annotation Completed in {elapsed:.1f}s; API calls: {self.total_api_calls}")
        print(f"Annotated: {len(out)}/{len(df)} ({(len(out)/max(1,len(df)))*100:.1f}%)")
        print(f"With GO:  {with_go}/{len(df)} ({(with_go/max(1,len(df)))*100:.1f}%)")
        return out

# -------- output helpers --------
def create_annotations_from_cascading(gene_annotations, summary_df, output_file):
    res = summary_df.copy()
    for col in ["protein_id","locus_tag","gene_id","product","go_ids",
                "go_process","go_function","go_component","annotation_source"]:
        if col not in res.columns: res[col] = ""

    for gene, payload in gene_annotations.items():
        idx = payload["idx"]
        ann = payload["annotation"]
        res.at[idx, "protein_id"] = payload["protein_id"]
        res.at[idx, "annotation_source"] = payload["source"]
        res.at[idx, "product"] = ann.get("Protein names", "")

        gp = ann.get("Gene Ontology (biological process)", "")
        gf = ann.get("Gene Ontology (molecular function)", "")
        gc = ann.get("Gene Ontology (cellular component)", "")
        res.at[idx, "go_process"] = gp
        res.at[idx, "go_function"] = gf
        res.at[idx, "go_component"] = gc

        go_ids = []
        for field in (gp, gf, gc):
            if isinstance(field, str) and field:
                go_ids.extend(re.findall(r"GO:\d+", field))
        if go_ids:
            seen = set()
            ordered = [x for x in go_ids if not (x in seen or seen.add(x))]
            res.at[idx, "go_ids"] = ";".join(ordered)

    res.to_csv(output_file, sep="\t", index=False)
    print(f"Annotations saved: {output_file}")
    return res

def separate_core_accessory_genes(all_summary_file, core_summary_file):
    all_df = pd.read_csv(all_summary_file, sep="\t")
    core_df = pd.read_csv(core_summary_file, sep="\t")
    core_set = set(core_df["gene"].tolist())
    accessory_df = all_df[~all_df["gene"].isin(core_set)]
    print("Gene separation:")
    print(f"  - Total genes: {len(all_df)}")
    print(f"  - Core genes: {len(core_df)}")
    print(f"  - Accessory genes: {len(accessory_df)}")
    return core_df, accessory_df

# -------- CLI --------
def main():
    p = argparse.ArgumentParser(description="Cascading UniProt annotation - optimized for large-scale")
    p.add_argument('--input-summary', required=True)
    p.add_argument('--core-summary', required=True)
    p.add_argument('--annotate-scope', choices=['core','all'], default='all')
    p.add_argument('--output-core-annotations', required=True)
    p.add_argument('--output-core-merged', required=True)
    p.add_argument('--output-accessory-annotations', required=True)
    p.add_argument('--output-accessory-merged', required=True)
    p.add_argument('--max-alternatives', type=int, default=4)
    p.add_argument('--delay', type=float, default=0.1)
    p.add_argument('--workers', type=int, default=6)
    p.add_argument('--chunk-size', type=int, default=5000)
    p.add_argument('--sample-size', type=int, default=None)
    p.add_argument('--min-go-rate', type=float, default=0.7)
    args = p.parse_args()

    print("Cascading UniProt Annotation - Batched")
    print("="*60)
    print(f"Annotation scope: {args.annotate_scope}")
    print(f"Max alternatives: {args.max_alternatives}")
    print(f"Workers (unused): {args.workers}  # batching avoids per-ID threads")
    print(f"Chunk size: {args.chunk_size}")
    if args.sample_size:
        print(f"Sample size: {args.sample_size}")
    print(f"Min GO rate for early stopping: {args.min_go_rate*100:.0f}%")

    cache_dir = os.path.join(os.path.dirname(args.output_core_annotations), '.uniprot_cache_cascading')

    core_df, accessory_df = separate_core_accessory_genes(args.input_summary, args.core_summary)
    annotator = CascadingUniProtAnnotator(cache_dir=cache_dir, delay=args.delay, max_workers=args.workers)

    print("\nProcessing Core Genes")
    print("-"*40)
    core_ann = annotator.annotate_genes_cascading(core_df, args.max_alternatives, is_accessory=False)
    create_annotations_from_cascading(core_ann, core_df, args.output_core_annotations)
    pd.read_csv(args.output_core_annotations, sep="\t").to_csv(args.output_core_merged, sep="\t", index=False)
    print(f"Core merged output: {args.output_core_merged}")

    if args.annotate_scope == "all" and len(accessory_df) > 0:
        print("\nProcessing Accessory Genes")
        print("-"*40)
        acc_ann = annotator.annotate_genes_cascading(
            accessory_df, args.max_alternatives, is_accessory=True,
            chunk_size=args.chunk_size, sample_size=args.sample_size, min_go_rate=args.min_go_rate
        )
        create_annotations_from_cascading(acc_ann, accessory_df, args.output_accessory_annotations)
        pd.read_csv(args.output_accessory_annotations, sep="\t").to_csv(args.output_accessory_merged, sep="\t", index=False)
        print(f"Accessory merged output: {args.output_accessory_merged}")
    else:
        print("\nSkipping accessory genes")
        pd.DataFrame().to_csv(args.output_accessory_annotations, sep="\t", index=False)
        pd.DataFrame().to_csv(args.output_accessory_merged, sep="\t", index=False)

    print("\nDone. Optimized for large-scale gene sets.")

if __name__ == "__main__":
    main()
