#!/usr/bin/env python3
"""
GO term clustering using GOATools with Resnik similarity and custom Information Content.
Replicates Revigo's approach but uses your Streptomyces term frequencies instead of database statistics.

This version:
- Keeps the original "combined" outputs for compatibility
- Also writes separate Core and Accessory outputs:
    <prefix>_core_goatools_{bp,mf}_clusters.tsv
    <prefix>_accessory_goatools_{bp,mf}_clusters.tsv
"""

import pandas as pd
import numpy as np
import argparse
from collections import Counter, defaultdict
import sys
import os
import re

# -------------------------
# Optional dependencies
# -------------------------
try:
    from goatools.obo_parser import GODag
    from goatools.semantic import TermCounts, get_info_content, resnik_sim, lin_sim
    from goatools.godag.go_tasks import get_go2ancestors, get_go2descendants
    GOATOOLS_AVAILABLE = True
except ImportError:
    print("ERROR: GOATools not installed. Install with: pip install goatools")
    GOATOOLS_AVAILABLE = False

try:
    from sklearn.cluster import AgglomerativeClustering
    SKLEARN_AVAILABLE = True
except ImportError:
    print("WARNING: scikit-learn not available. Will use simpler clustering.")
    SKLEARN_AVAILABLE = False

# -------------------------
# Settings / helpers
# -------------------------

RELNS = {'is_a', 'part_of'}  # Revigo-like

def go_ancestors_map(go_dag, relns=RELNS):
    """Return dict GO->set(ancestors) with API-compat across GOATOOLS versions."""
    terms = list(go_dag.values())  # GOTerm objects, not strings
    try:
        return get_go2ancestors(terms, relns)
    except TypeError:
        return get_go2ancestors(terms)

def agglomerative_precomputed(distance_matrix, n_clusters, linkage):
    """Compat wrapper for sklearn AgglomerativeClustering metric/affinity param."""
    try:
        model = AgglomerativeClustering(
            n_clusters=n_clusters, metric='precomputed', linkage=linkage
        )
    except TypeError:
        model = AgglomerativeClustering(
            n_clusters=n_clusters, affinity='precomputed', linkage=linkage
        )
    return model.fit_predict(distance_matrix)

# -------------------------
# Parsing / tallying
# -------------------------

def extract_go_ids(cell):
    """Return a set of normalized GO IDs (GO:XXXXXXX) from a mixed-format cell."""
    if pd.isna(cell):
        return set()
    s = str(cell).strip()
    if not s or s.lower() == 'nan':
        return set()
    ids = set(re.findall(r'(?:GO:)?(\d{7})', s))
    return {f"GO:{x}" for x in ids}

def tally_by_namespace(df, go_dag, dedup_within_row=True):
    """
    Build frequency dicts using ontology-specific columns:
      - go_process    -> BP
      - go_function   -> MF
      - go_component  -> CC
    """
    bp_freq = Counter(); mf_freq = Counter(); cc_freq = Counter()

    for _, row in df.iterrows():
        bp_ids = extract_go_ids(row.get('go_process', ''))
        mf_ids = extract_go_ids(row.get('go_function', ''))
        cc_ids = extract_go_ids(row.get('go_component', ''))
        if dedup_within_row:
            bp_ids, mf_ids, cc_ids = set(bp_ids), set(mf_ids), set(cc_ids)
        bp_freq.update(bp_ids); mf_freq.update(mf_ids); cc_freq.update(cc_ids)

    def _ns_filter(freq, ns):
        out = {}
        for go, c in freq.items():
            if go in go_dag and (getattr(go_dag[go], 'namespace', '') or '').lower() == ns:
                out[go] = c
        return out

    bp_freq = _ns_filter(bp_freq, 'biological_process')
    mf_freq = _ns_filter(mf_freq, 'molecular_function')
    cc_freq = _ns_filter(cc_freq, 'cellular_component')

    bp_names = {go: go_dag[go].name for go in bp_freq}
    mf_names = {go: go_dag[go].name for go in mf_freq}
    cc_names = {go: go_dag[go].name for go in cc_freq}
    return bp_freq, mf_freq, cc_freq, bp_names, mf_names, cc_names

# -------------------------
# OBO handling
# -------------------------

def download_go_obo(force_download=False):
    """Download GO ontology file if not present"""
    obo_file = "go-basic.obo"
    if not os.path.exists(obo_file) or force_download:
        print("Downloading GO ontology...")
        try:
            import urllib.request
            url = "https://purl.obolibrary.org/obo/go/go-basic.obo"
            urllib.request.urlretrieve(url, obo_file)
            print(f"Downloaded {obo_file}")
        except Exception as e:
            print(f"Download failed: {e}")
            print("Please download manually from: https://purl.obolibrary.org/obo/go/go-basic.obo")
            sys.exit(1)
    return obo_file

# -------------------------
# IC + Resnik similarity
# -------------------------

def calculate_custom_information_content(go_ids, term_frequencies, go_dag, ontology='BP'):
    """Calculate Information Content from your term frequencies (like Revigo but with your data)."""
    print(f"Calculating custom Information Content for {ontology}...")

    namespace_map = {
        'BP': 'biological_process',
        'MF': 'molecular_function',
        'CC': 'cellular_component',
        'biological_process': 'biological_process',
        'molecular_function': 'molecular_function',
        'cellular_component': 'cellular_component',
    }
    ont = namespace_map.get(ontology, ontology).lower()

    # Filter GO IDs by ontology
    filtered_go_ids = [go for go in go_ids
                       if go in go_dag and (go_dag[go].namespace or '').lower() == ont]
    print(f"Found {len(filtered_go_ids)} {ont} terms")
    if not filtered_go_ids:
        return {}

    # Frequencies from your data
    go_frequencies = {go_id: int(term_frequencies.get(go_id, 1)) for go_id in filtered_go_ids}

    total_freq = float(sum(go_frequencies.values()))
    if total_freq <= 0:
        return {}

    # IC(t) = -log( freq(t) / total_freq )
    ic_values = {}
    for go_id, freq in go_frequencies.items():
        prob = max(1e-15, freq / total_freq)
        ic_values[go_id] = -np.log(prob)

    # Propagate IC values up the hierarchy (ancestors less specific => not greater IC)
    go2anc = go_ancestors_map(go_dag)
    for go_id in filtered_go_ids:
        ancestors = go2anc.get(go_id, set())
        current_ic = ic_values.get(go_id, 0.0)
        for anc in ancestors:
            if anc in ic_values:
                ic_values[anc] = min(ic_values[anc], current_ic * 0.9)
            else:
                ic_values[anc] = current_ic * 0.9

    print(f"Calculated IC for {len(ic_values)} terms")
    print(f"IC range: {min(ic_values.values()):.3f} - {max(ic_values.values()):.3f}")
    return ic_values

def calculate_resnik_similarity_matrix(go_ids, ic_values, go_dag):
    """Calculate Resnik similarity matrix between all GO term pairs."""
    print(f"Calculating Resnik similarity matrix for {len(go_ids)} terms...")
    go_list = sorted(list(go_ids))
    n_terms = len(go_list)
    sim_matrix = np.zeros((n_terms, n_terms), dtype=float)

    go2anc = go_ancestors_map(go_dag)

    for i in range(n_terms):
        go1 = go_list[i]
        anc1 = go2anc.get(go1, set()) | {go1}
        ic1 = ic_values.get(go1, 0.0)
        sim_matrix[i, i] = ic1  # self similarity

        for j in range(i + 1, n_terms):
            go2 = go_list[j]
            anc2 = go2anc.get(go2, set()) | {go2}
            common = anc1 & anc2
            mica_ic = max((ic_values.get(a, 0.0) for a in common), default=0.0)
            sim_matrix[i, j] = mica_ic
            sim_matrix[j, i] = mica_ic

    print(f"Similarity matrix calculated. Range: {sim_matrix.min():.3f} - {sim_matrix.max():.3f}")
    return sim_matrix, go_list

# -------------------------
# Clustering
# -------------------------

def cluster_similarity_matrix(sim_matrix, go_list, term_frequencies, threshold=0.4, method='average'):
    """Cluster GO terms based on similarity matrix."""
    n = len(go_list)
    print(f"Clustering {n} terms with threshold {threshold}...")

    if n == 0:
        print("No terms; skipping clustering.")
        return []

    sim_max = float(sim_matrix.max())
    if sim_max <= 0:
        print("Similarity matrix is all zeros; skipping clustering.")
        return []

    # Convert similarity to distance in [0,1]; keep self-distance = 0
    distance_matrix = 1.0 - (sim_matrix / sim_max)
    np.fill_diagonal(distance_matrix, 0.0)

    if SKLEARN_AVAILABLE and n >= 2:
        max_clusters = max(2, int(n * (1 - threshold)))
        n_clusters = max(2, min(max_clusters, n))  # ensure 2..n
        linkage = method if method in {'average', 'complete', 'single'} else 'average'
        cluster_labels = agglomerative_precomputed(distance_matrix, n_clusters, linkage)
    else:
        print("Using simple threshold-based clustering...")
        cluster_labels = [-1] * n
        cluster_id = 0
        assigned = set()
        for i in range(n):
            if i in assigned:
                continue
            current = [i]
            assigned.add(i)
            for j in range(i + 1, n):
                if j not in assigned and sim_matrix[i, j] > threshold:
                    current.append(j)
                    assigned.add(j)
            for idx in current:
                cluster_labels[idx] = cluster_id
            cluster_id += 1
        for i in range(n):
            if cluster_labels[i] == -1:
                cluster_labels[i] = cluster_id
                cluster_id += 1

    # Group terms by cluster and select representatives (highest frequency)
    clusters = defaultdict(list)
    total_freq = max(1, sum(int(term_frequencies.get(go_id, 1)) for go_id in go_list))
    for idx, lab in enumerate(cluster_labels):
        go_id = go_list[idx]
        freq = int(term_frequencies.get(go_id, 1))
        clusters[lab].append((go_id, freq))

    print(f"Created {len(clusters)} clusters")

    cluster_data = []
    for cid, items in clusters.items():
        items.sort(key=lambda x: (-x[1], x[0]))  # by freq desc, then GO ID
        rep_go = items[0][0]
        for go_id, freq in items:
            cluster_data.append({
                'TermID': go_id,
                'Name': go_id,  # filled later
                'Representative': f'CLUSTER_{cid}_REP',
                'Representative_Name': rep_go,  # filled later
                'Value': freq / total_freq,
                'Cluster_Size': len(items),
                'Term_Frequency': freq,
                'cluster_id': cid
            })
    return cluster_data

def add_term_names(cluster_data, go_id_to_term, go_dag):
    """Add human-readable term names to cluster data."""
    for item in cluster_data:
        go_id = item['TermID']
        if go_id in go_id_to_term:
            item['Name'] = go_id_to_term[go_id]
        elif go_id in go_dag:
            item['Name'] = go_dag[go_id].name

        rep_go_id = item['Representative_Name']
        if rep_go_id in go_id_to_term:
            item['Representative_Name'] = go_id_to_term[rep_go_id]
        elif rep_go_id in go_dag:
            item['Representative_Name'] = go_dag[rep_go_id].name
    return cluster_data

# -------------------------
# Small runner for one group/ns
# -------------------------

def run_one_namespace(freq_dict, name_map, ontology, go_dag, threshold):
    """Run IC→Resnik→cluster for a single group & ontology; return (df, n_clusters)."""
    out_cols = ['Cluster', 'TermID', 'Name', 'Representative_Name', 'Value']
    if not freq_dict:
        return pd.DataFrame(columns=out_cols), 0

    ic = calculate_custom_information_content(freq_dict.keys(), freq_dict, go_dag, ontology)
    if not ic:
        return pd.DataFrame(columns=out_cols), 0

    sim, go_list = calculate_resnik_similarity_matrix(freq_dict.keys(), ic, go_dag)
    clusters = cluster_similarity_matrix(sim, go_list, freq_dict, threshold)
    for item in clusters:
        item['Name'] = name_map.get(item['TermID'], item['TermID'])
        rep = item['Representative_Name']
        item['Representative_Name'] = name_map.get(rep, rep)

    df = pd.DataFrame(clusters) if clusters else pd.DataFrame(columns=out_cols)
    if df.empty:
        return pd.DataFrame(columns=out_cols), 0

    df['Cluster'] = df['Representative'].str.extract(r'(\d+)').astype(int)
    df = df[out_cols].sort_values(['Cluster', 'Value'], ascending=[True, False])
    return df, df['Cluster'].nunique()

# -------------------------
# Main
# -------------------------

def main():
    parser = argparse.ArgumentParser(description='GOATools clustering with Resnik similarity and custom IC')
    parser.add_argument('--core-annotations', required=True, help='Core annotations file')
    parser.add_argument('--accessory-annotations', required=True, help='Accessory annotations file')
    parser.add_argument('--output-prefix', required=True, help='Output prefix (e.g., Strep_cdhit)')
    parser.add_argument('--method', required=True, help='Method name (for summary only)')
    parser.add_argument('--threshold', type=float, default=0.4, help='Similarity threshold')
    parser.add_argument('--download-go', action='store_true', help='Force download GO ontology')
    args = parser.parse_args()

    if not GOATOOLS_AVAILABLE:
        print("GOATools not available. Please install: pip install goatools")
        sys.exit(1)

    print("GOATools clustering with Resnik similarity and custom Information Content")
    print(f"Similarity threshold: {args.threshold}")

    obo_file = download_go_obo(args.download_go)

    print("Loading GO ontology...")
    go_dag = GODag(obo_file, optional_attrs={'relationship'})
    print(f"Loaded {len(go_dag)} GO terms")

    print("Reading annotations...")
    core_ann = pd.read_csv(args.core_annotations, sep='\t')
    acc_ann  = pd.read_csv(args.accessory_annotations, sep='\t')
    all_ann  = pd.concat([core_ann, acc_ann], ignore_index=True)
    print(f"Processing {len(all_ann)} annotations")

    # Tally per group
    c_bp, c_mf, c_cc, c_bp_names, c_mf_names, c_cc_names = tally_by_namespace(core_ann, go_dag, True)
    a_bp, a_mf, a_cc, a_bp_names, a_mf_names, a_cc_names = tally_by_namespace(acc_ann,  go_dag, True)
    t_bp, t_mf, t_cc, t_bp_names, t_mf_names, t_cc_names = tally_by_namespace(all_ann,  go_dag, True)

    print(f"Core  : {len(c_bp)} BP, {len(c_mf)} MF, {len(c_cc)} CC")
    print(f"Access: {len(a_bp)} BP, {len(a_mf)} MF, {len(a_cc)} CC")
    print(f"All   : {len(t_bp)} BP, {len(t_mf)} MF, {len(t_cc)} CC")

    # Run per ontology, per group
    # Combined (kept as the original filenames for Nextflow compatibility)
    bp_all_df, n_bp_all = run_one_namespace(t_bp, t_bp_names, 'BP', go_dag, args.threshold)
    mf_all_df, n_mf_all = run_one_namespace(t_mf, t_mf_names, 'MF', go_dag, args.threshold)

    bp_all_df.to_csv(f'{args.output_prefix}_goatools_bp_clusters.tsv', sep='\t', index=False)
    mf_all_df.to_csv(f'{args.output_prefix}_goatools_mf_clusters.tsv', sep='\t', index=False)

    # Core
    bp_core_df, n_bp_core = run_one_namespace(c_bp, c_bp_names, 'BP', go_dag, args.threshold)
    mf_core_df, n_mf_core = run_one_namespace(c_mf, c_mf_names, 'MF', go_dag, args.threshold)

    bp_core_df.to_csv(f'{args.output_prefix}_core_goatools_bp_clusters.tsv', sep='\t', index=False)
    mf_core_df.to_csv(f'{args.output_prefix}_core_goatools_mf_clusters.tsv', sep='\t', index=False)

    # Accessory
    bp_acc_df, n_bp_acc = run_one_namespace(a_bp, a_bp_names, 'BP', go_dag, args.threshold)
    mf_acc_df, n_mf_acc = run_one_namespace(a_mf, a_mf_names, 'MF', go_dag, args.threshold)

    bp_acc_df.to_csv(f'{args.output_prefix}_accessory_goatools_bp_clusters.tsv', sep='\t', index=False)
    mf_acc_df.to_csv(f'{args.output_prefix}_accessory_goatools_mf_clusters.tsv', sep='\t', index=False)

    # ---- Summary ----
    summary_lines = [
        "GOATOOLS RESNIK CLUSTERING WITH CUSTOM INFORMATION CONTENT",
        "==========================================================",
        f"Analysis Date: {pd.Timestamp.now()}",
        f"Method: {args.method}",
        f"Similarity Threshold: {args.threshold}",
        f"Algorithm: Resnik similarity with custom IC from term frequencies",
        "",
        "APPROACH:",
        "✓ Uses GOATools Resnik semantic similarity",
        "✓ Custom Information Content from your Streptomyces data",
        "✓ Hierarchical clustering of similarity matrix",
        "✓ Representative = most frequent term in each cluster",
        "",
        "INPUT (counts of unique GO IDs per group):",
        f"  Core:      BP={len(c_bp)}, MF={len(c_mf)}",
        f"  Accessory: BP={len(a_bp)}, MF={len(a_mf)}",
        f"  Combined:  BP={len(t_bp)}, MF={len(t_mf)}",
        "",
        "CLUSTERING RESULTS (number of clusters):",
        f"  Core:      BP={n_bp_core}, MF={n_mf_core}",
        f"  Accessory: BP={n_bp_acc},  MF={n_mf_acc}",
        f"  Combined:  BP={n_bp_all},  MF={n_mf_all}",
        "",
        "OUTPUT FILES (combined for pipeline compatibility):",
        f"  • {args.output_prefix}_goatools_bp_clusters.tsv",
        f"  • {args.output_prefix}_goatools_mf_clusters.tsv",
        "",
        "OUTPUT FILES (separate core/accessory for pangenome interpretation):",
        f"  • {args.output_prefix}_core_goatools_bp_clusters.tsv",
        f"  • {args.output_prefix}_core_goatools_mf_clusters.tsv",
        f"  • {args.output_prefix}_accessory_goatools_bp_clusters.tsv",
        f"  • {args.output_prefix}_accessory_goatools_mf_clusters.tsv",
        "",
        "NOTE:",
        "  The 'Cluster' column is a numeric cluster ID within each file. IDs are not comparable across files."
    ]

    with open(f'{args.output_prefix}_goatools_summary.txt', 'w') as f:
        f.write('\n'.join(summary_lines))

    print("\nGOATools clustering completed!")
    print(f"Core      : BP={n_bp_core} clusters, MF={n_mf_core} clusters")
    print(f"Accessory : BP={n_bp_acc} clusters,  MF={n_mf_acc} clusters")
    print(f"Combined  : BP={n_bp_all} clusters,  MF={n_mf_all} clusters")

if __name__ == '__main__':
    main()

