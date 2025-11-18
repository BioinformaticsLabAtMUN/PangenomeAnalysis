#!/usr/bin/env python3
"""
Cascading UniProt Annotator - Optimized for large-scale accessory genes
Handles hundreds of thousands of genes efficiently
"""

import argparse
import pandas as pd
import os
import requests
import time
import json
import hashlib
from pathlib import Path
from io import StringIO
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading
from tqdm import tqdm
import random

class CascadingUniProtAnnotator:
    """Smart cascading annotator that minimizes API calls"""
    
    def __init__(self, cache_dir=None, delay=0.1, max_workers=10):
        self.delay = delay
        self.max_workers = max_workers
        self.session = requests.Session()
        self.session.headers.update({
            'User-Agent': 'PangenomeAnnotator/2.0 (research purpose)'
        })
        
        # Thread-safe progress tracking
        self.lock = threading.Lock()
        self.total_api_calls = 0
        
        # Setup caching
        if cache_dir:
            self.cache_dir = Path(cache_dir)
            self.cache_dir.mkdir(parents=True, exist_ok=True)
        else:
            self.cache_dir = None
    
    def get_cache_key(self, uniprot_id):
        """Generate cache key for a single UniProt ID"""
        if not self.cache_dir:
            return None
        return f"{uniprot_id}.json"
    
    def load_from_cache(self, uniprot_id):
        """Load result from cache"""
        if not self.cache_dir:
            return None
        cache_file = self.cache_dir / self.get_cache_key(uniprot_id)
        if cache_file.exists():
            try:
                with open(cache_file, 'r') as f:
                    return json.load(f)
            except:
                pass
        return None
    
    def save_to_cache(self, uniprot_id, result):
        """Save result to cache"""
        if not self.cache_dir:
            return
        cache_file = self.cache_dir / self.get_cache_key(uniprot_id)
        try:
            with open(cache_file, 'w') as f:
                json.dump(result, f)
        except:
            pass
    
    def fetch_single_annotation(self, uniprot_id):
        """Fetch annotation for a single UniProt ID"""
        
        # Check cache first
        cached_result = self.load_from_cache(uniprot_id)
        if cached_result is not None:
            return cached_result
        
        try:
            with self.lock:
                self.total_api_calls += 1
            
            url = f"https://rest.uniprot.org/uniprotkb/{uniprot_id}"
            params = {
                'format': 'tsv',
                'fields': 'accession,id,protein_name,go_p,go_f,go_c,organism_name,gene_names'
            }
            
            response = self.session.get(url, params=params, timeout=30)
            
            if response.status_code == 200 and response.text.strip():
                df = pd.read_csv(StringIO(response.text), sep='\t')
                if len(df) > 0:
                    row = df.iloc[0]
                    annotation = {
                        'Entry': row.get('Entry', uniprot_id),
                        'Entry Name': row.get('Entry Name', ''),
                        'Protein names': row.get('Protein names', ''),
                        'Gene Ontology (biological process)': row.get('Gene Ontology (biological process)', ''),
                        'Gene Ontology (molecular function)': row.get('Gene Ontology (molecular function)', ''),
                        'Gene Ontology (cellular component)': row.get('Gene Ontology (cellular component)', ''),
                        'Organism': row.get('Organism', ''),
                        'Gene Names': row.get('Gene Names', '')
                    }
                    
                    self.save_to_cache(uniprot_id, annotation)
                    return annotation
            
            # Cache empty result
            self.save_to_cache(uniprot_id, None)
            return None
                    
        except Exception as e:
            self.save_to_cache(uniprot_id, None)
            return None
        
        finally:
            if self.delay > 0:
                time.sleep(self.delay)
    
    def has_go_annotations(self, annotation):
        """Check if annotation has GO terms"""
        if not annotation:
            return False
        
        go_fields = [
            'Gene Ontology (biological process)',
            'Gene Ontology (molecular function)', 
            'Gene Ontology (cellular component)'
        ]
        
        for field in go_fields:
            value = annotation.get(field, '')
            if value and str(value).strip() and str(value).strip() != 'nan':
                return True
        return False
    
    def get_annotation_quality_score(self, annotation):
        """Score annotation quality (higher = better)"""
        if not annotation:
            return 0
        
        score = 1  # Base score
        
        # Bonus for protein name
        protein_name = annotation.get('Protein names', '')
        if protein_name and protein_name.strip():
            score += 2
            
        # Big bonus for GO terms
        if self.has_go_annotations(annotation):
            score += 10
        
        return score
    
    def fetch_batch_parallel(self, protein_ids, desc="Fetching"):
        """Fetch annotations for a batch of protein IDs in parallel with progress bar"""
        
        if not protein_ids:
            return {}
        
        results = {}
        
        # Use tqdm for progress tracking when batch is large
        use_progress = len(protein_ids) > 100
        
        if use_progress:
            pbar = tqdm(total=len(protein_ids), desc=desc, unit="proteins")
        
        with ThreadPoolExecutor(max_workers=self.max_workers) as executor:
            future_to_id = {executor.submit(self.fetch_single_annotation, pid): pid 
                           for pid in protein_ids}
            
            for future in as_completed(future_to_id):
                pid = future_to_id[future]
                try:
                    annotation = future.result()
                    if annotation:
                        results[pid] = annotation
                except Exception as e:
                    pass  # Silently skip errors
                
                if use_progress:
                    pbar.update(1)
        
        if use_progress:
            pbar.close()
        
        return results
    
    def annotate_genes_cascading(self, summary_df, max_alternatives=4, is_accessory=False, 
                                chunk_size=5000, sample_size=None, min_go_rate=0.7):
        """
        Annotate genes using cascading approach to minimize API calls
        
        Args:
            summary_df: DataFrame with gene information
            max_alternatives: Maximum alternative proteins to try
            is_accessory: Whether this is accessory genes (enables optimizations)
            chunk_size: Process genes in chunks of this size
            sample_size: If set, only process a random sample (for testing)
            min_go_rate: Stop early if GO annotation rate exceeds this
        """
        
        # Sample if requested (useful for testing on large datasets)
        if sample_size and len(summary_df) > sample_size:
            print(f"Sampling {sample_size} genes from {len(summary_df)} total")
            summary_df = summary_df.sample(n=sample_size, random_state=42)
        
        print(f"\nCascading UniProt Annotation")
        print(f"Total genes: {len(summary_df)}")
        print(f"Max alternatives: {max_alternatives}")
        if is_accessory and len(summary_df) > 10000:
            print(f"Chunk size: {chunk_size}")
            print(f"Early stopping GO rate: {min_go_rate * 100:.0f}%")
        print("=" * 60)
        
        start_time = time.time()
        
        # Initialize result tracking
        gene_annotations = {}  # gene -> best annotation info
        
        # Process in chunks if accessory and large
        if is_accessory and len(summary_df) > chunk_size:
            return self._annotate_genes_chunked(summary_df, max_alternatives, chunk_size, min_go_rate)
        
        # Original logic for small datasets
        return self._annotate_genes_standard(summary_df, max_alternatives)
    
    def _annotate_genes_chunked(self, summary_df, max_alternatives, chunk_size, min_go_rate):
        """Process large gene sets in chunks with early stopping"""
        
        start_time = time.time()  # Add this line to define start_time
        gene_annotations = {}
        total_chunks = (len(summary_df) + chunk_size - 1) // chunk_size
        
        print(f"\nProcessing {len(summary_df)} genes in {total_chunks} chunks of {chunk_size}")
        print("-" * 60)
        
        for chunk_idx in range(total_chunks):
            start_idx = chunk_idx * chunk_size
            end_idx = min((chunk_idx + 1) * chunk_size, len(summary_df))
            chunk_df = summary_df.iloc[start_idx:end_idx]
            
            print(f"\nChunk {chunk_idx + 1}/{total_chunks}: Genes {start_idx + 1}-{end_idx}")
            
            # Round 1: Dominant proteins only
            dominant_proteins = []
            gene_to_dominant = {}
            
            for idx, row in chunk_df.iterrows():
                gene = row['gene']
                dominant_str = str(row.get('dominant_proteins', ''))
                
                if dominant_str and dominant_str != 'nan':
                    dominant_list = dominant_str.split(',') if ',' in dominant_str else [dominant_str]
                    dominant_protein = dominant_list[0].strip()
                    
                    if dominant_protein:
                        dominant_proteins.append(dominant_protein)
                        gene_to_dominant[gene] = (dominant_protein, idx)
            
            # Fetch dominant proteins
            desc = f"Chunk {chunk_idx + 1}/{total_chunks} - Dominant"
            dominant_annotations = self.fetch_batch_parallel(list(set(dominant_proteins)), desc=desc)
            
            # Apply results
            chunk_annotated = 0
            chunk_with_go = 0
            
            for gene, (protein, idx) in gene_to_dominant.items():
                if protein in dominant_annotations:
                    annotation = dominant_annotations[protein]
                    gene_annotations[gene] = {
                        'annotation': annotation,
                        'protein_id': protein,
                        'source': 'dominant_protein',
                        'idx': idx
                    }
                    chunk_annotated += 1
                    if self.has_go_annotations(annotation):
                        chunk_with_go += 1
            
            # Calculate GO rate for this chunk
            go_rate = chunk_with_go / len(chunk_df) if len(chunk_df) > 0 else 0
            print(f"  Annotated: {chunk_annotated}/{len(chunk_df)} - GO rate: {go_rate * 100:.1f}%")
            
            # Early stopping: if GO rate is good enough, skip alternatives for this chunk
            if go_rate >= min_go_rate:
                print(f"  Skipping alternatives (GO rate {go_rate * 100:.1f}% >= {min_go_rate * 100:.0f}%)")
                continue
            
            # Round 2: Try cluster representatives for genes without GO
            genes_needing_go = []
            for gene in gene_to_dominant:
                if gene in gene_annotations and not self.has_go_annotations(gene_annotations[gene]['annotation']):
                    genes_needing_go.append(gene)
            
            if genes_needing_go and len(genes_needing_go) <= chunk_size * 0.3:  # Only if <30% need GO
                cluster_representatives = []
                gene_to_cluster = {}
                
                for gene in genes_needing_go:
                    idx = gene_annotations[gene]['idx']
                    row = summary_df.loc[idx]
                    cluster_rep = str(row.get('cluster_representative', '')).strip()
                    
                    if cluster_rep and cluster_rep != 'nan':
                        cluster_representatives.append(cluster_rep)
                        gene_to_cluster[gene] = cluster_rep
                
                if cluster_representatives:
                    desc = f"Chunk {chunk_idx + 1}/{total_chunks} - Representatives"
                    cluster_annotations = self.fetch_batch_parallel(list(set(cluster_representatives)), desc=desc)
                    
                    for gene, protein in gene_to_cluster.items():
                        if protein in cluster_annotations:
                            annotation = cluster_annotations[protein]
                            if self.get_annotation_quality_score(annotation) > \
                               self.get_annotation_quality_score(gene_annotations[gene]['annotation']):
                                gene_annotations[gene] = {
                                    'annotation': annotation,
                                    'protein_id': protein,
                                    'source': 'cluster_representative',
                                    'idx': gene_annotations[gene]['idx']
                                }
        
        elapsed = time.time() - start_time
        
        # Calculate final statistics
        total_annotated = len(gene_annotations)
        go_annotated = sum(1 for ga in gene_annotations.values() 
                          if self.has_go_annotations(ga['annotation']))
        
        print("=" * 60)
        print(f"Chunked Annotation Completed")
        print(f"Total time: {elapsed:.1f}s ({elapsed/60:.1f} minutes)")
        print(f"Total API calls: {self.total_api_calls}")
        print(f"Final results:")
        print(f"  - Annotated genes: {total_annotated}/{len(summary_df)} ({total_annotated/len(summary_df)*100:.1f}%)")
        print(f"  - With GO terms: {go_annotated}/{len(summary_df)} ({go_annotated/len(summary_df)*100:.1f}%)")
        if elapsed > 0:
            print(f"Speed: {self.total_api_calls/elapsed:.1f} API calls/second")
        
        return gene_annotations
    
    def _annotate_genes_standard(self, summary_df, max_alternatives):
        """Original annotation logic for smaller gene sets"""
        
        start_time = time.time()  # Add this line to define start_time
        gene_annotations = {}
        
        # Round 1: Try dominant proteins for all genes
        print("\nRound 1: Dominant Proteins")
        print("-" * 30)
        
        dominant_proteins = []
        gene_to_dominant = {}
        
        for idx, row in summary_df.iterrows():
            gene = row['gene']
            dominant_str = str(row.get('dominant_proteins', ''))
            
            if dominant_str and dominant_str != 'nan':
                dominant_list = dominant_str.split(',') if ',' in dominant_str else [dominant_str]
                dominant_protein = dominant_list[0].strip()
                
                if dominant_protein:
                    dominant_proteins.append(dominant_protein)
                    gene_to_dominant[gene] = (dominant_protein, idx)
        
        # Fetch dominant protein annotations
        print(f"  Fetching {len(set(dominant_proteins))} unique proteins...")
        dominant_annotations = self.fetch_batch_parallel(list(set(dominant_proteins)))
        
        # Apply dominant protein results
        annotated_genes = set()
        
        for gene, (protein, idx) in gene_to_dominant.items():
            if protein in dominant_annotations:
                annotation = dominant_annotations[protein]
                gene_annotations[gene] = {
                    'annotation': annotation,
                    'protein_id': protein,
                    'source': 'dominant_protein',
                    'idx': idx
                }
                annotated_genes.add(gene)
        
        success_rate = len(annotated_genes) / len(summary_df) * 100
        go_count = sum(1 for g in annotated_genes if self.has_go_annotations(gene_annotations[g]['annotation']))
        go_rate = go_count / len(annotated_genes) * 100 if annotated_genes else 0
        
        print(f"  Results: {len(annotated_genes)}/{len(summary_df)} ({success_rate:.1f}%)")
        print(f"  With GO: {go_count} ({go_rate:.1f}%)")
        
        # Round 2: Try cluster representatives for genes without GO annotations
        genes_needing_go = []
        for gene in gene_annotations:
            if not self.has_go_annotations(gene_annotations[gene]['annotation']):
                genes_needing_go.append(gene)
        
        if genes_needing_go:
            print(f"\nRound 2: Cluster Representatives (for {len(genes_needing_go)} genes without GO)")
            print("-" * 30)
            
            cluster_representatives = []
            gene_to_cluster = {}
            
            for gene in genes_needing_go:
                idx = gene_annotations[gene]['idx']
                row = summary_df.iloc[idx]
                cluster_rep = str(row.get('cluster_representative', '')).strip()
                
                if cluster_rep and cluster_rep != 'nan':
                    cluster_representatives.append(cluster_rep)
                    gene_to_cluster[gene] = cluster_rep
            
            if cluster_representatives:
                print(f"  Fetching {len(set(cluster_representatives))} unique proteins...")
                cluster_annotations = self.fetch_batch_parallel(list(set(cluster_representatives)))
                
                # Update genes with better cluster representative annotations
                improved_genes = 0
                for gene, protein in gene_to_cluster.items():
                    if protein in cluster_annotations:
                        annotation = cluster_annotations[protein]
                        current_score = self.get_annotation_quality_score(gene_annotations[gene]['annotation'])
                        new_score = self.get_annotation_quality_score(annotation)
                        
                        if new_score > current_score:
                            gene_annotations[gene] = {
                                'annotation': annotation,
                                'protein_id': protein,
                                'source': 'cluster_representative',
                                'idx': gene_annotations[gene]['idx']
                            }
                            improved_genes += 1
                
                print(f"  Improved annotations: {improved_genes} genes")
        
        # Round 3+: Try alternative proteins for genes still without GO annotations
        for alt_round in range(max_alternatives):
            genes_still_needing_go = []
            for gene in gene_annotations:
                if not self.has_go_annotations(gene_annotations[gene]['annotation']):
                    genes_still_needing_go.append(gene)
            
            if not genes_still_needing_go:
                break
                
            print(f"\nRound {3 + alt_round}: Alternative Protein #{alt_round + 1} (for {len(genes_still_needing_go)} genes)")
            print("-" * 30)
            
            alt_proteins = []
            gene_to_alt = {}
            
            for gene in genes_still_needing_go:
                idx = gene_annotations[gene]['idx']
                row = summary_df.iloc[idx]
                alt_str = str(row.get('alternative_proteins', ''))
                
                if alt_str and alt_str != 'nan':
                    alt_list = alt_str.split(',')
                    if len(alt_list) > alt_round:
                        alt_protein = alt_list[alt_round].strip()
                        if alt_protein:
                            alt_proteins.append(alt_protein)
                            gene_to_alt[gene] = alt_protein
            
            if alt_proteins:
                print(f"  Fetching {len(set(alt_proteins))} unique proteins...")
                alt_annotations = self.fetch_batch_parallel(list(set(alt_proteins)))
                
                # Update genes with better alternative annotations
                improved_genes = 0
                for gene, protein in gene_to_alt.items():
                    if protein in alt_annotations:
                        annotation = alt_annotations[protein]
                        current_score = self.get_annotation_quality_score(gene_annotations[gene]['annotation'])
                        new_score = self.get_annotation_quality_score(annotation)
                        
                        if new_score > current_score:
                            gene_annotations[gene] = {
                                'annotation': annotation,
                                'protein_id': protein,
                                'source': f'alternative_protein_{alt_round + 1}',
                                'idx': gene_annotations[gene]['idx']
                            }
                            improved_genes += 1
                
                print(f"  Improved annotations: {improved_genes} genes")
            else:
                print("  No alternative proteins available")
        
        elapsed = time.time() - start_time
        
        # Calculate final statistics
        total_annotated = len(gene_annotations)
        go_annotated = sum(1 for ga in gene_annotations.values() 
                          if self.has_go_annotations(ga['annotation']))
        
        print("=" * 60)
        print(f"Cascading Annotation Completed")
        print(f"Total time: {elapsed:.1f}s ({elapsed/60:.1f} minutes)")
        print(f"Total API calls: {self.total_api_calls}")
        print(f"Final results:")
        print(f"  - Annotated genes: {total_annotated}/{len(summary_df)} ({total_annotated/len(summary_df)*100:.1f}%)")
        print(f"  - With GO terms: {go_annotated}/{len(summary_df)} ({go_annotated/len(summary_df)*100:.1f}%)")
        if elapsed > 0:
            print(f"Speed: {self.total_api_calls/elapsed:.1f} API calls/second")
        
        return gene_annotations

def create_annotations_from_cascading(gene_annotations, summary_df, output_file):
    """Create final annotation file from cascading results"""
    
    result_df = summary_df.copy()
    
    # Add annotation columns
    annotation_columns = {
        'protein_id': '',
        'locus_tag': '',
        'gene_id': '',
        'product': '',
        'go_ids': '',
        'go_process': '',
        'go_function': '',
        'go_component': '',
        'annotation_source': ''
    }
    
    for col in annotation_columns:
        result_df[col] = annotation_columns[col]
    
    # Fill in results
    for gene, gene_data in gene_annotations.items():
        idx = gene_data['idx']
        annotation = gene_data['annotation']
        
        result_df.at[idx, 'protein_id'] = gene_data['protein_id']
        result_df.at[idx, 'annotation_source'] = gene_data['source']
        result_df.at[idx, 'product'] = annotation.get('Protein names', '')
        
        # GO terms
        go_process = annotation.get('Gene Ontology (biological process)', '')
        go_function = annotation.get('Gene Ontology (molecular function)', '')
        go_component = annotation.get('Gene Ontology (cellular component)', '')
        
        result_df.at[idx, 'go_process'] = go_process
        result_df.at[idx, 'go_function'] = go_function
        result_df.at[idx, 'go_component'] = go_component
        
        # Extract GO IDs
        go_ids = []
        for go_field in [go_process, go_function, go_component]:
            if pd.notna(go_field) and go_field:
                import re
                go_matches = re.findall(r'GO:\d+', str(go_field))
                go_ids.extend(go_matches)
        
        result_df.at[idx, 'go_ids'] = ';'.join(list(dict.fromkeys(go_ids))) if go_ids else ''
    
    # Save results
    result_df.to_csv(output_file, sep='\t', index=False)
    print(f"Annotations saved: {output_file}")
    
    return result_df

def separate_core_accessory_genes(all_summary_file, core_summary_file):
    """Separate genes into core and accessory"""
    all_df = pd.read_csv(all_summary_file, sep='\t')
    core_df = pd.read_csv(core_summary_file, sep='\t')
    core_genes = set(core_df['gene'].tolist())
    accessory_df = all_df[~all_df['gene'].isin(core_genes)]
    
    print(f"Gene separation:")
    print(f"  - Total genes: {len(all_df)}")
    print(f"  - Core genes: {len(core_df)}")
    print(f"  - Accessory genes: {len(accessory_df)}")
    
    return core_df, accessory_df

def main():
    parser = argparse.ArgumentParser(description='Cascading UniProt annotation - optimized for large-scale')
    parser.add_argument('--input-summary', required=True, help='Input all genes summary TSV file')
    parser.add_argument('--core-summary', required=True, help='Input core genes summary TSV file')
    parser.add_argument('--annotate-scope', choices=['core', 'all'], default='all', help='Annotation scope')
    parser.add_argument('--output-core-annotations', required=True, help='Output core annotations file')
    parser.add_argument('--output-core-merged', required=True, help='Output core merged file')
    parser.add_argument('--output-accessory-annotations', required=True, help='Output accessory annotations file')
    parser.add_argument('--output-accessory-merged', required=True, help='Output accessory merged file')
    parser.add_argument('--max-alternatives', type=int, default=4, help='Maximum alternative proteins to try')
    parser.add_argument('--delay', type=float, default=0.1, help='Delay between requests (seconds)')
    parser.add_argument('--workers', type=int, default=10, help='Number of parallel workers')
    parser.add_argument('--chunk-size', type=int, default=5000, help='Chunk size for accessory genes')
    parser.add_argument('--sample-size', type=int, default=None, help='Sample size for testing (optional)')
    parser.add_argument('--min-go-rate', type=float, default=0.7, help='Minimum GO rate for early stopping')
    
    args = parser.parse_args()
    
    print("Cascading UniProt Annotation - Optimized for Large-Scale")
    print("=" * 60)
    print(f"Annotation scope: {args.annotate_scope}")
    print(f"Max alternatives: {args.max_alternatives}")
    print(f"Parallel workers: {args.workers}")
    print(f"Chunk size: {args.chunk_size}")
    if args.sample_size:
        print(f"Sample size: {args.sample_size}")
    print(f"Min GO rate for early stopping: {args.min_go_rate * 100:.0f}%")
    
    # Setup cache directory
    cache_dir = os.path.join(os.path.dirname(args.output_core_annotations), '.uniprot_cache_cascading')
    
    # Separate core and accessory genes
    core_df, accessory_df = separate_core_accessory_genes(args.input_summary, args.core_summary)
    
    # Initialize annotator
    annotator = CascadingUniProtAnnotator(cache_dir=cache_dir, delay=args.delay, max_workers=args.workers)
    
    # Process core genes
    print("\nProcessing Core Genes")
    print("-" * 40)
    
    core_gene_annotations = annotator.annotate_genes_cascading(
        core_df, 
        args.max_alternatives,
        is_accessory=False
    )
    create_annotations_from_cascading(core_gene_annotations, core_df, args.output_core_annotations)
    
    # Create merged output
    core_df_annotated = pd.read_csv(args.output_core_annotations, sep='\t')
    core_df_annotated.to_csv(args.output_core_merged, sep='\t', index=False)
    print(f"Core merged output: {args.output_core_merged}")
    
    # Process accessory genes (if scope is "all")
    if args.annotate_scope == "all" and len(accessory_df) > 0:
        print("\nProcessing Accessory Genes")
        print("-" * 40)
        
        accessory_gene_annotations = annotator.annotate_genes_cascading(
            accessory_df, 
            args.max_alternatives,
            is_accessory=True,
            chunk_size=args.chunk_size,
            sample_size=args.sample_size,
            min_go_rate=args.min_go_rate
        )
        create_annotations_from_cascading(accessory_gene_annotations, accessory_df, args.output_accessory_annotations)
        
        # Create merged output
        accessory_df_annotated = pd.read_csv(args.output_accessory_annotations, sep='\t')
        accessory_df_annotated.to_csv(args.output_accessory_merged, sep='\t', index=False)
        print(f"Accessory merged output: {args.output_accessory_merged}")
    else:
        print("\nSkipping accessory genes")
        pd.DataFrame().to_csv(args.output_accessory_annotations, sep='\t', index=False)
        pd.DataFrame().to_csv(args.output_accessory_merged, sep='\t', index=False)
    
    print("\nCascading Annotation Completed")
    print("Optimized approach for large-scale gene sets!")

if __name__ == "__main__":
    main()