#!/bin/bash -ue
echo "Using NCBI annotation with local GTFs"
python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/annotate_ncbi_with_gtf.py \
    --input-summary Strep_cdhit_all_dominant_summary.tsv \
    --core-summary Strep_cdhit_core_dominant_summary.tsv \
    --gtf-directory /Users/saba/Documents/GitHub/PangenomeAnalysis/GTFs \
    --annotate-scope core \
    --output-core-annotations Strep_cdhit_core_genes_annotated.tsv \
    --output-core-merged Strep_cdhit_core_genes_with_annotations.tsv \
    --output-accessory-annotations Strep_cdhit_accessory_genes_annotated.tsv \
    --output-accessory-merged Strep_cdhit_accessory_genes_with_annotations.tsv
