#!/bin/bash -ue
export PYTHONUNBUFFERED=1

python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/validate_pangenome.py \
    --gene-matrix Strep_cdhit_strain_by_gene.npz \
    --allele-matrix Strep_cdhit_strain_by_allele.npz \
    --input-dir FASTA_inputs \
    --allele-names Strep_cdhit_allele_names.tsv \
    --workers 1 \
    --batch-size 50 \
    --output-dir . 2>&1 | tee Strep_cdhit_validation_summary.txt
