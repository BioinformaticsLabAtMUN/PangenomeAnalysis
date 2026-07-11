#!/bin/bash -ue
python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/extract_dominant_alleles.py \
    --allele-matrix Strep_cdhit_strain_by_allele.npz \
    --allele-labels Strep_cdhit_strain_by_allele.npz.labels.txt \
    --allele-names Strep_cdhit_allele_names.tsv \
    --allele-faa Strep_cdhit_renamed.fasta \
    --core-genes Strep_core_genes.txt \
    --all-dominant-faa Strep_cdhit_all_dominant_alleles.faa \
    --all-dominant-summary Strep_cdhit_all_dominant_summary.tsv \
    --core-dominant-faa Strep_cdhit_core_dominant_alleles.faa \
    --core-dominant-summary Strep_cdhit_core_dominant_summary.tsv \
    --max-alternatives 3
