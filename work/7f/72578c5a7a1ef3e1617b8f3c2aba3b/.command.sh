#!/bin/bash -ue
mkdir -p input_fasta

find FASTA_inputs -follow -type f \( -iname "*.fa" -o -iname "*.faa" -o -iname "*.fasta" -o -iname "*.fas" \) 2>/dev/null | \
  while read file; do
    base=$(basename "$file")
    if [[ "$base" != *"renamed"* && "$base" != *"cdhit"* ]]; then
        ln -sf "$(readlink -f "$file")" input_fasta/
    fi
  done

python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/build_pangenome_tables.py \
    --input-dir input_fasta \
    --output-dir . \
    --name Strep_cdhit \
    --allele-names Strep_cdhit_allele_names.tsv \
    --shared-headers /Users/saba/Documents/GitHub/PangenomeAnalysis/consolidated/shared_headers.tsv
