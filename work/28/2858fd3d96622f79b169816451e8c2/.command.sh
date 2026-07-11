#!/bin/bash -ue
python /Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/consolidate_sequences.py \
    FASTA_inputs \
    consolidated.faa \
    consolidated_swift.faa \
    shared_headers.tsv \
    missing_headers.txt \
    descriptions.tsv \
    organism_map.tsv
