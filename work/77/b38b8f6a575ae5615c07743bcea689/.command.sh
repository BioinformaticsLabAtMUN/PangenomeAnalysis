#!/bin/bash -ue
SCRIPT_PATH="/Users/saba/Documents/GitHub/PangenomeAnalysis/scripts/unified_rename_sequences_simple_representative.py"

echo -e "UniProt_ID\tPDB_file\tStatus" > dummy_mapping.tsv

if [[ "cdhit" != "foldseek" ]]; then
    MAPPING_FILE=dummy_mapping.tsv
else
    MAPPING_FILE=dummy_mapping.tsv
fi

python $SCRIPT_PATH \
    cdhit_output.faa.clstr \
    consolidated.faa \
    Strep_cdhit_renamed.fasta \
    Strep_cdhit_allele_names.tsv \
    Strep \
    cdhit \
    shared_headers.tsv \
    $MAPPING_FILE
