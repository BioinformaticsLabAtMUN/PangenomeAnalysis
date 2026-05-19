README.md
# AlphaFold DB PDB Downloader

Downloads AlphaFold DB PDB files for UniProt accessions using the EBI API.

## Requirements
- bash
- curl
- python3

## Input format
A tab-separated file with UniProt accessions in column 1.
Header row is allowed.

Example:

uniprot_id
A0A0F7VP09
A0A0F7VRE9

## Usage
```bash
download_pdbs.sh <metadata.tsv> <base_out_dir> [batch_size]

Example:
download_pdbs.sh test_metadata.tsv /tmp/pdb_test_run_api 1000

Output
PDB files:
<base_out_dir>/alphafold_pdbs/AF-<ACC>-F1-model_v<latest>.pdb
Mapping TSV:
<base_out_dir>/uniprot_to_pdb_mapping.tsv
Log:
<base_out_dir>/alphafold_pdbs/logs/download.log

 



