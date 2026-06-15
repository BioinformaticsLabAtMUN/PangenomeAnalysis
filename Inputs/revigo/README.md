# Revigo Input Files

This directory contains sample input files for the Revigo-based functional analysis step of the pipeline (`analyzeRevigoFunctionalCore` process).

## What Revigo is used for

Revigo (https://revigo.irb.hr/) takes a list of GO (Gene Ontology) terms and summarizes/clusters them by removing redundant terms, grouping similar ones, and visualizing them based on semantic similarity. The pipeline uses Revigo's output table to interpret the functional categories of the core genome.

## How to generate the input files

1. **Get a list of GO terms with their frequencies.**
   From your core genome annotations (e.g. `*_core_genes_annotated.tsv`), extract the GO terms associated with Biological Process (BP) and Molecular Function (MF) categories, and count how many times each GO term appears across the core genes.

   Your input list should look like:
   ```
   GO:0006412	15
   GO:0008152	9
   GO:0003735	7
   ...
   ```
   (GO term ID, tab, frequency/count)

2. **Submit to Revigo.**
   - Go to https://revigo.irb.hr/
   - Paste or upload your list of GO terms (with frequencies as the numeric value)
   - Run Revigo separately for **Biological Process** and **Molecular Function** terms
   - Once finished, download/export the results table for each (CSV or TSV)

3. **Save the output tables.**
   These exported tables are the files the pipeline needs as input for the Revigo functional core analysis step.

## Files needed by the pipeline

For each clustering method (CD-HIT, SwiftOrtho, Foldseek), the pipeline expects two Revigo tables:

| Parameter (in `nextflow.config`) | Description |
|---|---|
| `revigo_cdhit_bp_file`       | Revigo output table for CD-HIT core genome — Biological Process GO terms |
| `revigo_cdhit_mf_file`       | Revigo output table for CD-HIT core genome — Molecular Function GO terms |
| `revigo_swiftortho_bp_file`  | Revigo output table for SwiftOrtho core genome — Biological Process GO terms |
| `revigo_swiftortho_mf_file`  | Revigo output table for SwiftOrtho core genome — Molecular Function GO terms |
| `revigo_foldseek_bp_file`    | Revigo output table for Foldseek core genome — Biological Process GO terms |
| `revigo_foldseek_mf_file`    | Revigo output table for Foldseek core genome — Molecular Function GO terms |

If these files are not present, the pipeline will log a warning and skip the Revigo functional analysis step for that clustering method. Everything else will still run normally.

## Sample files

This directory contains example Revigo output tables, so you can see the expected format and test the analysis step.


## Setting the file paths

Once you have your own Revigo output tables, point to them in `nextflow.config`:

```groovy
params {
    revigo_cdhit_bp_file = "${params.baseDir}/Inputs/revigo/my_cdhit_bp_output.tsv"
    revigo_cdhit_mf_file = "${params.baseDir}/Inputs/revigo/my_cdhit_mf_output.tsv"
    // ... and similarly for swiftortho / foldseek if used
}
```