# Pangenome Analysis Pipeline

A Nextflow pipeline for pangenome analysis of bacterial proteomes using CD-HIT, Foldseek, or SwiftOrtho clustering, with UniProt functional annotation and core genome analysis.

---

## What this pipeline does

1. Merges all your input protein FASTA files into one consolidated file
2. Clusters similar proteins across all strains using your chosen method
3. Builds presence/absence matrices (which genes are in which strains)
4. Identifies core genes — those present in all or most strains
5. Annotates core and accessory genes using UniProt
6. Produces plots, tables, and summary statistics

---

## Requirements

You need two things installed before running anything:

**Nextflow**
```bash
curl -s https://get.nextflow.io | bash
sudo mv nextflow /usr/local/bin/
nextflow -version   # should print a version number
```

**Conda** (Miniconda is fine)
- Download from: https://docs.conda.io/en/latest/miniconda.html
- All Python packages are installed automatically by the pipeline on first run. You do not need to install anything else manually.

---

## Pipeline directory structure

Your pipeline folder should look like this before running:

```
pipeline/
├── main.nf
├── modules.nf
├── nextflow.config        ← this is where you change settings
├── envs/
│   ├── pangenome_env.yml
│   ├── goatools_env.yml
│   └── swiftortho_env.yml
├── scripts/               ← all Python scripts live here
├── fasta_foldseek/        ← your input FASTA files go here
└── my_proteome_metadata.tsv
```

---

## Input files you need to provide

### 1. Protein FASTA files (required)

One `.fasta`, `.faa`, or `.fa` file per strain, all in the same directory. The pipeline reads all files in that directory automatically.

```
fasta_foldseek/
├── UP000000377.fasta
├── UP000000428.fasta
└── ...
```

### 2. Proteome metadata file (required)

A tab-separated file with one row per strain. Example:

```
Proteome_ID    Organism                     Strain
UP000000377    Streptomyces coelicolor      A3(2)
UP000000428    Streptomyces avermitilis     MA-4680
```

---

## How to change settings

**All settings are in `nextflow.config`.** Open it in any text editor. The most important ones are near the top under `params`:

```groovy
params {
    name_prefix            = 'Strep'          // prefix on all output files
    raw_input_directory    = "${params.baseDir}/fasta_foldseek"
    proteome_metadata_file = "${params.baseDir}/my_proteome_metadata.tsv"
    clustering_method      = "cdhit"          // cdhit | foldseek | swiftortho | all
    database               = "uniprot"        // uniprot | ncbi
    threads                = 12               // CPU threads to use
    cdhit_identity         = 0.65             // CD-HIT identity threshold (0–1)
    cdhit_coverage         = 0.75             // CD-HIT coverage threshold (0–1)
}
```

You can also override any setting on the command line without editing the file — just add `--parameter_name value`:

```bash
nextflow run main.nf --threads 24 --cdhit_identity 0.70
```

Command line values always take priority over `nextflow.config`.

---

## How conda environments work

You do **not** need to create or activate any conda environment yourself. When you run the pipeline for the first time, Nextflow reads the `.yml` files in `envs/` and builds the environments automatically. This happens once and takes about 5–10 minutes. Every subsequent run reuses the cached environments instantly.

The environments are saved to `.conda_cache/` inside your working directory. If you delete that folder, Nextflow rebuilds them on the next run.

Three environments are used:
- `pangenome_env` — used by almost all steps (Python analysis, CD-HIT)
- `goatools_env` — used for GO term clustering
- `swiftortho_env` — used only if you run SwiftOrtho

---

## Running the pipeline

Navigate into the pipeline directory first:

```bash
cd /path/to/pipeline
```

**Basic run (CD-HIT, UniProt annotation):**
```bash
nextflow run main.nf
```

**Specify a clustering method:**
```bash
nextflow run main.nf --clustering_method cdhit
```

**Custom input paths:**
```bash
nextflow run main.nf \
  --raw_input_directory /path/to/your/fastas \
  --proteome_metadata_file /path/to/metadata.tsv
```

**Resume an interrupted run:**
```bash
nextflow run main.nf -resume
```

`-resume` tells Nextflow to skip any steps that already finished successfully. You can always resume safely — it never redoes completed work.

---
## Output files

Everything is written to `output/uniprot_output/` in your working directory:

```
output/
└── uniprot_output/
    ├── annotations/cdhit/
    │   ├── Strep_cdhit_core_genes_annotated.tsv       ← core gene annotations with GO terms
    │   └── Strep_cdhit_accessory_genes_annotated.tsv
    ├── core_genome/cdhit/
    │   ├── Strep_core_genes.txt                       ← list of core gene IDs
    │   └── Strep_beta_binomial_fit.png
    ├── gene_structure_analysis/cdhit/
    │   ├── Strep_cdhit_pangenome_composition_simplified.png
    │   └── Strep_cdhit_missing_core_genes_summary_table.tsv
    ├── heaps_analysis/cdhit/
    │   └── Strep_heaps_law_plot.png                   ← pangenome openness plot
    └── pangenome_tables/cdhit/
        └── Strep_cdhit_strain_by_gene.npz             ← presence/absence matrix
```

The most useful files:
- `core_genes_annotated.tsv` — your annotated core genome
- `pangenome_composition_simplified.png` — quick visual summary
- `heaps_law_plot.png` — shows whether the pangenome is open or closed

---

## Troubleshooting

**"Input directory does not exist"**
Check that `raw_input_directory` in `nextflow.config` points to your FASTA files, or override it:
```bash
nextflow run main.nf --raw_input_directory /absolute/path/to/fastas
```

**"conda: command not found"**
Conda is not on your PATH. Close and reopen your terminal, or run `source ~/.bashrc`.

**UniProt annotation is slow or failing**
You may be hitting the UniProt API rate limit. Reduce parallel requests:
```bash
nextflow run main.nf --uniprot_max_workers 4 --uniprot_batch_size 200
```

**Run crashed partway through**
Re-run with `-resume` — Nextflow picks up where it left off:
```bash
nextflow run main.nf -resume
```

**Out of memory**
Edit the `process` block in `nextflow.config` and increase memory for the failing step. The step name is shown in the error message.

---

## Citation

If you use this pipeline, please cite:

> Sadeghi Najabadi, S. (2026). *Comparison of sequence-based vs structure-based pangenome analyses of Streptomyces*. MSc thesis, Memorial University of Newfoundland.
