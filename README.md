# Pangenome Analysis Pipeline

A comprehensive Nextflow pipeline for pangenomic analysis using multiple clustering methods (CD-HIT, SwiftOrtho, and Foldseek) with automated functional annotation via UniProt API.

## Table of Contents
- [Overview](#overview)
- [Features](#features)
- [Requirements](#requirements)
- [Installation](#installation)
- [Input Data](#input-data)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Output Structure](#output-structure)
- [Clustering Methods](#clustering-methods)
- [Troubleshooting](#troubleshooting)

---

## Overview

This pipeline performs comprehensive pangenome analysis on protein sequences from multiple strains/genomes. It supports three clustering methods and integrates automated functional annotation through the UniProt API.

**Key capabilities:**
- Multi-method protein clustering (CD-HIT, SwiftOrtho, Foldseek)
- Automated consolidation of input sequences
- Core and accessory genome identification
- Functional annotation via UniProt
- GO term enrichment analysis with Revigo
- Heaps' Law analysis for pangenome dynamics
- Comprehensive visualization and validation

---

## Features

### Clustering Methods
1. **CD-HIT**: Sequence similarity-based clustering (fast, sequence-only)
2. **SwiftOrtho**: Orthology-based clustering (balanced approach)
3. **Foldseek**: Structure-based clustering (requires PDB files)

### Analysis Modules
- Core genome identification and validation
- Accessory genome characterization
- Pangenome size dynamics (Heaps' Law)
- Functional annotation using UniProt API
- GO term enrichment analysis
- Gene structure analysis
- Allele diversity assessment

---

## Requirements

### Software Dependencies
- **Nextflow** (≥22.0)
- **Conda/Mamba** (for environment management)
- **CD-HIT** (for sequence clustering)
- **SwiftOrtho** (for orthology-based clustering)
- **Foldseek** (for structure-based clustering)
- **Python** (≥3.8)

### Python Packages
Create a `requirements.txt` file with:
```
biopython>=1.79
numpy>=1.21.0
pandas>=1.3.0
scipy>=1.7.0
matplotlib>=3.4.0
seaborn>=0.11.0
requests>=2.26.0
tqdm>=4.62.0
fastcluster>=1.2.0
```

### Conda Environments
The pipeline requires two conda environments:

**pangenome_env** (for most analyses):
```bash
conda create -n pangenome_env python=3.9
conda activate pangenome_env
pip install -r requirements.txt
```

**swiftortho_env** (for SwiftOrtho only):
```bash
conda create -n swiftortho_env python=2.7
conda activate swiftortho_env
conda install -c bioconda blast
# Install SwiftOrtho from: https://github.com/Rinoahu/SwiftOrtho
```

---

## Installation

### 1. Clone the Repository
```bash
git clone <your-repository-url>
cd pangenome-pipeline
```

### 2. Install Dependencies
```bash
# Install CD-HIT
conda install -c bioconda cd-hit

# Install Foldseek
wget https://mmseqs.com/foldseek/foldseek-linux-avx2.tar.gz
tar xvfz foldseek-linux-avx2.tar.gz
sudo cp foldseek/bin/foldseek /usr/local/bin/

# Install SwiftOrtho
git clone https://github.com/Rinoahu/SwiftOrtho.git tools/SwiftOrtho
```

### 3. Create Conda Environments
```bash
# Create pangenome environment
conda create -n pangenome_env python=3.9 -y
conda activate pangenome_env
pip install -r requirements.txt

# Create SwiftOrtho environment
conda create -n swiftortho_env python=2.7 -y
conda activate swiftortho_env
conda install -c bioconda blast -y
```

### 4. Set Up Directory Structure
```bash
mkdir -p {input,output,gtf,cache,revigo,tools}
mkdir -p scriptse envs bin
```

---

## Input Data

### Required Input Files

#### 1. Protein FASTA Files
Place your protein sequence files in a directory:
```
input_directory/
├── strain_001.faa
├── strain_002.faa
├── strain_003.faa
└── ...
```

**Example FASTA format:**
```
>UniProtID_001
MKLLVLGLPGAGKGTQAQFITEKFPVCHQTYIEGHLSATGSTFRMKVLVVGESGNGGKE
ITRIGKVVFNEPTAAALATKVAKLLSKELLRNIDRNVVVILDEVDMALAEILAKNPNSG
>UniProtID_002
MSTDFVIADRGEIDTPDEHVKKVLVEEGVDILVFGYPAFYYTPQVIELFTRNGFIQLDE
```

#### 2. Proteome Metadata File
Create a TSV file with strain/genome information:

**File:** `my_proteome_metadata.tsv`
```tsv
Accession	Organism	Strain	Assembly
GCF_000001.1	Streptomyces coelicolor	A3(2)	GCA_000001.1
GCF_000002.1	Streptomyces griseus	NBRC 13350	GCA_000002.1
GCF_000003.1	Streptomyces venezuelae	ATCC 10712	GCA_000003.1
```

#### 3. PDB Files (For Foldseek Only)
If using structure-based clustering:
```
pdb_directory/
├── UniProtID_001.pdb
├── UniProtID_002.pdb
└── ...
```

---

## Configuration

### Basic Configuration
Edit `pangenome_uniprot.nf` or create a `nextflow.config` file:

```groovy
params {
    // Basic settings
    name_prefix = 'MyProject'
    baseDir = "$PWD"
    
    // Input directories
    raw_input_directory = "${baseDir}/input_fasta"
    proteome_metadata_file = "${baseDir}/my_proteome_metadata.tsv"
    
    // Clustering method
    clustering_method = "cdhit"  // Options: "cdhit", "swiftortho", "foldseek", "all"
    
    // Database for annotation
    database = "uniprot"  // Options: "uniprot", "ncbi"
    
    // CD-HIT parameters
    cdhit_identity = 0.65     // 65% sequence identity
    cdhit_coverage = 0.75     // 75% coverage
    threads = 12              // Number of CPU threads
    
    // Output directory
    output_dir = "${baseDir}/output"
}
```

### Advanced Configuration

#### CD-HIT Parameters
```groovy
params {
    cdhit_identity = 0.90      // Sequence identity threshold (0.4-1.0)
    cdhit_coverage = 0.80      // Coverage threshold (0.0-1.0)
    cdhit_mode = "accurate"    // "fast" or "accurate"
    threads = 16               // CPU threads
}
```

#### Foldseek Parameters
```groovy
params {
    foldseek_path = '/usr/local/bin/foldseek'
    pdb_directory = "${baseDir}/pdb_files"
    foldseek_sensitivity = 9.5
    foldseek_min_seq_id = 0.3
    foldseek_threads = 10
}
```

#### UniProt API Settings
```groovy
params {
    uniprot_batch_size = 500      // Batch size for API requests (100-500)
    uniprot_max_workers = 8       // Concurrent API workers (4-12)
}
```

---

## Running the Pipeline

### Quick Start

#### 1. Basic Run with CD-HIT
```bash
nextflow run pangenome_uniprot.nf \
    --name_prefix MyProject \
    --raw_input_directory ./input_fasta \
    --proteome_metadata_file ./metadata.tsv \
    --clustering_method cdhit \
    --threads 12
```

#### 2. Run All Methods
```bash
nextflow run pangenome_uniprot.nf \
    --name_prefix MyProject \
    --raw_input_directory ./input_fasta \
    --proteome_metadata_file ./metadata.tsv \
    --clustering_method all \
    --threads 12
```

#### 3. Structure-Based Clustering with Foldseek
```bash
nextflow run pangenome_uniprot.nf \
    --name_prefix MyProject \
    --raw_input_directory ./pdb_directory \
    --proteome_metadata_file ./metadata.tsv \
    --clustering_method foldseek \
    --pdb_directory ./pdb_files \
    --foldseek_threads 10
```

### Step-by-Step Execution

#### Step 1: Prepare Input Data
```bash
# Create directory structure
mkdir -p input_fasta metadata

# Copy your FASTA files
cp /path/to/your/fastas/*.faa input_fasta/

# Create metadata file
cat > metadata.tsv << EOF
Accession	Organism	Strain	Assembly
strain001	Streptomyces_sp	001	GCA_001
strain002	Streptomyces_sp	002	GCA_002
EOF
```

#### Step 2: Test with Small Dataset
```bash
# Run with 3 genomes first
nextflow run pangenome_uniprot.nf \
    --name_prefix Test \
    --raw_input_directory ./input_fasta \
    --proteome_metadata_file ./metadata.tsv \
    --clustering_method cdhit \
    --threads 4 \
    -resume
```

#### Step 3: Full Production Run
```bash
# Run complete pipeline
nextflow run pangenome_uniprot.nf \
    --name_prefix Production \
    --raw_input_directory ./input_fasta \
    --proteome_metadata_file ./metadata.tsv \
    --clustering_method all \
    --threads 16 \
    --cdhit_identity 0.70 \
    --cdhit_coverage 0.80 \
    -with-report execution_report.html \
    -with-timeline timeline.html \
    -with-dag flowchart.png \
    -resume
```

### Resume Failed Runs
```bash
# Nextflow automatically resumes from last successful step
nextflow run pangenome_uniprot.nf -resume
```

---

## Output Structure

### Directory Layout
```
output/
├── uniprot_output/
│   ├── cdhit/
│   │   ├── cdhit_output.faa              # Clustered sequences
│   │   └── cdhit_output.faa.clstr        # Cluster file
│   ├── swiftortho/
│   │   └── clusters/
│   ├── foldseek/
│   │   ├── database/
│   │   └── clusters/
│   ├── renamed_sequences/
│   │   ├── cdhit/
│   │   │   ├── MyProject_cdhit_renamed.fasta
│   │   │   └── MyProject_cdhit_allele_names.tsv
│   │   ├── swiftortho/
│   │   └── foldseek/
│   ├── pangenome_tables/
│   │   ├── cdhit/
│   │   │   ├── MyProject_cdhit_strain_by_allele.npz
│   │   │   ├── MyProject_cdhit_strain_by_gene.npz
│   │   │   ├── MyProject_cdhit_strain_by_allele.npz.labels.txt
│   │   │   └── MyProject_cdhit_strain_by_gene.npz.labels.txt
│   │   ├── swiftortho/
│   │   └── foldseek/
│   ├── validation/
│   │   └── cdhit/
│   │       └── MyProject_cdhit_validation_summary.txt
│   ├── pangenome_visualizations/
│   │   └── cdhit/
│   │       └── MyProject_cdhit_gene_heatmap.png
│   ├── heaps_analysis/
│   │   └── cdhit/
│   │       ├── MyProject_heaps_law_plot.png
│   │       └── MyProject_heaps_law_results.csv
│   ├── core_genome_analysis/
│   │   └── cdhit/
│   │       ├── core_genome_genes.txt
│   │       ├── core_genome_summary.json
│   │       └── missing_core_genes_top50_strains.png
│   ├── dominant_alleles/
│   │   └── cdhit/
│   │       ├── MyProject_dominant_alleles_with_core.faa
│   │       └── MyProject_dominant_summary_with_core.tsv
│   ├── annotations/
│   │   └── cdhit/
│   │       ├── annotations_merged.tsv
│   │       └── annotations_summary.txt
│   ├── gene_structure/
│   │   └── cdhit/
│   │       ├── gene_structure_analysis.txt
│   │       └── gene_structure_plots.png
│   └── functional_core/
│       └── cdhit/
│           ├── functional_core_analysis.txt
│           └── go_enrichment_results.tsv
└── cache/                                 # Cached intermediate files
```

### Key Output Files

#### 1. Cluster Files
**File:** `cdhit_output.faa.clstr`
```
>Cluster 0
0	5366aa, >strain001|protein001... *
1	5352aa, >strain002|protein045... at 95.2%
2	5344aa, >strain003|protein123... at 94.8%
>Cluster 1
0	2134aa, >strain001|protein002... *
1	2130aa, >strain004|protein089... at 96.1%
```

#### 2. Renamed Sequences
**File:** `MyProject_cdhit_renamed.fasta`
```
>MyProject_gene_0001_allele_001
MKLLVLGLPGAGKGTQAQFITEKFPVCHQTYIEGHLSATGSTFRMKVLVV...
>MyProject_gene_0001_allele_002
MKLLVLGLPGAGKGTQAQFITEKFPVCHQTYIEGHLSATGSTFRMKVLVV...
>MyProject_gene_0002_allele_001
MSTDFVIADRGEIDTPDEHVKKVLVEEGVDILVFGYPAFYYTPQVIELFTR...
```

#### 3. Allele Names Mapping
**File:** `MyProject_cdhit_allele_names.tsv`
```
Gene_ID	Allele_ID	Original_ID	Strain
MyProject_gene_0001	MyProject_gene_0001_allele_001	strain001|protein001	strain001
MyProject_gene_0001	MyProject_gene_0001_allele_002	strain002|protein045	strain002
MyProject_gene_0002	MyProject_gene_0002_allele_001	strain001|protein002	strain001
```

#### 4. Core Genome Summary
**File:** `core_genome_summary.json`
```json
{
  "total_genes": 1534,
  "core_genes": 987,
  "accessory_genes": 547,
  "core_percentage": 64.3,
  "total_strains": 50,
  "core_threshold": 0.95
}
```

#### 5. Functional Annotations
**File:** `annotations_merged.tsv`
```
Gene_ID	Protein_Name	GO_Terms	KEGG_Pathway	EC_Number	Function
MyProject_gene_0001	DNA polymerase III	GO:0003887;GO:0006260	ko00230	2.7.7.7	DNA replication
MyProject_gene_0002	RNA polymerase	GO:0003899;GO:0006351	ko03020	2.7.7.6	Transcription
```

#### 6. Heaps' Law Results
**File:** `MyProject_heaps_law_results.csv`
```
Genomes,Pan_genes,Core_genes,Accessory_genes,Heaps_kappa,Heaps_gamma
5,2345,1200,1145,0.65,0.42
10,3456,1150,2306,0.67,0.43
15,4123,1120,3003,0.68,0.44
```

---

## Clustering Methods

### CD-HIT (Recommended for Most Cases)
**Pros:**
- Fast execution
- Low memory requirements
- Well-established algorithm
- Good for large datasets (>100 genomes)

**Cons:**
- Only considers sequence similarity
- May miss distant homologs

**Best for:** Standard pangenome analysis, large datasets

**Usage:**
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method cdhit \
    --cdhit_identity 0.70 \
    --cdhit_coverage 0.80 \
    --threads 16
```

### SwiftOrtho (Best for Orthology)
**Pros:**
- Identifies orthologous groups
- Considers evolutionary relationships
- Accurate for cross-species analysis

**Cons:**
- Slower than CD-HIT
- Higher memory usage
- Requires Python 2.7

**Best for:** Cross-species comparisons, orthology analysis

**Usage:**
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method swiftortho \
    --min_similarity 0.65 \
    --swiftortho_coverage 0.75
```

### Foldseek (Structure-Based)
**Pros:**
- Identifies structural homologs
- Can detect distant relationships
- Independent of sequence similarity

**Cons:**
- Requires PDB structure files
- Computationally intensive
- Limited by structure availability

**Best for:** Structure-function relationships, remote homology

**Usage:**
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method foldseek \
    --pdb_directory ./pdb_files \
    --foldseek_sensitivity 9.5 \
    --foldseek_threads 10
```

### Compare All Methods
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method all \
    --threads 16
```

---

## Troubleshooting

### Common Issues

#### 1. Out of Memory Error
**Problem:** `java.lang.OutOfMemoryError`

**Solution:**
```bash
# Increase Nextflow memory
export NXF_OPTS='-Xms1g -Xmx4g'

# Or reduce batch size
nextflow run pangenome_uniprot.nf \
    --uniprot_batch_size 200 \
    --threads 8
```

#### 2. Conda Environment Not Found
**Problem:** `conda: command not found`

**Solution:**
```bash
# Add conda to PATH
export PATH="/path/to/conda/bin:$PATH"

# Or use full path in detect_conda.sh
echo 'export PATH="/path/to/conda/bin:$PATH"' > bin/detect_conda.sh
```

#### 3. UniProt API Rate Limiting
**Problem:** `429 Too Many Requests`

**Solution:**
```bash
# Reduce concurrent workers
nextflow run pangenome_uniprot.nf \
    --uniprot_batch_size 200 \
    --uniprot_max_workers 4
```

#### 4. Foldseek Database Creation Failed
**Problem:** `Database files not found`

**Solution:**
```bash
# Verify PDB files exist
ls -lh pdb_directory/*.pdb | head

# Check Foldseek installation
foldseek version

# Manually create database
foldseek createdb pdb_directory/ testDB
```

#### 5. Empty Output Files
**Problem:** No genes in core genome

**Solution:**
```bash
# Lower core genome threshold
# Edit the analysis script or check your input quality

# Verify input files are not empty
for file in input_fasta/*.faa; do
    echo "$file: $(grep -c '>' $file) sequences"
done
```

### Debug Mode
```bash
# Run with debug output
nextflow run pangenome_uniprot.nf \
    -with-trace \
    -with-report debug_report.html \
    -with-timeline debug_timeline.html \
    -process.echo true
```

### Check Logs
```bash
# View execution logs
less .nextflow.log

# Check specific process logs
cat work/*/*/.command.log
```

---

## Performance Optimization

### For Large Datasets (>100 Genomes)
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method cdhit \
    --cdhit_mode fast \
    --threads 32 \
    --uniprot_batch_size 500 \
    --uniprot_max_workers 12
```

### For Small Datasets (<20 Genomes)
```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method all \
    --cdhit_mode accurate \
    --threads 8
```

### Resource Requirements
| Genomes | Method | CPU | RAM | Time |
|---------|--------|-----|-----|------|
| 10 | CD-HIT | 4 | 8GB | 15min |
| 50 | CD-HIT | 8 | 16GB | 1h |
| 100 | CD-HIT | 16 | 32GB | 3h |
| 50 | SwiftOrtho | 8 | 24GB | 4h |
| 50 | Foldseek | 16 | 32GB | 6h |

---

## Citation

If you use this pipeline, please cite:
- **CD-HIT:** Fu et al. (2012) Bioinformatics
- **SwiftOrtho:** Hu & Friedberg (2019) GigaScience
- **Foldseek:** van Kempen et al. (2023) Nature Biotechnology
- **UniProt:** The UniProt Consortium (2023) Nucleic Acids Research

---

## Support

For issues and questions:
- Check the [Troubleshooting](#troubleshooting) section
- Review logs in `.nextflow.log` and `work/` directory
- Open an issue on GitHub (if applicable)

---

## License

[Specify your license here]

---

## Version History

**v1.0.0** - Initial release
- Multi-method clustering support
- UniProt annotation integration
- Comprehensive analysis modules
