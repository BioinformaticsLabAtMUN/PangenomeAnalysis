# Pangenome Analysis Pipeline

A comprehensive Nextflow pipeline for pangenome analysis using multiple clustering methods (CD-HIT, SwiftOrtho, Foldseek) with functional annotation and core genome analysis.

## Table of Contents
- [Overview](#overview)
- [Prerequisites](#prerequisites)
- [Installation](#installation)
- [Input Files](#input-files)
- [Configuration](#configuration)
- [Running the Pipeline](#running-the-pipeline)
- [Output Files](#output-files)
- [Detailed Examples](#detailed-examples)
- [Troubleshooting](#troubleshooting)

## Overview

This pipeline performs:
1. **Sequence consolidation** - Merges multiple FASTA files from different strains
2. **Clustering** - Groups similar proteins using CD-HIT, SwiftOrtho, or Foldseek
3. **Pangenome table generation** - Creates presence/absence matrices
4. **Core genome analysis** - Identifies core genes present in all strains
5. **Functional annotation** - Annotates proteins using UniProt or NCBI databases
6. **GO term analysis** - Performs Gene Ontology enrichment and clustering
7. **Visualization** - Generates plots and summary statistics

## Prerequisites

### Software Requirements
- Nextflow (≥21.04.0)
- Python (≥3.8)
- CD-HIT (≥4.8.1)
- Foldseek (for structure-based clustering)
- Conda or Mamba (for SwiftOrtho environment)

### Python Dependencies
Install required Python packages:
```bash
pip install -r requirements.txt
```

**Expected requirements.txt contents:**
```
scipy>=1.7.0
numpy>=1.21.0
pandas>=1.3.0
biopython>=1.79
requests>=2.26.0
goatools>=1.2.3
matplotlib>=3.4.0
seaborn>=0.11.0
```

## Input Files

### Required Files

#### 1. Protein FASTA Files
- **Location**: Directory containing one FASTA file per strain
- **Format**: `.fasta`, `.faa`, or `.fa` extensions
- **Example structure**:
```
input_directory/
├── strain_001.faa
├── strain_002.faa
├── strain_003.faa
└── ...
```

**Sample FASTA file** (`strain_001.faa`):
```
>A0A4P6TQ31 hypothetical protein
MKLVFSLSLLSALVAAHAAEQKTEDLSFKIGEMESRLTKYSSSIMMSDIFKQIQKFYKQ
TNISKIESLFKQIESIFDKFMEAIPDFVKKGIIDEEEDSLLQ
>Q8RN19 DNA polymerase
MTKQILRELQVAQPEDQVEAILKFVVDQNWKDKELATELLKEIEDFKDPAYITIQDIGE
SSSLREKVDKFVQAIQKDRQVAEAKKRRKSGFRFDRQGKLYTDEKSNTVNVKNFKK
>P12345 RNA polymerase subunit
MSEIEQLKELEVQARERMQQQLQSMEQRARRRAEQARKQEEAERLQQALQDTMEREQRR
MNEEMARERLSQLQMEMRQQQVQLQVQNQRVQMQQQMQQQMQQQMQQQPQQQMMQQQ
```

#### 2. Proteome Metadata File
- **File**: `proteome_metadata.tsv` or similar
- **Format**: Tab-separated values with headers
- **Example**:
```tsv
Strain_ID	Organism	Assembly_Accession	Genome_Size	GC_Content	Source
strain_001	Streptomyces coelicolor	GCA_000203835.1	8667507	72.1	NCBI
strain_002	Streptomyces avermitilis	GCA_000069185.1	9025608	70.7	NCBI
strain_003	Streptomyces venezuelae	GCA_000583135.1	8227074	73.2	NCBI
```

#### 3. Shared Headers File (Optional)
- **File**: `shared_headers.tsv`
- **Purpose**: Map synonymous protein IDs across strains
- **Format**: Tab-separated, one line per synonym group
- **Example**:
```tsv
A0A4P6TQ31	WP_003959657.1	YP_001825635.1
Q8RN19	WP_011028374.1	NP_624051.1
```

### Optional Files (for Foldseek)

#### PDB Structure Files
- **Location**: Directory containing AlphaFold PDB files
- **Format**: `.pdb` files named by protein accession
- **Example**:
```
pdb_files/
├── A0A4P6TQ31.pdb
├── Q8RN19.pdb
└── P12345.pdb
```

#### PDB Mapping File
- **File**: `pdb_mapping.tsv`
- **Format**:
```tsv
Protein_ID	PDB_File	Confidence
A0A4P6TQ31	A0A4P6TQ31.pdb	95.2
Q8RN19	Q8RN19.pdb	87.4
```

### Optional Files (for Revigo Analysis)

#### Revigo GO Term Tables
- **Location**: `revigo/` directory
- **Files**: 
  - `Revigo_BP_Table.tsv` (Biological Process)
  - `Revigo_MF_Table.tsv` (Molecular Function)
- **Format**:
```tsv
Representative	Name	Dispensability	Frequency	Plot_X	Plot_Y
GO:0006810	transport	0.000	12.5	-4.2	2.1
GO:0055085	transmembrane transport	0.150	8.3	-3.8	1.9
```

## Configuration

### Basic Configuration

Edit the main parameters in `pangenome_uniprot.nf` or create a custom config file:

```groovy
// Basic settings
params.name_prefix = 'MyProject'
params.baseDir = "$PWD"
params.raw_input_directory = "${params.baseDir}/input_fastas"
params.proteome_metadata_file = "${params.baseDir}/metadata.tsv"

// Clustering method selection
params.clustering_method = "cdhit"  // Options: "cdhit", "swiftortho", "foldseek", "all"

// Database for annotation
params.database = "uniprot"  // Options: "ncbi", "uniprot"

// CD-HIT parameters
params.cdhit_identity = 0.65   // Sequence identity threshold (0-1)
params.cdhit_coverage = 0.75   // Alignment coverage threshold (0-1)
params.threads = 12            // Number of CPU threads

// Foldseek parameters (if using structural clustering)
params.foldseek_min_seq_id = 0.3
params.foldseek_threads = 10
```

### Advanced Configuration Options

```groovy
// UniProt annotation optimization
params.uniprot_batch_size = 500    // Batch size for API requests
params.uniprot_max_workers = 8     // Concurrent API workers

// SwiftOrtho parameters
params.swiftortho_evalue = "1e-5"
params.swiftortho_coverage = 0.75

// Core genome analysis
params.annotate_scope = "core"     // Options: "core", "all"
params.max_alternatives_per_cluster = 3

// Visualization
params.max_strains_missing_plot = 50  // Max strains in missing genes plot
```

## Running the Pipeline

### Step 1: Prepare Your Environment

```bash
# Navigate to pipeline directory
cd /path/to/pangenome_pipeline

# Verify Nextflow installation
nextflow -version

# Check Python dependencies
python3 -c "import scipy, numpy, pandas, Bio; print('Dependencies OK')"
```

### Step 2: Organize Input Files

```bash
# Create directory structure
mkdir -p input_fastas
mkdir -p output
mkdir -p cache

# Place your FASTA files in input_fastas/
cp /path/to/your/fastas/*.faa input_fastas/

# Create or copy metadata file
cp /path/to/metadata.tsv .
```

### Step 3: Configure the Pipeline

Edit `pangenome_uniprot.nf` or create `nextflow.config`:

```groovy
params {
    name_prefix = 'Streptomyces'
    raw_input_directory = "${baseDir}/input_fastas"
    proteome_metadata_file = "${baseDir}/metadata.tsv"
    clustering_method = "cdhit"
    database = "uniprot"
    threads = 12
}
```

### Step 4: Run the Pipeline

#### Basic Run (CD-HIT clustering with UniProt annotation)
```bash
nextflow run pangenome_uniprot.nf
```

#### Run with Specific Clustering Method
```bash
# CD-HIT only
nextflow run pangenome_uniprot.nf --clustering_method cdhit

# Foldseek only (requires PDB files)
nextflow run pangenome_uniprot.nf --clustering_method foldseek \
    --pdb_directory ./pdb_files

# Run all methods
nextflow run pangenome_uniprot.nf --clustering_method all
```

#### Run with Custom Parameters
```bash
nextflow run pangenome_uniprot.nf \
    --name_prefix MyGenera \
    --raw_input_directory ./my_fastas \
    --cdhit_identity 0.70 \
    --cdhit_coverage 0.80 \
    --threads 16 \
    --database uniprot
```

#### Resume Failed Run
```bash
nextflow run pangenome_uniprot.nf -resume
```

### Step 5: Monitor Progress

Nextflow will display progress in real-time:
```
N E X T F L O W  ~  version 21.04.0
Launching `pangenome_uniprot.nf` [peaceful_galileo] - revision: abc123

executor >  local (12)
[ab/cd1234] process > consolidateSequences      [100%] 1 of 1 ✔
[ef/gh5678] process > runCDHIT                   [100%] 1 of 1 ✔
[ij/kl9012] process > renameCDHITSequences       [100%] 1 of 1 ✔
[mn/op3456] process > generateCDHITPangenomeTables [100%] 1 of 1 ✔
[qr/st7890] process > analyzeAndValidateCoreGenomeCDHIT [  0%] 0 of 1
```

## Output Files

### Directory Structure

After successful completion:
```
output/
├── uniprot_output/
│   ├── cdhit/
│   │   ├── cdhit_output.faa          # Clustered representative sequences
│   │   ├── cdhit_output.faa.clstr    # Cluster assignments
│   │   ├── renamed_sequences/
│   │   │   └── allele_names.tsv      # Original ID → Allele ID mapping
│   │   ├── pangenome/
│   │   │   ├── Strep_strain_by_allele.npz    # Allele presence/absence matrix
│   │   │   ├── Strep_strain_by_gene.npz      # Gene presence/absence matrix
│   │   │   ├── allele_list.txt               # List of all alleles
│   │   │   └── gene_list.txt                 # List of all genes
│   │   ├── core_genome/
│   │   │   ├── core_genes.txt                # Core genes (present in all)
│   │   │   ├── core_analysis_report.txt      # Statistics
│   │   │   ├── missing_core_genes.png        # Visualization
│   │   │   └── core_sequences.faa            # Core gene sequences
│   │   ├── annotations/
│   │   │   ├── dominant_alleles_core_summary.tsv   # Core gene annotations
│   │   │   ├── uniprot_annotation_results.tsv      # Full UniProt data
│   │   │   └── annotation_stats.txt                # Coverage statistics
│   │   ├── functional_analysis/
│   │   │   ├── go_term_enrichment.tsv
│   │   │   ├── functional_categories.png
│   │   │   └── pathway_analysis.tsv
│   │   └── heaps_law/
│   │       ├── heaps_analysis.png
│   │       └── pangenome_stats.txt
│   ├── swiftortho/                   # Similar structure for SwiftOrtho
│   └── foldseek/                     # Similar structure for Foldseek
└── consolidated/
    ├── consolidated.faa              # All proteins combined
    ├── shared_headers.tsv            # Synonym mapping
    └── metadata.tsv                  # Strain metadata
```

### Key Output Files Explained

#### 1. Cluster File (`cdhit_output.faa.clstr`)
Shows how sequences cluster together:
```
>Cluster 0
0	450aa, >A0A4P6TQ31... *
1	448aa, >WP_003959657.1... at +/99.33%
2	452aa, >YP_001825635.1... at +/98.89%

>Cluster 1
0	523aa, >Q8RN19... *
1	520aa, >WP_011028374.1... at +/97.52%
```
- `*` indicates the representative sequence
- Percentage shows sequence identity to representative

#### 2. Allele Names File (`allele_names.tsv`)
Maps original protein IDs to standardized allele names:
```tsv
Original_ID	New_ID
A0A4P6TQ31	G000001A1
WP_003959657.1	G000001A2
Q8RN19	G000002A1
WP_011028374.1	G000002A2
```

#### 3. Core Genes File (`core_genes.txt`)
Lists genes present in all strains:
```
G000001
G000015
G000027
G000042
...
```

#### 4. Core Analysis Report (`core_analysis_report.txt`)
```
=== CORE GENOME ANALYSIS REPORT ===

Dataset: Streptomyces pangenome
Total strains: 150
Total genes: 45,892
Core genes: 1,247 (2.72% of pangenome)

Core gene presence:
- Perfect core (100% strains): 1,247 genes
- Soft core (≥95% strains): 1,523 genes

Missing core gene statistics:
- Genes missing from 0 strains: 1,247
- Genes missing from 1 strain: 0
- Genes missing from 2 strains: 0

Top strains with most missing core genes:
1. strain_042: 3 missing genes
2. strain_087: 2 missing genes
3. strain_123: 1 missing gene
```

#### 5. Dominant Alleles Annotation (`dominant_alleles_core_summary.tsv`)
Core genes with functional annotations:
```tsv
Gene_ID	Representative_Allele	Protein_Name	Gene_Name	GO_Terms	Organism	Function
G000001	G000001A1	DNA polymerase III subunit alpha	dnaE	GO:0003887;GO:0006260	Streptomyces	DNA replication
G000002	G000002A3	RNA polymerase sigma factor	rpoD	GO:0003677;GO:0006352	Streptomyces	Transcription initiation
G000003	G000003A1	50S ribosomal protein L1	rplA	GO:0003735;GO:0006412	Streptomyces	Translation
```

#### 6. Heaps Law Analysis (`heaps_analysis.png`)
Visualization showing pangenome openness:
- X-axis: Number of genomes sampled
- Y-axis: Number of gene families
- Curve showing core vs. accessory genome growth

#### 7. GO Term Enrichment (`go_term_enrichment.tsv`)
```tsv
GO_ID	GO_Name	Category	Gene_Count	P_Value	FDR
GO:0006260	DNA replication	BP	45	1.2e-15	3.5e-13
GO:0006412	translation	BP	78	2.3e-25	8.9e-23
GO:0003677	DNA binding	MF	112	5.4e-18	1.2e-15
```

## Detailed Examples

### Example 1: Basic CD-HIT Analysis

**Input**: 50 bacterial genome FASTA files

```bash
# Step 1: Setup
mkdir -p strep_analysis/input_fastas
cd strep_analysis

# Step 2: Copy FASTA files
cp /data/genomes/*.faa input_fastas/

# Step 3: Create metadata
cat > metadata.tsv << EOF
Strain_ID	Organism	Source
strain_001	Streptomyces coelicolor	Lab_strain
strain_002	Streptomyces avermitilis	Type_strain
EOF

# Step 4: Run pipeline
nextflow run /path/to/pangenome_uniprot.nf \
    --name_prefix Streptomyces \
    --raw_input_directory ./input_fastas \
    --proteome_metadata_file ./metadata.tsv \
    --cdhit_identity 0.65 \
    --threads 8
```

**Expected Runtime**: ~30 mins for 140 genomes with ~8,000 proteins each

**Expected Output**:
- Consolidated FASTA: ~400,000 total protein sequences
- After clustering: ~45,000 gene families
- Core genome: ~1,200-1,500 genes (present in all 50 strains)
- Output directory: ~5-10 GB

### Example 2: Foldseek Structure-Based Clustering

**Requires**: PDB structure files from AlphaFold

```bash
# Step 1: Prepare PDB files
mkdir -p pdb_files
# Copy or download AlphaFold predictions for your proteins

# Step 2: Run with Foldseek
nextflow run pangenome_uniprot.nf \
    --clustering_method foldseek \
    --pdb_directory ./pdb_files \
    --foldseek_min_seq_id 0.3 \
    --foldseek_threads 10
```

**Output differences**:
- May cluster more distantly related proteins with similar structures
- Useful for identifying functional homologs with low sequence identity

### Example 3: Comprehensive Analysis (All Methods)

```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method all \
    --database uniprot \
    --annotate_scope all \
    --threads 16
```

**Output**: Generates three separate analysis directories:
- `cdhit/` - Sequence-based clustering
- `swiftortho/` - Orthology-based clustering  
- `foldseek/` - Structure-based clustering

Each with independent core genome and annotation results for comparison.

### Example 4: Using Existing Foldseek Results

If you already have Foldseek clustering results:

```bash
nextflow run pangenome_uniprot.nf \
    --clustering_method foldseek \
    --use_existing_foldseek_db true \
    --existing_foldseek_db_path ./existing_db \
    --use_existing_foldseek_clusters true \
    --existing_foldseek_clusters_path ./clusters.tsv
```

## Troubleshooting

### Common Issues

#### 1. "No FASTA files found"
**Problem**: Pipeline can't find input files

**Solution**:
```bash
# Check file extensions
ls input_fastas/*.{fasta,faa,fa}

# Verify directory parameter
nextflow run pangenome_uniprot.nf --raw_input_directory $(pwd)/input_fastas
```

#### 2. "CD-HIT command not found"
**Problem**: CD-HIT not in PATH

**Solution**:
```bash
# Install CD-HIT
conda install -c bioconda cd-hit

# Or build from source
wget https://github.com/weizhongli/cdhit/releases/download/V4.8.1/cd-hit-v4.8.1-2019-0228.tar.gz
tar xzf cd-hit-v4.8.1-2019-0228.tar.gz
cd cd-hit-v4.8.1-2019-0228
make
sudo make install
```

#### 3. "UniProt annotation failed"
**Problem**: API requests timing out or rate limited

**Solution**:
```bash
# Reduce batch size and workers
nextflow run pangenome_uniprot.nf \
    --uniprot_batch_size 200 \
    --uniprot_max_workers 4
```


#### 4. Pipeline Stops Unexpectedly
**Problem**: Process fails without clear error

**Solution**:
```bash
# Enable detailed logging
nextflow run pangenome_uniprot.nf -with-trace -with-report report.html

# Check the trace file
cat trace.txt

# Resume from last checkpoint
nextflow run pangenome_uniprot.nf -resume
```

### Validation Steps

#### Check Input Files
```bash
# Count FASTA files
ls input_fastas/*.faa | wc -l

# Check FASTA format
head -n 20 input_fastas/strain_001.faa

# Verify metadata
cat metadata.tsv | column -t
```

#### Monitor Resource Usage
```bash
# While pipeline is running
top -u $USER          # CPU and memory
df -h                 # Disk space
```

#### Verify Outputs
```bash
# Check core genome size
wc -l output/uniprot_output/cdhit/core_genome/core_genes.txt

# View top annotated genes
head -n 20 output/uniprot_output/cdhit/annotations/dominant_alleles_core_summary.tsv

# Check for errors in log
grep -i error .nextflow.log
```

## Performance Optimization

### For Large Datasets (>100 genomes)

```groovy
// Optimize CD-HIT
params.cdhit_mode = "fast"        // Trade accuracy for speed
params.threads = 32               // Use all available cores

// Optimize annotation
params.uniprot_batch_size = 500
params.uniprot_max_workers = 12

// Limit annotation scope
params.annotate_scope = "core"    // Only annotate core genes
```

### For High-Quality Analysis

```groovy
// Use accurate mode
params.cdhit_mode = "accurate"
```

## Citations

If you use this pipeline, please cite:

---

**Last Updated**: January 2026  
**Pipeline Version**: 2.0  
**Nextflow DSL**: 2
