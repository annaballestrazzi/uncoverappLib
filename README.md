# uncoverappLib

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17524340.svg)](https://doi.org/10.5281/zenodo.17524340)

**Interactive web application for clinical assessment and annotation of coverage gaps in target genes.**

Read more: [bioRxiv paper](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1)

---

## 📋 Table of Contents

* [Quick Start Guide](#-quick-start-guide)
* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Introduction](#introduction)
* [Usage Modes](#usage-modes)
  - [Interactive Mode (Shiny App)](#interactive-mode-shiny-app)
  - [Batch Mode (Command-Line)](#batch-mode-command-line)
* [Input Files](#input-files)
* [Examples](#examples)

---

## 🚀 Quick Start Guide

**Get started in 3 steps:**

### Step 1: Install Package
```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("annaballestrazzi/uncoverappLib")
```

### Step 2: Download Annotation Files
```r
library(uncoverappLib)

# Setup annotation cache (run once)
setup_uncoverapp()

# Verify annotations are ready
check_annotations()
```

### Step 3: Launch App
```r
# Launch interactive app
run.uncoverapp()
```

---

## Prerequisites

- **R** ≥ 4.0.0
- **Java** installed (for Excel export)
- **Annotation files** (automatically downloaded via `setup_uncoverapp()`)
  - `sorted_hg19.bed.gz` + index (hg19 genome)
  - `sorted_hg38.bed.gz` + index (hg38 genome)
  - Files from [Zenodo](https://zenodo.org/records/17524340)

---

## Installation

### From GitHub (Development Version)
```r
install.packages("devtools")
devtools::install_github("annaballestrazzi/uncoverappLib")
```

### Verify Installation
```r
library(uncoverappLib)
packageVersion("uncoverappLib")  # Should show 1.10.0 or higher
```

---

## Introduction

**unCOVERApp** is an interactive web application for graphical inspection of sequence coverage within gene regions, designed for clinical diagnostics.

### Key Features

✅ **Coverage Analysis**
- Visual inspection of coverage gaps at base-pair resolution
- Multiple filtering modes (gene, chromosome, region, genome-wide)
- Support for both BAM and BED coverage files
- Handles both 0-based (BED) and 1-based (IGV/VCF) coordinate systems

✅ **Clinical Annotations**
- Integrates [dbNSFP v4.0](https://sites.google.com/site/jpopgen/dbNSFP) annotations
- Prediction scores: MutationAssessor, SIFT, CADD, M-CAP, Polyphen2
- Population frequencies: gnomAD
- Clinical databases: ClinVar, OMIM
- HGVS nomenclature

✅ **Statistical Tools**
- Binomial probability calculator
- [Maximum credible allele frequency (maxAF)](http://cardiodb.org/allelefrequencyapp/) for rare diseases
- Coverage quality metrics

✅ **Export Options**
- Formatted Excel files with conditional coloring
- Publication-ready plots with Gviz
- Tab-separated text files for pipelines

---

## Usage Modes

### Interactive Mode (Shiny App)

**Best for:** Exploratory analysis, single samples, visual inspection

```r
library(uncoverappLib)

# Launch in browser (recommended)
run.uncoverapp()

# Or specify window
run.uncoverapp(where = "browser")  # Default browser
run.uncoverapp(where = "window")   # RStudio window
run.uncoverapp(where = "viewer")   # RStudio viewer pane
```

**App Features:**
- Real-time filtering and annotation
- Interactive plots (Gviz gene tracks)
- Binomial probability calculator
- maxAF calculator
- Excel export with formatting

**Filtering Options:**
1. **By Gene**: Enter HGNC symbol (e.g., `BRCA1`, `TP53`)
2. **By Chromosome**: Select chr1-chr22, chrX, chrY, chrM
3. **By Region**: Genomic coordinates (e.g., `chr17:41196312-41277500`)
4. **All Chromosomes**: Genome-wide analysis

---

### Batch Mode (Command-Line)

**Best for:** Automated pipelines, multiple samples, reproducibility

#### Workflow

```
Input Files → buildInput() → Coverage BED → buildAnnotation() → Annotated Variants
```

#### Step 1: Generate Coverage Data

Process BAM or BED coverage files:

```r
library(uncoverappLib)

# Example: Process BAM files for gene panel
buildInput(
  geneList = "gene_panel.txt",      # Gene list (one per line)
  sampleList = "samples.list",      # BAM/BED file paths
  genome = "hg38",                   # Reference: hg19 or hg38
  chromosome_notation = "chr",       # "chr" or "number"
  type_input = "genes",              # "genes" or "target" (BED)
  type_coverage = "bam",             # "bam" or "bed"
  outDir = "./results",
  MAPQ.min = 20,                     # Minimum mapping quality
  base.quality = 20,                 # Minimum base quality
  input_coord_system = "0-based"     # For BED input: "0-based" or "1-based"
)

# Output:
# - results/output/DATE.bed (coverage matrix)
# - results/output/DATE_statistical_summary.txt
```

**Input file formats:**

`gene_panel.txt`:
```
BRCA1
BRCA2
TP53
PTEN
```

`samples.list`:
```
/path/to/patient1.bam
/path/to/patient2.bam
/path/to/patient3.bam
```

**Or use target BED file:**

`targets.bed`:
```
chr17	41196312	41277500	BRCA1
chr13	32889617	32973809	BRCA2
```

```r
buildInput(
  geneList = "targets.bed",
  type_input = "target",   # Use BED regions
  # ... other parameters
)
```

---

#### Step 2: Annotate Low-Coverage Positions

Identify variants at positions with low coverage:

```r
# Annotate positions below threshold
buildAnnotation(
  sample_data = "results/output/Mon_Nov_11_2024.bed",
  target_sample = "patient1",          # Column name in BED file
  coverage_threshold = 20,              # Coverage cutoff
  genome = "hg38",
  output_intersect = "patient1_lowcov.tsv",
  output_formatted = "patient1_lowcov.xlsx"
)

# Output files:
# - patient1_lowcov.tsv (tab-separated, all columns)
# - patient1_lowcov.xlsx (formatted Excel with colors)
```

**Excel file includes:**
- 🔴 Red highlighting: Pathogenic predictions, ClinVar entries
- 🟢 Green highlighting: Benign predictions
- 🟡 Yellow highlighting: High-impact variants

---

#### Batch Processing Example

Process multiple samples in a cohort:

```r
# Get coverage file from buildInput
coverage_file <- "results/output/Mon_Nov_11_2024.bed"

# List of samples (column names in coverage file)
samples <- c("patient1", "patient2", "patient3")

# Process each sample
for (sample in samples) {
  cat("Processing:", sample, "\n")
  
  buildAnnotation(
    sample_data = coverage_file,
    target_sample = sample,
    coverage_threshold = 20,
    genome = "hg38",
    output_formatted = paste0("output/", sample, "_annotated.xlsx")
  )
}
```

---

## Input Files

### For Processing (buildInput)

**Gene List** (`genes.txt`):
```
BRCA1
BRCA2
TP53
```

**Target BED** (`targets.bed`):
```
chr17	41196312	41277500	BRCA1
chr13	32889617	32973809	BRCA2
chr17	7571720	7590868	TP53
```

**Sample List** (`samples.list`):
```
/full/path/to/sample1.bam
/full/path/to/sample2.bam
```

**Coverage BED** (if using `type_coverage = "bed"`):
```
chr1	100	200	50
chr1	200	300	45
chr1	300	400	30
```
Format: chromosome, start (0-based), end, coverage depth

---

### For Interactive App

**Coverage BED** (`.bed`, `.bedGraph`, `.bed.gz`):
```
chromosome	start	end	sample_1	sample_2
chr1	100	200	50	48
chr1	200	300	45	42
```

**Configuration:**
- Reference genome: hg19 or hg38
- Coverage threshold: e.g., 20x
- Coordinate system: 0-based (BED) or 1-based (IGV)
- Chromosome notation: "chr1" or "1"

---

## Examples

### Example 1: Quick Test with Included Data
```r
library(uncoverappLib)

# Get example files (included in package)
gene_file <- system.file("Test/example_gene.txt", package = "uncoverappLib")
sample_list <- system.file("Test/example.list", package = "uncoverappLib")

# Check what type of files are in the list
readLines(sample_list)  # Shows BED or BAM files

# Process example - adjust type_coverage based on your files
buildInput(
  geneList = gene_file,
  sampleList = sample_list,
  genome = "hg19",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",           # ✅ Change to "bed" if using BED files
  input_coord_system = "0-based",  # ✅ Specify for BED files
  outDir = tempdir()
)

# Or if using BAM files:
# buildInput(
#   geneList = gene_file,
#   sampleList = sample_list,
#   genome = "hg19",
#   chromosome_notation = "chr",
#   type_input = "genes",
#   type_coverage = "bam",
#   MAPQ.min = 20,
#   base.quality = 20,
#   outDir = tempdir()
# )

# Check output
list.files(file.path(tempdir(), "output"))
```

---

### Example 2: Clinical Gene Panel (Cardiology)
```r
# Create gene list
cardio_genes <- c("MYH7", "MYBPC3", "TNNT2", "TNNI3", "TPM1", "ACTC1")
writeLines(cardio_genes, "cardio_panel.txt")

# Create BAM list
bam_files <- c(
  "/data/patient_001.bam",
  "/data/patient_002.bam",
  "/data/patient_003.bam"
)
writeLines(bam_files, "bams.list")

# Process cohort with BAM files
buildInput(
  geneList = "cardio_panel.txt",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",           # ✅ BAM files
  MAPQ.min = 30,                   # ✅ Required for BAM
  base.quality = 20,               # ✅ Required for BAM
  outDir = "./cardio_results"
)

# Annotate low-coverage for each patient
coverage_file <- list.files("cardio_results/output", pattern = "\\.bed$", full.names = TRUE)[1]

for (patient in c("count_patient_001", "count_patient_002", "count_patient_003")) {
  buildAnnotation(
    sample_data = coverage_file,
    target_sample = patient,
    coverage_threshold = 30,
    genome = "hg38",
    output_formatted = paste0("cardio_results/", patient, "_lowcov.xlsx")
  )
}
```

---

### Example 3: Using Pre-computed BED Coverage Files
```r
# If you already have coverage BED files (faster!)
bed_files <- c(
  "/data/sample1_coverage.bed",
  "/data/sample2_coverage.bed"
)
writeLines(bed_files, "coverage.list")

# BED format example (0-based):
# chr1  100  200  50
# chr1  200  300  45

buildInput(
  geneList = "genes.txt",
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",           # ✅ BED coverage files
  input_coord_system = "0-based",  # ✅ Specify coordinate system
  outDir = "./results"
)
```

---

### Example 4: Interactive App - Region-Specific Analysis (BRCA1 Exon 11)
```r
# Launch app for manual analysis
run.uncoverapp()

# In the app interface:
# 1. Load your coverage file (WGS.bed.gz)
# 2. Select "Filter by: Region coordinates"
# 3. Enter region: chr17:41243452-41246877 (BRCA1 exon 11)
# 4. Set threshold: 20x
# 5. Click "Calculate Low Coverage"
# 6. Click "Calculate Annotations"
# 7. Explore binomial probability for specific positions
# 8. Export results to Excel
```

---

## Troubleshooting

### Common Issues

**"No annotation files found"**
```r
# Re-run setup
setup_uncoverapp()
check_annotations()
```

**"Gene not found"**
- Use official HGNC symbols (e.g., `TP53` not `p53`)
- Check genome version matches your data (hg19 vs hg38)

**"No coverage data"**
- Verify BAM/BED chromosome naming matches ("chr1" vs "1")
- Check coordinate system (0-based vs 1-based)

**"Cannot open file"**
- Use absolute paths in `.list` files
- Check file permissions

---

## Citation

If you use uncoverappLib in your research, please cite:

```
Ballestraz et al. (2020). unCOVERApp: a web application for 
clinical assessment of sequence coverage at the base-pair level.
bioRxiv. doi: 10.1101/2020.02.10.939769
```

---

## Documentation

- **Vignettes**: `browseVignettes("uncoverappLib")`
- **Function help**: `?buildInput`, `?buildAnnotation`
- **GitHub Issues**: [Report bugs](https://github.com/annaballestrazzi/uncoverappLib/issues)

---

## License

GPL-3 | file LICENSE

---

## Version History

**v1.10.0** (Current)
- ✅ Improved gene coordinate retrieval (hg19/hg38 specific)
- ✅ Fixed filter modes (gene/chromosome/region/all)
- ✅ Added BED coverage file support
- ✅ Coordinate system selection (0-based/1-based)
- ✅ Enhanced maxAF calculator compatibility
- ✅ Parameter name updates: `sampleList` (was `bamList`), `chromosome_notation` (was `type_bam`)

**v0.1.0** (Initial release)
- Interactive Shiny app
- BAM file processing
- dbNSFP annotation integration
