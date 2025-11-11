# uncoverappLib

This repository is home of R package uncoverappLib launching *unCOVERApp*, 
a web application for clinical assessment and annotation of coverage gaps in
target genes. Read more about unCOVERApp on [biorxiv](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1)

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17524340.svg)](https://doi.org/10.5281/zenodo.17524340)

# Table of content

* [Prerequisites](#prerequisites)
* [Installation](#installation)
* [Introduction](#introduction)
* [Download annotation files](#download-annotation-files)
* [Input](#input)
* [Usage](#usage)

## Prerequisites

This app requires following dependencies:

- R >= 4.0.0
- Java installed 
- Annotation files (`sorted_hg19.bed.gz`, `sorted_hg19.bed.gz.tbi`, `sorted_hg38.bed.gz`, `sorted_hg38.bed.gz.tbi`) downloadable via [Zenodo](https://zenodo.org/records/17524340) with the `getAnnotationFiles()` function of *uncoverappLib*.

## Installation

To install this package, start R and enter: 

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install("uncoverappLib")     
```

The development version can be installed directly from Github:

```r
install.packages("devtools")
devtools::install_github("Manuelaio/uncoverappLib")
```

## Introduction

The rapid spread of NGS-technology has been coupled since its beginning with 
development of bioinformatic tools for data analysis and interpretation. 
However, despite increasing accuracy of available approaches, the need to 
assess sequencing quality of the analysis targets at the base-pair resolution 
poses growing challenges especially in the clinical settings.  

In diagnostics indeed, meticulous investigation of every single target base is 
often required to exclude that pathogenic events across the gene of interest 
may be missed due to uneven sequence coverage.

**unCOVERApp** is an interactive web-application for graphical inspection of sequence coverage within gene regions.

unCOVERApp highlights low coverage genomic positions, according to the coverage
threshold specified by the user, providing *dbNSFP-based annotation*s for 
clinical assessment of low coverage regions. 

It implements basic statistical tools such as binomial probability calculation 
that genomic positions are adequately covered, and 
[maximum credible allele frequency](http://cardiodb.org/allelefrequencyapp/). 

## Download annotation files

To associate low-coverage sites with functional and clinical annotations, 
unCOVERApp uses [dbNSFP](https://sites.google.com/site/jpopgen/dbNSFP) 
version 4.0 stored in four files:

**For hg19 genome:**
* `sorted_hg19.bed.gz`: a genomically-sorted, TABIX-indexed, BGZipped BED file containing selected columns from dbNSFP version v4.0
* `sorted_hg19.bed.gz.tbi`: TABIX-indexed file for hg19

**For hg38 genome:**
* `sorted_hg38.bed.gz`: a genomically-sorted, TABIX-indexed, BGZipped BED file containing selected columns from dbNSFP version v4.0
* `sorted_hg38.bed.gz.tbi`: TABIX-indexed file for hg38

Those files are stored on Zenodo at the following
[link](https://zenodo.org/records/17524340) for downloading.

The annotation files enclose prediction scores (MutationAssessor, SIFT, CADD, 
M-CAP, Polyphen2-HVAR), allele frequencies observed in 
gnomAD data, dbsnp accession number, HGVS notations and clinical annotation 
information from ClinVar and OMIM. Loading these files allows the annotation 
of each low coverage genomic position user-defined.

Run following commands to correctly create a unique cache for **uncoverappLib** 
using **BiocFileCache** package:

```r
library(uncoverappLib)
getAnnotationFiles()
```

The annotation files will be downloaded and cached in:
```r
rappdirs::user_cache_dir("uncoverapp")
```

## Input

### Processing and Statistical Summary page

For data preprocessing, unCOVERApp accepts:

**Input filtering file** (gene list or target regions):
- Text file with `.txt` extension containing HGNC official gene names (one per row)
- BED file with target regions (chromosome, start, end coordinates)

**Genomic data file** (coverage information, BAM or BED):
- `.list` file containing absolute paths to BAM files (one per row)
- `.list` file containing absolute paths to BED coverage files (one per row)

**Note**: In the output, the first file in the list corresponds to sample 1, and so forth.

### Coverage Analysis page

For coverage analysis and visualization, unCOVERApp accepts:

**Input coverage file**:
- BED file (`.bed`, `.bedGraph`) with columns: chromosome, start, end, coverage (by sample)
- Compressed formats: `.bed.gz`, `.bedGraph.gz`, `.zip`
- CSV file with coverage data

**Filtering options**:
- Filter by **Gene name**: Enter HGNC gene symbol (e.g., POLG)
- Filter by **Chromosome**: Select specific chromosome (chr1-chr22, chrX, chrY)
- Filter by **All chromosomes**: Analyze entire genome
- Filter by **Region coordinates**: Specify genomic region (e.g., chr7:37453745-45627643)

**Configuration parameters**:
- Reference genome: hg19 or hg38
- Coverage threshold: Minimum coverage depth (default: 20x)
- Coordinate system: 0-based (BED standard) or 1-based (IGV, VCF)
- Chromosome notation: with 'chr' prefix or numeric only
- Minimum mapping quality (MAPQ)
- Minimum base quality

For more details on working with unCOVERApp see Vignette.
## Examples

Example data files are provided in `inst/Test/`:
- `example_gene.txt` - Sample gene list (POLG)
- `example.list` - Sample BAM file list
- `example_POLG.bam` - Example BAM file

### Quick Start Example
```r
library(uncoverappLib)

# Get example files
gene_file <- system.file("Test/example_gene.txt", package = "uncoverappLib")
bam_list <- system.file("Test/example.list", package = "uncoverappLib")

# Process example data
buildInput(
  geneList = gene_file,
  bamList = bam_list,
  genome = "hg19",
  type_bam = "chr",
  type_input = "genes",
  type_coverage = "bam",
  outDir = tempdir()
)
```


## Usage

uncoverappLib supports **two usage modes**:

### Interactive Mode (Shiny App)

Load library and launch the interactive web application:
```r
library(uncoverappLib)

# Launch in RStudio window
run.uncoverapp(where = "window")

# Or in default browser
run.uncoverapp(where = "browser")

# Or in RStudio viewer pane
run.uncoverapp(where = "viewer")
```

**Features:**
- Real-time coverage analysis and visualization
- Interactive filtering by gene, chromosome, or region
- Dynamic annotation with dbNSFP
- Binomial probability calculator
- maxAF calculator for rare diseases
- Visual plots with Gviz
- Export to formatted Excel files

For more details, see the package vignettes and examples above.
---

### Batch Mode (Command-Line)

For automated pipelines and batch processing of multiple samples.

#### Step 1: Generate Coverage Data

Process BAM files or BED coverage files for target genes:
```r
library(uncoverappLib)

# Process BAM files
buildInput(
  geneList = "genes.txt",           # Gene list or BED file
  bamList = "samples.list",          # List of BAM/BED files
  genome = "hg19",                   # hg19 or hg38
  type_bam = "chr",                  # Chromosome notation
  type_input = "genes",              # "genes" or "target"
  type_coverage = "bam",             # "bam" or "bed"
  outDir = "./output",
  MAPQ.min = 20,
  base.quality = 20
)

# Output files:
# - output/output/DATE.bed (coverage data)
# - output/output/DATE_statistical_summary.txt
```

**Input file formats:**

`genes.txt`:
```
BRCA1
BRCA2
TP53
```

`samples.list`:
```
/path/to/sample1.bam
/path/to/sample2.bam
```

#### Step 2: Annotate Low-Coverage Positions

Identify and annotate all genomic positions below coverage threshold:
```r
# Annotate low-coverage positions genome-wide
annotate_all_lowcov(
  sample_data = "output/output/Mon_Nov_11_2024.bed",
  target_sample = "sample1",
  coverage_threshold = 20,
  genome = "hg19",
  output_intersect = "annotated.tsv",
  output_formatted = "annotated.xlsx"
)

# Output files:
# - annotated.tsv (tab-separated data)
# - annotated.xlsx (formatted Excel with conditional coloring)
```

#### Batch Processing Multiple Samples
```r
# Process cohort
samples <- c("sample1", "sample2", "sample3")
coverage_file <- "output/output/Mon_Nov_11_2024.bed"

for (sample in samples) {
  annotate_all_lowcov(
    sample_data = coverage_file,
    target_sample = sample,
    coverage_threshold = 20,
    genome = "hg19",
    output_formatted = paste0(sample, "_annotated.xlsx")
  )
}
```

For detailed batch processing examples, see the package vignettes.
```
---
User can define where uncoverapp will be launched with the `where` option:

- `browser` option will open `uncoverapp` in your default browser
- `viewer` option will open `uncoverapp` in RStudio viewer pane
- `window` option will open `uncoverapp` in RStudio window

For more details on working with uncoverapp see Vignette or [Documentation.pdf](https://github.com/Manuelaio/unCOVERApp/blob/master/Documentation.pdf) on Github.
