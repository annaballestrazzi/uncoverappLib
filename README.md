# uncoverappLib

This repository is home of R package uncoverappLib launching *unCOVERApp*, 
a web application for clinical assessment and annotation of coverage gaps in
target genes. Read more about unCOVERApp on [biorxiv](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1)

[![DOI](https://zenodo.org/badge/254597958.svg)](https://zenodo.org/badge/latestdoi/254597958)

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
- Annotation files (`sorted_hg19.bed.gz`, `sorted_hg19.bed.gz.tbi`, `sorted_hg38.bed.gz`, `sorted_hg38.bed.gz.tbi`) downloadable via [Zenodo](https://zenodo.org/record/3747448#.XpBmnVMzbOR) with the `getAnnotationFiles()` function of *uncoverappLib*.

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

Major informations about unCOVERApp R dependences are 
[here](https://github.com/Manuelaio/test_dependence).

[![Build Status](https://travis-ci.com/Manuelaio/test_dependence.svg?branch=master)](https://travis-ci.com/Manuelaio/test_dependence)

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
[link](https://zenodo.org/record/3747448#.XpBmnVMzbOR) for downloading. 

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
## Test
I test automatici di questo package si trovano nella cartella tests/testthat, originariamente sviluppati nella cartella inst/Test
- example.list
- example_gene.txt


## Usage

Load library and set up R environment with annotation files as following. 
The way to launch unCOVERApp is with the `run.uncoverapp(where="window")` function:

```r
library(uncoverappLib)
run.uncoverapp(where="window")
```

User can define where uncoverapp will be launched with the `where` option:

- `browser` option will open `uncoverapp` in your default browser
- `viewer` option will open `uncoverapp` in RStudio viewer pane
- `window` option will open `uncoverapp` in RStudio window

For more details on working with uncoverapp see Vignette or [Documentation.pdf](https://github.com/Manuelaio/unCOVERApp/blob/master/Documentation.pdf) on Github.
