# uncoverappLib <img src="shiny-dir/www/logo.png" align="right" height="90"/>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17524340.svg)](https://doi.org/10.5281/zenodo.17524340)
[![R ≥ 4.0](https://img.shields.io/badge/R-%3E%3D4.0-blue)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow)](https://opensource.org/licenses/MIT)

> **Interactive coverage analysis and variant annotation for targeted clinical sequencing.**  
> Identify gaps. Annotate variants. Support clinical decisions — without uploading a single file to an external server.

[Read the paper](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1) · [Quick Start](inst/extdata/QUICK_START.md) · [User Guide](inst/extdata/USER_GUIDE.md) · [Citation](#citation)

---

## What uncoverappLib does

uncoverappLib is an R package for the analysis of sequencing coverage in targeted gene panels, whole-exome, and whole-genome experiments. Starting from raw BAM files or pre-computed coverage files, it identifies genomic positions with insufficient coverage, retrieves functional and clinical annotations for all known variants at those positions, and provides statistical tools to support clinical variant interpretation.

The package has two complementary modes of use:

**Interactive app** — a Shiny application that runs locally in your browser. No data leaves your computer. Designed for exploratory analysis, visual inspection of coverage, and interactive filtering.

**Standalone functions** — `buildInput()` and `buildAnnotation()` can be called from the R console or within automated pipelines, making them suitable for batch processing of large cohorts or integration into existing bioinformatics workflows.

The two modes are fully interoperable: the coverage BED file produced by `buildInput()` can be loaded directly into the interactive app for visual exploration, and results obtained interactively can be reproduced exactly using `buildAnnotation()` on the same file.

---

## What you can do

### Processing and Statistical Summary tab

Upload your gene list (or a custom target BED with chromosomal coordinates) and a list of BAM or pre-computed BED coverage files. After processing:

- A **statistical summary** is generated per gene per sample — width-weighted mean and median coverage, total bases below threshold, OMIM disease annotation.
- A **filter panel** appears that lets you select samples and set a coverage threshold. Positions where all selected samples are above the threshold are hidden; values above the threshold are displayed as `>N` rather than the exact number.
- A **Download XLSX** button exports the data: if you applied the filter, only the filtered subset is exported; if you did not, the full coverage table is downloaded.

### Coverage Analysis tab

Load a pre-processed BED file (from `buildInput()` or the Processing tab) and apply filters **before** calculating coverage, so that only the data you care about is processed. The available filter modes are:

- **Gene name** — enter an HGNC symbol; the app retrieves genomic coordinates automatically from the UCSC database
- **Chromosome** — select from a dropdown (chr1–chr22, chrX, chrY, chrM)
- **Genomic region** — enter coordinates in `chr:start-end` format
- **All chromosomes** — analyse the entire file at once

After choosing filter mode and parameters, click *Calculate Low Coverage Regions* to identify positions below the threshold, then *Calculate Annotations on Low Coverage* to annotate them. The **Download Annotations** button (visible only after annotation is complete) exports the annotated table as a colour-coded Excel file.

### Variant annotation

For every low-coverage position, uncoverappLib queries a local tabix-indexed annotation database (dbNSFP/VEP) and returns all known variants at that position, with: ClinVar classification, CADD Phred score, gnomAD allele frequency, MutationAssessor functional impact (H/M/L), SIFT, PolyPhen-2, M-CAP, and HGVS nomenclature from VEP canonical transcripts.

### Maximum credible allele frequency (maxAF)

A calculator for the maximum credible population allele frequency tailored to a specific disease model. Input the inheritance pattern, disease prevalence, allelic and genetic heterogeneity, and penetrance — and the app computes a gene-specific threshold used to filter the annotation table. Download as Excel.

### Binomial probability

Given the observed coverage depth at a position and an expected allele fraction, compute the binomial probability distribution and 95% confidence interval for detecting at least N variant-supporting reads. Includes a visual interpretation indicating whether the current coverage is sufficient.

### Gene coverage plot

A multi-track Gviz visualisation of coverage along a full gene, showing exon/intron structure from the UCSC transcript database and coverage coloured red (below threshold) or blue (above threshold). Downloadable as PNG.

### OMIM highlighting

Genes in the curated neurodevelopmental disease gene list (`sys_ndd_2025_subset.tsv`) are automatically flagged in annotation tables: dark blue for OMIM genes with a ClinVar pathogenic variant, light blue for OMIM genes without pathogenic annotation.

---

## Input files — what you need

### 1. How to specify target regions

You can define genomic regions to analyse in two ways:

**Option A — Gene list** (`.txt`): a plain text file with one official HGNC gene symbol per line. uncoverappLib translates symbols to genomic coordinates automatically using the UCSC transcript database.

```
BRCA1
BRCA2
TP53
POLG
```

Requirements: official HGNC symbols (e.g. `TP53` not `p53`), one per line, no header.

**Option B — Target BED** (`.bed`): a four-column tab-separated file with custom genomic regions. Use this for custom capture designs, amplicon panels, or specific exon subsets.

```
chr17    41196312    41277500    BRCA1
chr13    32889617    32973809    BRCA2
chr15    89859516    89876985    POLG
```

Columns: chromosome, start, end, region name. Both `chr`-prefixed (`chr1`) and numeric (`1`) chromosome notation are supported — specify which with the `chromosome_notation` parameter.

---

### 2. How to provide coverage data

Coverage data is always provided via a `.list` file — a plain text file with one absolute file path per line, one file per sample:

```
/absolute/path/to/sample1.bam
/absolute/path/to/sample2.bam
```

**Option A — BAM files**: raw aligned sequencing data. uncoverappLib computes per-base coverage using `Rsamtools::pileup()` with configurable minimum mapping quality (MAPQ) and minimum base quality. BAM files must be indexed — a `.bam.bai` index file must be present in the same directory.

**Option B — Pre-computed BED coverage files**: per-base or per-interval coverage already computed by bedtools, mosdepth, samtools depth, or GATK. This is 10–50× faster than BAM processing and is recommended for large cohorts.

The BED coverage format requires at minimum four tab-separated columns: chromosome, start, end, coverage depth.

```
# 0-based format (bedtools genomecov -bg, mosdepth)
chr17    41196311    41196312    45

# 1-based format (samtools depth)
chr17    41196312    41196312    45
```

You must specify the coordinate system (`"0-based"` or `"1-based"`) when using BED input — the app converts everything to 1-based internally.

| Tool that generated the BED | `input_coord_system` |
|-----------------------------|----------------------|
| `bedtools genomecov -bg` | `"0-based"` |
| `mosdepth` | `"0-based"` |
| `samtools depth` | `"1-based"` |

**How to generate BED coverage:**
```bash
# bedtools (0-based)
bedtools genomecov -ibam sample.bam -bg > sample_coverage.bed

# mosdepth (0-based, fast)
mosdepth --no-abbrev sample_prefix sample.bam

# samtools depth (1-based)
samtools depth -a sample.bam > sample_depth.bed
```

---

### 3. The output BED — bridge between standalone and app

The primary output of `buildInput()` is a tab-separated BED file with one row per genomic position and one coverage column per sample:

```
chromosome    start      end        SYMBOL    count_sample1    count_sample2
chr17         41196312   41196312   BRCA1     45               38
chr15         89859516   89859516   POLG      12               67
```

**This file can be loaded directly into the Coverage Analysis tab of the interactive app.** This is the recommended workflow for large datasets: run `buildInput()` once from the R console to generate the coverage matrix for all samples, then load the BED into the app for interactive exploration, visualisation, annotation, and export.

---

## Installation

### Step 1 — Install Bioconductor packages first

Several dependencies come from Bioconductor. Install them before installing uncoverappLib to ensure correct version resolution:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c(
  "Gviz",
  "Homo.sapiens",
  "org.Hs.eg.db",
  "TxDb.Hsapiens.UCSC.hg19.knownGene",
  "TxDb.Hsapiens.UCSC.hg38.knownGene",
  "EnsDb.Hsapiens.v75",
  "EnsDb.Hsapiens.v86",
  "OrganismDbi",
  "Rsamtools",
  "GenomicRanges",
  "AnnotationDbi",
  "GenomeInfoDb",
  "IRanges",
  "S4Vectors",
  "BiocStyle"
))
```

### Step 2 — Install uncoverappLib

```r
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")

devtools::install_github("annaballestrazzi/uncoverappLib")
```

All remaining CRAN dependencies (`shiny`, `shinyjs`, `shinyBS`, `shinyWidgets`, `shinycssloaders`, `DT`, `openxlsx`, `condformat`, `stringr`, `rappdirs`, `rlist`, `processx`, `dplyr`, `tidyselect`, `waiter`, `markdown`, `knitr`, `rmarkdown`, `testthat`) are installed automatically.

### Step 3 — Download annotation files (one-time, requires internet)

```r
library(uncoverappLib)

setup_uncoverapp()    # downloads dbNSFP/VEP annotation files — ~1 GB, 10–20 min
check_annotations()   # verifies files are present and valid
```

After this step everything works offline. Files are cached at `~/.cache/uncoverappLib/` (Linux/macOS) or `%APPDATA%/uncoverappLib/` (Windows).

**For production or air-gapped environments**, set environment variables to point to local annotation files and skip the download entirely:

```r
Sys.setenv(UNCOVERAPP_HG19_ANNOTATION = "/path/to/sorted_hg19.bed.gz")
Sys.setenv(UNCOVERAPP_HG38_ANNOTATION = "/path/to/sorted_hg38.bed.gz")
```

The corresponding `.tbi` tabix index must be in the same directory as the `.bed.gz` file.

---

## Quick start

### Interactive app

```r
library(uncoverappLib)
uncoverappLib.run()
```

**Typical workflow:**

1. **Processing tab** — upload gene list / target BED + BAM or BED list → set parameters → *Process Coverage Files*. Use the filter panel that appears to subset by sample and threshold, then *Download XLSX* (filtered or full).
2. **Coverage Analysis tab** — load the output BED → select genome, threshold, sample, and filter mode → *Calculate Low Coverage Regions* → *Calculate Annotations on Low Coverage* → *Download Annotations*
3. **Calculate AF tab** — set disease parameters to compute maxAF and filter rare variants
4. **Binomial distribution tab** — enter a genomic position to model detection probability

### Standalone pipeline

```r
library(uncoverappLib)

# Generate coverage matrix from BAM files
buildInput(
  geneList            = "genes.txt",     # gene list OR target BED
  sampleList          = "samples.list",  # paths to BAM or BED coverage files
  genome              = "hg38",
  chromosome_notation = "chr",
  type_input          = "genes",         # "genes" or "target"
  type_coverage       = "bam",           # "bam" or "bed"
  MAPQ.min            = 20,              # BAM only
  base.quality        = 20,             # BAM only
  outDir              = "./results"
)

# Annotate low-coverage positions for one sample
buildAnnotation(
  sample_data        = "./results/output/DATE.bed",
  target_sample      = "sample1",
  coverage_threshold = 20,
  genome             = "hg38",
  output_formatted   = "./results/sample1_annotated.xlsx"
)
```

For all scenarios combining gene list / target BED × BAM / BED coverage × single / multiple samples, with complete parameter explanations, see `inst/extdata/USER_GUIDE.md`.

---

## Outputs

### From the Processing tab / `buildInput()`

| File | Description |
|------|-------------|
| `DATE.bed` | Coverage matrix: chromosome, start, end, SYMBOL, `count_<sample>` per column. Loadable in the app. |
| `DATE_statistical_summary.txt` | Per-gene per-sample statistics: weighted mean/median coverage, bases below threshold, OMIM annotation. |
| `coverage_data_DATE.xlsx` | Optional XLSX from the filter panel — filtered subset (if *Apply Filter* was used) or full table. |

### From the Coverage Analysis tab / `buildAnnotation()`

| File | Description |
|------|-------------|
| `annotated_all_lowcov.xlsx` | Annotated variants with conditional colour formatting: ClinVar (red/green), CADD ≥20 (red), MutationAssessor H=red/M=yellow, M-CAP D (red), AF\_gnomAD <0.01 (red), OMIM genes (blue), high-impact variants (yellow row). |
| `annotated_all_lowcov.tsv` | Same data as tab-separated text (if `output_intersect` specified in `buildAnnotation()`). |

### Colour coding in the annotated Excel

| Colour | Meaning | Condition |
|--------|---------|-----------|
| 🔴 Red | Pathogenic signal | AF_gnomAD < 0.01, CADD ≥ 20, ClinVar pathogenic, M-CAP D |
| 🟢 Green | Benign signal | AF_gnomAD > 0.01, SIFT Tolerated |
| 🟡 Yellow | High-impact variant | MutationAssessor H/M + ClinVar + AF < 0.01 |
| 🔵 Blue (dark) | OMIM gene with ClinVar pathogenic variant | OMIM NDD gene list |
| 🔵 Blue (light) | OMIM gene without pathogenic annotation | OMIM NDD gene list |
---


## Documentation

| Document | Description |
|----------|-------------|
| `inst/extdata/QUICK_START.md` | Installation and app launch — get up and running in minutes |
| `inst/extdata/USER_GUIDE.md` | Complete documentation: every parameter, every tab, all scenarios for both interactive and standalone mode, troubleshooting, FAQ |
| `inst/extdata/prep_input_full.md` | Detailed guide to the Processing tab and `buildInput()` parameters |
| `vignettes/uncoverapp_introduction.Rmd` | Package introduction vignette |
| `vignettes/uncoverapp_technical.Rmd` | Technical reference vignette |
| `?buildInput`, `?buildAnnotation` | R help pages |

---

## Performance

| Scenario | BAM input | BED input |
|----------|-----------|-----------|
| 1 sample, 10 genes | 2–5 min | 10–30 sec |
| 1 sample, 100 genes | 10–20 min | 1–2 min |
| 10 samples, 50 genes | 1–2 hours | 5–10 min |

For large analyses, pre-compute coverage once with `bedtools genomecov -bg` or `mosdepth` and reuse the BED file across multiple sessions and gene panels.

---

## Troubleshooting

**"No annotation files found"** — run `setup_uncoverapp()` then `check_annotations()`.

**"Gene not found"** — use official HGNC symbols (`TP53` not `p53`). Check that the genome (hg19/hg38) matches your data.

**"No chromosome overlap"** — check notation: `chr1` → `chromosome_notation = "chr"`, `1` → `chromosome_notation = "number"`.

**"Coordinate mismatch / positions off by 1"** — set `input_coord_system = "0-based"` for bedtools/mosdepth, `"1-based"` for samtools depth.

**"App is slow"** — use BED coverage files instead of BAM (10–50× faster). Run `buildInput()` from the console and load the output BED into the app.

Full troubleshooting: `inst/extdata/USER_GUIDE.md`  
Report bugs: https://github.com/annaballestrazzi/uncoverappLib/issues

---

## Citation

```bibtex
@article{iovino2020uncoverapp,
  title   = {unCOVERApp: an interactive graphical application for clinical
             assessment of sequence coverage at the base-pair level},
  author  = {Iovino, Emanuela and Pippucci, Tommaso and Ballestrazzi, Anna},
  journal = {bioRxiv},
  year    = {2020},
  doi     = {10.1101/2020.02.10.939769},
  url     = {https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1}
}
```

---

## Project structure

```
uncoverappLib/
├── R/
│   ├── buildInput.R               Standalone coverage pipeline
│   ├── buildAnnotation.R          Standalone annotation pipeline
│   ├── getAnnotationFiles.R       Annotation file management
│   ├── run.uncoverapp.R           App launcher
│   └── setup.R                    First-time setup
├── inst/extdata/
│   ├── QUICK_START.md             Eight scenario commands
│   ├── USER_GUIDE.md              Complete user documentation
│   ├── prep_input_full.md         Processing tab detailed guide
│   ├── intro.md                   App home page
│   ├── sys_ndd_2025_subset.tsv    OMIM NDD gene database
│   ├── example_POLG.bam           Example BAM file
│   ├── example_bam_POLG.list      Example BAM list
│   ├── example_bed_POLG.list      Example BED list
│   ├── example_genes.txt          Example gene list
│   └── example.bed                Example coverage BED
├── shiny-dir/www/
│   ├── logo.png
│   └── modern-style.css
├── server.R / ui.R                Shiny app
├── compute-*.R                    Shiny compute modules (8 files)
├── waiter-helpers.R               Loading spinner
└── vignettes/
    ├── uncoverapp_introduction.Rmd
    └── uncoverapp_technical.Rmd
```

---

## License

MIT · See `LICENSE` file

## Contributors

**Emanuela Iovino** (emanuela.iovino@unibo.it) — creator and author  
**Tommaso Pippucci** (tommaso.pippucci@unibo.it) — author  
**Anna Ballestrazzi** (anna.ballestrazzi@studio.unibo.it) — author  
**Affiliation:** University of Bologna

**Acknowledgements:** dbNSFP · UCSC Genome Browser · gnomAD · ClinVar
