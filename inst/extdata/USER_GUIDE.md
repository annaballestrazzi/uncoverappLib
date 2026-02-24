# uncoverappLib — User Guide

**Complete documentation for the interactive app and standalone functions.**

---

## Table of Contents

1. [Overview](#1-overview)
2. [Input files — formats and requirements](#2-input-files-formats-and-requirements)
3. [Interactive app — Processing and Statistical Summary tab](#3-interactive-app-processing-and-statistical-summary-tab)
4. [Interactive app — Coverage Analysis tab](#4-interactive-app-coverage-analysis-tab)
5. [Interactive app — Calculate AF tab](#5-interactive-app-calculate-af-tab)
6. [Interactive app — Binomial Distribution tab](#6-interactive-app-binomial-distribution-tab)
7. [Standalone — buildInput()](#7-standalone-buildinput)
8. [Standalone — buildAnnotation()](#8-standalone-buildannotation)
9. [Connecting standalone and app](#9-connecting-standalone-and-app)
10. [Output files — formats and content](#10-output-files-formats-and-content)
11. [Excel colour coding](#11-excel-colour-coding)
12. [Troubleshooting](#12-troubleshooting)
13. [FAQ](#13-faq)

---

## 1. Overview

uncoverappLib identifies genomic positions with insufficient sequencing coverage and annotates them with functional and clinical variant information. It can be used interactively through a Shiny application or programmatically through two standalone R functions.

**Interactive app** — load a coverage file, set filters, identify low-coverage positions, retrieve variant annotations, and export results to Excel. Everything runs locally; no data is transmitted to external servers.

**Standalone functions** — `buildInput()` generates a coverage matrix from BAM or BED files; `buildAnnotation()` annotates low-coverage positions for a specific sample. These functions are suitable for batch processing and pipeline integration.

The BED coverage file produced by `buildInput()` can be loaded directly into the app. This is the recommended workflow for large datasets: generate the coverage matrix once from the R console, then use the app for interactive exploration.

---

## 2. Input files: formats and requirements

### Gene list

Plain text file, one official HGNC gene symbol per line, no header.

```
BRCA1
BRCA2
TP53
POLG
MLH1
```

Use official HGNC symbols — `TP53` not `p53`, `BRCA1` not `Brca1`. The app and `buildInput()` translate gene symbols to genomic coordinates automatically using the UCSC transcript database.

### Target BED

Four-column tab-separated file specifying custom genomic regions. Use this instead of a gene list when you have a custom capture design, amplicon panel, or want to define specific exon subsets.

```
chr17    41196312    41277500    BRCA1
chr13    32889617    32973809    BRCA2
chr15    89859516    89876985    POLG
```

Columns: chromosome, start position, end position, region name (gene symbol or amplicon ID). Both `chr`-prefixed (`chr1`) and numeric (`1`) chromosome notation are supported — set the `chromosome_notation` parameter accordingly.

### Sample list

Plain text file with one absolute file path per line, one file per sample. Used for both BAM and BED coverage input.

```
/absolute/path/to/sample1.bam
/absolute/path/to/sample2.bam
```

Paths must be absolute (not relative). No header. Recommended extension: `.list`.

### BAM files

Raw aligned sequencing data. Requirements:
- Indexed: a `.bam.bai` file must exist in the same directory
- Aligned to hg19 or hg38
- Chromosome naming must be consistent with the reference genome

Quality filters applied during pileup (configurable): minimum mapping quality (MAPQ) and minimum base quality (Phred score). Default for both is 20.

### BED coverage files

Per-base or per-interval coverage pre-computed by an external tool. Minimum four tab-separated columns: chromosome, start, end, coverage depth.

```
chr17    41196311    41196312    45
chr17    41196312    41196313    48
```

You must specify the coordinate system of the BED file:

| Tool | Output format | `input_coord_system` |
|------|--------------|----------------------|
| `bedtools genomecov -bg` | 0-based | `"0-based"` |
| `mosdepth` | 0-based | `"0-based"` |
| `samtools depth` | 1-based | `"1-based"` |
| GATK DepthOfCoverage | check output | check output |

uncoverappLib converts all input to 1-based coordinates internally. Specifying the wrong system shifts all positions by 1.

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

## 3. Interactive app: Processing and Statistical Summary tab

This tab takes your raw data (BAM or BED coverage files + gene list or target BED) and produces the coverage matrix needed for all downstream analysis.

### Step-by-step

**1. Select reference genome**  
Choose `hg19` or `hg38`. Must match the genome your BAM files were aligned to.

**2. Select chromosome notation**  
Choose `chr` if your files use `chr1`-style naming, `number` if they use `1`-style.

**3. Upload the input filtering file**  
Upload your gene list (`.txt`, one HGNC symbol per line) or target BED (`.bed`, four columns). Select the corresponding radio button: *List of gene names* or *Target BED*.

**4. Select coverage file type**  
Choose *BAM file* or *BED coverage file*. This reveals additional parameters:
- For BAM: set minimum mapping quality (MAPQ) and minimum base quality
- For BED: set coordinate system (0-based or 1-based)

**5. Upload the sample list**  
Upload your `.list` file with absolute paths to BAM or BED coverage files.

**6. Click Process Coverage Files**  
A loading overlay appears during processing. When complete, the coverage table appears in the *input for uncoverapp* tab of the main panel, and the statistical summary download button becomes active.

**7. Download the statistical summary (optional)**  
Click *Statistical Summary* to download a tab-separated file with per-gene per-sample statistics: width-weighted mean and median coverage, total bases below 20x threshold, and OMIM disease annotation.

### Optional: filter and download the coverage table

After processing, a filter panel appears above the coverage table. This lets you explore and export a subset of the data before opening the full Coverage Analysis tab.

- **Coverage threshold** — select a value (1–100). Positions where all selected samples are above the threshold are removed from the view. Values above threshold are displayed as `>N` rather than the exact number, to focus attention on the positions that matter.
- **Filter on samples** — select one or more samples using the multi-select picker. Only the selected samples are shown.
- **Apply Filter** — click to apply the selected threshold and sample filters to the table.
- **Download XLSX** — exports the data to Excel. If *Apply Filter* was clicked, the exported file contains only the filtered rows. If the filter was never applied, the full coverage table is exported.

---

## 4. Interactive app: Coverage Analysis tab

This is the main analysis tab. Filters are set **before** calculating — only data matching your filter criteria is processed.

### Load your data

In the main panel (right side):
- Click **Select input file** and choose your coverage BED file
- Check or uncheck **Header** depending on whether your file has a header row (default: checked)
- Click **load input file**

The *bed file* sub-tab shows the loaded data. Verify that sample column names are correct — you will need them to select the sample.

### Configure parameters (sidebar)

**Reference Genome** — select `hg19` or `hg38` to match your coverage data.

**Coverage threshold** — the minimum acceptable coverage depth. Positions with depth strictly below this value are considered low-coverage.

**Sample** — select a sample from the dropdown. The dropdown is populated automatically from the column names of your loaded file.

**Filter by** — choose one of four modes. Only the relevant input controls appear:

| Mode | What it does | Additional input required |
|------|-------------|--------------------------|
| Gene name | Restricts analysis to the genomic interval of a single gene | Enter HGNC symbol + click *Lookup UCSC Gene* |
| Chromosome | Restricts to a single chromosome | Select from dropdown (chr1–chr22, chrX, chrY, chrM) |
| Region coordinates | Restricts to a custom interval | Enter coordinates as `chr:start-end` |
| All chromosomes | Analyses the entire loaded file | None |

### Gene name mode — required step before calculating

If you select *Gene name*:
1. Enter the HGNC gene symbol in the **Gene name** field
2. Click **Lookup UCSC Gene** — this retrieves the gene's genomic coordinates from the UCSC database and populates the *UCSC gene* sub-tab
3. Verify the gene was found correctly in the *UCSC gene* sub-tab before proceeding

This step is mandatory. Without it, the low-coverage calculation and gene plot cannot run correctly.

### Calculate low-coverage positions

Click **Calculate Low Coverage Regions**.

The app filters positions matching your selected mode and threshold. Results appear in the *Low-coverage positions* sub-tab. An empty table means no positions are below the threshold — which is a good result.

### Calculate annotations

Click **Calculate Annotations on Low Coverage**.

The app queries the local tabix-indexed annotation database for all known variants overlapping the low-coverage positions. This takes 1–5 minutes depending on the number of positions. Results appear in the *Annotations on low-coverage positions* sub-tab.

The annotation table includes the following columns:

| Column | Description |
|--------|-------------|
| seqnames, start, end | Genomic coordinates |
| coverage, counts | Coverage depth and nucleotide counts |
| REF, ALT | Reference and alternate allele |
| dbsnp | dbSNP identifier |
| GENENAME | Gene symbol (highlighted if OMIM) |
| PROTEIN\_ensembl | Ensembl protein ID |
| MutationAssessor | Functional impact: H (high), M (medium), L (low), N (neutral) |
| SIFT | Protein function prediction |
| Polyphen2 | Damaging score |
| M\_CAP | Pathogenicity prediction: D (deleterious), T (tolerated) |
| CADD\_PHED | CADD Phred deleteriousness score (>20 = likely deleterious) |
| AF\_gnomAD | Population allele frequency in gnomAD |
| ClinVar | ClinVar pathogenicity classification |
| HGVSc\_VEP, HGVSp\_VEP | HGVS variant nomenclature from VEP canonical transcripts |

### Download Annotations

The **Download Annotations** button appears only after annotation has completed. Click it to export the annotated table as an Excel file with conditional colour formatting (see Section 11: Excel colour coding).

### Sub-tabs in the main panel

| Sub-tab | Content |
|---------|---------|
| bed file | The full loaded coverage file |
| Low-coverage positions | Positions below threshold, filtered by your selected mode |
| UCSC gene | Gene coordinates from UCSC (gene mode only; populated after *Lookup UCSC Gene*) |
| Gene coverage | Gviz coverage plot (gene mode only; see below) |
| Annotations on low-coverage positions | Annotated variant table |

### Gene coverage plot (gene mode only)

After completing the gene lookup and low-coverage calculation:
1. Go to the **Gene coverage** sub-tab
2. Click **Generate Gene Coverage Plot**
3. Wait 1–2 minutes — large genes take longer
4. The plot shows: chromosome ideogram, genomic coordinates, coverage track (blue = above threshold, red = below), gene structure with exons and introns
5. Click **Download Plot** to save as PNG

The plot is only available in *Gene name* filter mode, after *Lookup UCSC Gene* has been clicked and *Calculate Low Coverage Regions* has been completed.

---

## 5. Interactive app: Calculate AF tab

Computes the maximum credible population allele frequency (maxAF) for a specific disease model, then filters the annotation table to show only variants with gnomAD AF below that threshold.

### Parameters

| Parameter | Description |
|-----------|-------------|
| Inheritance | Monoallelic (dominant) or biallelic (recessive) |
| Prevalence | Disease prevalence expressed as 1 in N people |
| Allelic heterogeneity | Proportion of cases attributable to this specific gene (0–1) |
| Genetic heterogeneity | Proportion of cases with a genetic vs non-genetic cause (0–1) |
| Penetrance | Probability that a carrier develops the disease (0–1) |

### Formulas

**Monoallelic:** `maxAF = (1/2) × (1/prevalence) × hetA × hetG × (1/penetrance)`

**Biallelic:** `maxAF = √(1/prevalence) × hetA × √hetG × (1/√penetrance)`

### Output

- The computed maxAF value is displayed
- The annotation table is filtered to show only variants with `AF_gnomAD < maxAF`
- Click **Download maxAF** to export as Excel with colour formatting

This tab requires annotation to have been calculated first in the Coverage Analysis tab.

---

## 6. Interactive app: Binomial Distribution tab

Models the probability of detecting a variant at a specific genomic position given the observed coverage and an expected allele fraction.

### Parameters

| Parameter | Description |
|-----------|-------------|
| Allele fraction (p) | Expected fraction of variant-supporting reads (e.g. 0.05 for a 5% somatic variant) |
| Variant reads | Minimum number of variant-supporting reads you require to call the variant |
| Genomic position | A position from the low-coverage or annotation table |

### How it works

The app looks up the position in the annotation data. If an exact match is found, that position's coverage is used. If not, the closest annotated position is used and the distance is reported. The coverage value feeds the binomial model.

### Output

- **Binomial distribution plot** — bar chart of probability for each possible number of supporting reads, with 95% confidence interval bounds marked
- **Cumulative distribution plot** — P(X ≤ k) curve with your variant reads threshold marked
- **Interpretation box** — states whether your threshold falls within the 95% CI (green = detectable with confidence, red = likely below detection threshold), including the position used, gene, variant, gnomAD AF, and the probability of observing at least N reads

---

## 7. Standalone: buildInput()

Generates a coverage matrix from BAM or BED coverage files for target genes or regions. Output is a BED file loadable into the app.

### All scenarios

#### Scenario A — Gene list + BAM files

```r
buildInput(
  geneList            = "genes.txt",
  sampleList          = "samples.list",
  genome              = "hg38",
  chromosome_notation = "chr",
  type_input          = "genes",
  type_coverage       = "bam",
  MAPQ.min            = 20,
  base.quality        = 20,
  outDir              = "./results"
)
```

#### Scenario B — Gene list + BED coverage files

```r
buildInput(
  geneList            = "genes.txt",
  sampleList          = "coverage.list",
  genome              = "hg38",
  chromosome_notation = "chr",
  type_input          = "genes",
  type_coverage       = "bed",
  input_coord_system  = "0-based",    # or "1-based" depending on tool
  outDir              = "./results"
)
```

#### Scenario C — Target BED + BAM files

```r
buildInput(
  geneList            = "targets.bed",
  sampleList          = "samples.list",
  genome              = "hg38",
  chromosome_notation = "chr",
  type_input          = "target",     # ← "target" not "genes"
  type_coverage       = "bam",
  MAPQ.min            = 20,
  base.quality        = 20,
  outDir              = "./results"
)
```

#### Scenario D — Target BED + BED coverage files

```r
buildInput(
  geneList            = "targets.bed",
  sampleList          = "coverage.list",
  genome              = "hg38",
  chromosome_notation = "chr",
  type_input          = "target",
  type_coverage       = "bed",
  input_coord_system  = "0-based",
  outDir              = "./results"
)
```

### Parameter reference

| Parameter | Required | Values | Description |
|-----------|----------|--------|-------------|
| `geneList` | Yes | file path | Gene list (.txt) or target BED (.bed) |
| `sampleList` | Yes | file path | .list file with absolute paths to BAM or BED files |
| `genome` | Yes | `"hg19"`, `"hg38"` | Reference genome |
| `chromosome_notation` | Yes | `"chr"`, `"number"` | Chromosome naming convention in your files |
| `type_input` | Yes | `"genes"`, `"target"` | Type of region specification |
| `type_coverage` | Yes | `"bam"`, `"bed"` | Type of coverage input |
| `outDir` | Yes | directory path | Output directory (created if absent) |
| `MAPQ.min` | BAM only | integer, default 20 | Minimum mapping quality |
| `base.quality` | BAM only | integer, default 20 | Minimum base quality (Phred score) |
| `input_coord_system` | BED only | `"0-based"`, `"1-based"` | Coordinate system of input BED coverage files |
| `coverage_threshold` | No | integer | If set, exports an additional XLSX of positions below this threshold |

### Output

Two files are always generated in `outDir/output/`:
- `DATE.bed` — coverage matrix (chromosome, start, end, SYMBOL, `count_<sample>` columns)
- `DATE_statistical_summary.txt` — per-gene per-sample statistics

If `coverage_threshold` is specified, an additional file is generated:
- `DATE_low_coverage_<N>x.xlsx` — positions below the threshold with colour formatting

---

## 8. Standalone: buildAnnotation()

Annotates low-coverage positions for a single sample from a coverage BED file. Can be run on any BED file produced by `buildInput()` or loaded from the app.

### Basic usage

```r
buildAnnotation(
  sample_data        = "./results/output/DATE.bed",
  target_sample      = "sample1",
  coverage_threshold = 20,
  genome             = "hg38",
  output_formatted   = "./sample1_annotated.xlsx"
)
```

### Parameter reference

| Parameter | Required | Description |
|-----------|----------|-------------|
| `sample_data` | Yes | Path to coverage BED file (from buildInput or app) |
| `target_sample` | Yes | Sample name — must match the column name in the BED file |
| `coverage_threshold` | Yes | Coverage cutoff; positions below this value are annotated |
| `genome` | Yes | Reference genome: `"hg19"` or `"hg38"` |
| `output_intersect` | No | Path for tab-separated output file |
| `output_formatted` | No | Path for formatted Excel output with colour coding |

The `target_sample` value must match the sample column in the BED file. If the BED column is `count_patient1.bam`, use `target_sample = "patient1.bam"`.

### Batch annotation — multiple samples

After a single `buildInput()` run, annotate each sample in a loop:

```r
coverage_bed <- list.files("./results/output", pattern = "\\.bed$", full.names = TRUE)[1]

samples <- c("patient1", "patient2", "patient3")

for (s in samples) {
  message("Annotating: ", s)
  buildAnnotation(
    sample_data        = coverage_bed,
    target_sample      = s,
    coverage_threshold = 20,
    genome             = "hg38",
    output_formatted   = paste0("./annotations/", s, "_annotated.xlsx")
  )
}
```

---

## 9. Connecting standalone and app

The two modes are designed to complement each other. Common combined workflows:

### Generate coverage once, explore interactively

```r
# R console — generate coverage matrix for all samples
buildInput(
  geneList = "genes.txt", sampleList = "all_samples.list",
  genome = "hg38", chromosome_notation = "chr",
  type_input = "genes", type_coverage = "bam",
  MAPQ.min = 20, base.quality = 20, outDir = "./results"
)
```

Then in the app: open Coverage Analysis → load `results/output/DATE.bed` → explore interactively per gene and per sample.

### Pre-compute for a large cohort, annotate programmatically

```bash
# Shell — generate per-sample BED coverage files (fast)
for bam in /data/*.bam; do
  sample=$(basename $bam .bam)
  bedtools genomecov -ibam $bam -bg > /coverage/${sample}.bed
done
ls /coverage/*.bed | xargs realpath > coverage.list
```

```r
# R — build matrix and annotate all samples
buildInput(
  geneList = "panel.txt", sampleList = "coverage.list",
  genome = "hg38", chromosome_notation = "chr",
  type_input = "genes", type_coverage = "bed",
  input_coord_system = "0-based", outDir = "./cohort"
)

coverage_bed <- list.files("./cohort/output", pattern = "\\.bed$", full.names = TRUE)[1]
samples <- gsub("count_", "", grep("count_", colnames(read.table(coverage_bed, header=TRUE, nrow=1)), value=TRUE))

for (s in samples) {
  buildAnnotation(
    sample_data = coverage_bed, target_sample = s,
    coverage_threshold = 20, genome = "hg38",
    output_formatted = paste0("./annotations/", s, "_annotated.xlsx")
  )
}
```

---

## 10. Output files: formats and content

### Coverage BED (`DATE.bed`)

Tab-separated. One row per genomic position, one `count_<sample>` column per sample. Coordinates are 1-based regardless of input format.

```
chromosome    start      end        SYMBOL    count_sample1    count_sample2
chr17         41196312   41196312   BRCA1     45               38
chr17         41196313   41196313   BRCA1     48               42
chr15         89859516   89859516   POLG      12               67
```

### Statistical summary (`DATE_statistical_summary.txt`)

Tab-separated. One row per gene per sample.

| Column | Description |
|--------|-------------|
| SYMBOL | Gene name |
| sample | Sample identifier |
| Total\_bases | Total bases in target regions |
| Mean\_coverage | Width-weighted mean coverage |
| Median\_coverage | Width-weighted median coverage |
| bases\_under\_threshold | Total bases below the threshold |
| percentage\_under\_threshold | Percentage of bases below threshold |
| OMIM | OMIM disease annotation if gene is in NDD database |

### Annotated Excel (`*_annotated.xlsx`)

Excel workbook with conditional colour formatting. Contains all annotation columns listed in Section 4. See Section 11 for full colour coding rules.

### Low-coverage XLSX from Processing tab filter (`coverage_data_DATE.xlsx`)

Exported from the optional filter panel in the Processing tab. Contains the coverage matrix: if *Apply Filter* was clicked before downloading, the file contains only rows where at least one selected sample is below the threshold, with above-threshold values shown as `>N`. If no filter was applied, the full coverage table is exported.

---

## 11. Excel colour coding

### Annotation tables (from Coverage Analysis tab and buildAnnotation())

| Column | Condition | Colour |
|--------|-----------|--------|
| ClinVar | Value ≠ "." (variant present) | Red background |
| ClinVar | Value = "." (no entry) | Green background |
| CADD\_PHED | ≥ 20 | Red background |
| CADD\_PHED | < 20 | Green background |
| MutationAssessor | H (high impact) | Red background |
| MutationAssessor | M (medium impact) | Yellow background |
| MutationAssessor | other | Green background |
| M\_CAP | D (deleterious) | Red background |
| M\_CAP | other | Green background |
| AF\_gnomAD | < 0.01 | Red background |
| AF\_gnomAD | ≥ 0.01 | Green background |
| GENENAME | OMIM gene + pathogenic ClinVar | Dark blue background, white bold text |
| GENENAME | OMIM gene only | Light blue background |
| start, end | MutationAssessor H or M **AND** ClinVar ≠ "." **AND** AF < 0.01 | Yellow row highlight |

### maxAF table (from Calculate AF tab download)

Same rules as above, except AF\_gnomAD threshold is `< 0.05` instead of `< 0.01`.

---

## 12. Troubleshooting

**"No annotation files found"**  
Run `setup_uncoverapp()` then `check_annotations()`. If files are already downloaded but not found, set environment variables pointing to their location:
```r
Sys.setenv(UNCOVERAPP_HG19_ANNOTATION = "/path/to/sorted_hg19.bed.gz")
Sys.setenv(UNCOVERAPP_HG38_ANNOTATION = "/path/to/sorted_hg38.bed.gz")
```

**"Gene not found" / empty UCSC gene tab**  
Use official HGNC symbols (`TP53` not `p53`, `PRKN` not `PARK2`). Verify that the selected genome (hg19/hg38) matches your data. Some genes have different coordinates in hg19 vs hg38.

**"No chromosome overlap"**  
Check chromosome naming in your input files. If chromosomes are named `chr1`, set `chromosome_notation = "chr"`. If named `1`, set `chromosome_notation = "number"`. Use `samtools idxstats file.bam` to inspect BAM chromosome names.

**"Positions off by 1"**  
Wrong coordinate system for BED coverage input. Use `input_coord_system = "0-based"` for bedtools/mosdepth output, `"1-based"` for samtools depth.

**"Cannot open BAM file"**  
Ensure the `.bam.bai` index file exists in the same directory. Use absolute paths in the `.list` file. Test with `samtools view file.bam | head`.

**"Download Annotations button not appearing"**  
The button is shown only after *Calculate Annotations on Low Coverage* has completed successfully. If the annotation step produced 0 results, the button does not appear — check that low-coverage positions were found first.

**"App is slow / pileup takes too long"**  
Use pre-computed BED coverage files instead of BAM (10–50× faster). For large gene panels, run `buildInput()` from the R console and load the output BED into the app.

**"Plot does not appear"**  
In gene mode only: ensure you clicked *Lookup UCSC Gene* and verified the result in the *UCSC gene* sub-tab **before** clicking *Calculate Low Coverage Regions*. All three steps must be completed in order.

**"Binomial tab shows wrong position"**  
If the exact position you entered is not in the annotation data, the closest annotated position is used. The distance is reported in the output. Positions must come from the annotated low-coverage table.

---

## 13. FAQ

**Can I use uncoverappLib offline?**  
Yes. After running `setup_uncoverapp()` once, the app and all functions work without internet.

**Which genome builds are supported?**  
hg19 (GRCh37) and hg38 (GRCh38).

**Can I analyse whole-exome or whole-genome data?**  
Yes. Use pre-computed BED coverage files (from bedtools or mosdepth) for speed. Processing whole-exome BAM files directly via pileup is very slow.

**What is the difference between the gene list and the target BED input?**  
The gene list uses HGNC symbols and retrieves coordinates from the UCSC database — useful when you want to cover entire gene bodies including UTRs. The target BED lets you define exact intervals — useful for amplicon panels or when you want only specific exons.

**Can I process multiple samples together?**  
Yes. List all sample files in your `.list` file. `buildInput()` processes them and produces a single BED with one column per sample. You can then annotate each sample independently with `buildAnnotation()`.

**The output BED has column names like `count_sample1.bam` — is that correct?**  
Yes. Sample names are derived from the input file names. When running `buildAnnotation()`, use the full name as it appears in the column header (e.g. `target_sample = "sample1.bam"`).

**Can I load any BED coverage file into the app, or only files from buildInput()?**  
You can load any BED file that has chromosome, start, end, and one or more numeric coverage columns. The sample dropdown in Coverage Analysis is populated automatically from the column names.

**Does my data get uploaded anywhere?**  
No. Everything runs locally. No data leaves your machine.

**How do I cite uncoverappLib?**  
Iovino E, Pippucci T, Ballestrazzi A (2020). unCOVERApp: an interactive graphical application for clinical assessment of sequence coverage at the base-pair level. *bioRxiv*. doi: 10.1101/2020.02.10.939769

**Where do I report bugs?**  
https://github.com/annaballestrazzi/uncoverappLib/issues
