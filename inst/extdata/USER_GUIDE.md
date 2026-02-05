# uncoverappLib User Guide

**Complete documentation for clinical coverage analysis.**

---

## Table of Contents

- [Installation](#installation)
- [Getting Started](#getting-started)
- [Interactive Mode](#interactive-mode)
- [Batch Mode](#batch-mode)
- [Input Files](#input-files)
- [Output Files](#output-files)
- [Advanced Features](#advanced-features)
- [Troubleshooting](#troubleshooting)
- [FAQ](#faq)

---

## Installation

### Requirements

**System:**
- R ≥ 4.0.0
- Java Runtime Environment (for Excel export)
- 4 GB RAM minimum (8 GB recommended for large cohorts)
- 2 GB free disk space (for annotation cache)

**Internet:**
- Required for initial setup only (annotation download)
- App runs fully offline after setup

### Install Package

```r
# Install devtools if needed
install.packages("devtools")

# Install uncoverappLib from GitHub
devtools::install_github("annaballestrazzi/uncoverappLib")
```

### Download Annotations

**First-time setup (requires internet):**

```r
library(uncoverappLib)

# Download dbNSFP annotation files (~1GB)
# This takes 10-20 minutes depending on connection
setup_uncoverapp()

# Verify files downloaded correctly
check_annotations()
```

**What gets downloaded:**
- `sorted_hg19.bed.gz` + index (hg19 genome annotations)
- `sorted_hg38.bed.gz` + index (hg38 genome annotations)
- Files cached in: `~/.cache/uncoverappLib/` (Linux/Mac) or `%APPDATA%/uncoverappLib/` (Windows)

**Offline operation:**
After setup completes, the app works without internet. Annotations are queried from local cache.

---

## Getting Started

### Launch Interactive App

```r
library(uncoverappLib)

# Start the app in your default browser
uncoverappLib.run()
```

The app opens in a new browser window. Everything runs locally on your computer—no data is uploaded to external servers.

### Quick Test

Try the included example data:

```r
# Get example files
gene_file <- system.file("Test/example_gene.txt", package = "uncoverappLib")
sample_list <- system.file("Test/example.list", package = "uncoverappLib")

# Process example
buildInput(
  geneList = gene_file,
  sampleList = sample_list,
  genome = "hg19",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = tempdir()
)

# Check output
list.files(file.path(tempdir(), "output"))
```

---

## Interactive Mode

### Workflow Overview

```
1. Load coverage file
   ↓
2. Set parameters (genome, threshold, sample)
   ↓
3. Filter data (gene/chromosome/region/all)
   ↓
4. Calculate low coverage regions
   ↓
5. Calculate annotations (optional)
   ↓
6. Export results
```

### Step-by-Step Guide

#### 1. Load Your Data

**Option A: Use pre-processed BED file**
- Click **"Select input file"** button
- Choose your coverage BED file
- Click **"load input file"**

**Option B: Generate from BAM/BED files**
- Go to **"Processing and Statistical Summary"** page
- Upload gene list + sample list
- Click **"Process Coverage Files"**
- Wait for processing to complete
- Go to **"Coverage Analysis"** page
- Load the generated BED file

#### 2. Configure Parameters

**Reference Genome:**
- Select `hg19` or `hg38` to match your data

**Coverage Threshold:**
- Enter minimum coverage depth (e.g., `20`)
- Positions below threshold = "low coverage"

**Sample:**
- Select sample name from dropdown
- Matches column name in coverage file (e.g., `count_patient1`)

#### 3. Choose Filter Mode

**Gene name** (most common):
1. Enter HGNC gene symbol (e.g., `BRCA1`, `TP53`)
2. Click **"Lookup UCSC Gene"** button (🔍)
3. Verify gene coordinates in "UCSC gene" tab
4. Return to main view

**Chromosome:**
- Select specific chromosome (chr1-chr22, chrX, chrY, chrM)

**Region coordinates:**
- Enter format: `chr17:41243452-41277500`
- Or without chr: `17:41243452-41277500`

**All chromosomes:**
- Analyzes entire file (can be slow for large datasets)

#### 4. Calculate Low Coverage

Click **"Calculate Low Coverage Regions"** button (green)

**What happens:**
- Filters positions by threshold
- Switches to "Low-coverage positions" tab
- Shows table of gaps

**Empty table?**
- Good news! No coverage gaps found
- Your data has adequate coverage

#### 5. Annotate Variants (Optional)

Click **"Calculate Annotations"** button (orange)

**What happens:**
- Queries dbNSFP database for variants at low-coverage positions
- Takes 1-5 minutes depending on number of positions
- Results appear in "Annotations on low-coverage positions" tab

**Annotation columns:**
- **ClinVar:** Pathogenicity classification
- **CADD_PHED:** Deleteriousness score (>20 = likely deleterious)
- **MutationAssessor:** Functional impact (H=High, M=Medium, L=Low)
- **M_CAP:** Pathogenicity prediction (D=Deleterious, T=Tolerated)
- **SIFT:** Protein function prediction
- **Polyphen2:** Damaging score
- **AF_gnomAD:** Population allele frequency
- **HGVSc/HGVSp:** Variant nomenclature

**Color coding:**
- 🔴 Red background = Potentially pathogenic
- 🟢 Green background = Likely benign
- 🟡 Yellow highlight = High-impact variants (rare + pathogenic + high impact)
- 🔵 Blue background = OMIM genes (disease-associated)

#### 6. Generate Gene Coverage Plot (Gene Mode Only)

**Requirements:**
- Filter mode = "Gene name"
- Gene coordinates verified (clicked "Lookup UCSC Gene")
- Low coverage calculated

**Steps:**
1. Go to **"Gene coverage"** tab
2. Click **"Generate Gene Coverage Plot"** button
3. Wait 1-2 minutes
4. Plot appears showing:
   - Chromosome ideogram
   - Gene structure (exons/introns)
   - Coverage tracks (red = low, blue = high)
   - Genomic coordinates

**Export plot:**
- Right-click on image → "Save image as..."
- PNG format recommended

#### 7. Export Results

**Download Annotations:**
- Click **"Download Annotations"** button
- Saves Excel file with:
  - Conditional formatting (colors preserved)
  - All annotation columns
  - Filterable/sortable data

**Download Statistical Summary:**
- Go to "Processing and Statistical Summary" page
- Click **"Statistical Summary"** button
- Saves tab-separated file with coverage metrics per gene

---

## Batch Mode

### When to Use Batch Mode

✅ **Use batch mode for:**
- Multiple samples (>5)
- Automated pipelines
- Reproducible workflows
- Large gene lists (>100 genes)
- Scheduled/repeated analyses

❌ **Use interactive mode for:**
- Single samples
- Exploratory analysis
- Visual inspection
- Teaching/training

### Workflow

```
BAM/BED files → buildInput() → Coverage BED → buildAnnotation() → Annotated Excel
```

### buildInput() - Generate Coverage Matrix

**Purpose:** Process BAM or BED coverage files for target genes/regions

**Basic usage:**
```r
buildInput(
  geneList = "genes.txt",           # Input: gene list or target BED
  sampleList = "samples.list",      # Input: paths to BAM/BED files
  genome = "hg38",                   # Reference: hg19 or hg38
  chromosome_notation = "chr",       # "chr" or "number"
  type_input = "genes",              # "genes" or "target"
  type_coverage = "bam",             # "bam" or "bed"
  outDir = "./results"               # Output directory
)
```

**Parameters:**

| Parameter | Required | Options | Description |
|-----------|----------|---------|-------------|
| `geneList` | Yes | file path | Gene list (.txt) or target BED |
| `sampleList` | Yes | file path | List of BAM/BED files (.list) |
| `genome` | Yes | `"hg19"`, `"hg38"` | Reference genome |
| `chromosome_notation` | Yes | `"chr"`, `"number"` | Chromosome naming |
| `type_input` | Yes | `"genes"`, `"target"` | Input type |
| `type_coverage` | Yes | `"bam"`, `"bed"` | Coverage file type |
| `outDir` | Yes | directory path | Where to save output |
| `MAPQ.min` | BAM only | integer (default: 20) | Minimum mapping quality |
| `base.quality` | BAM only | integer (default: 20) | Minimum base quality |
| `input_coord_system` | BED only | `"0-based"`, `"1-based"` | Coordinate system |

**Output files:**
```
outDir/
└── output/
    ├── DATE.bed                        # Coverage matrix
    └── DATE_statistical_summary.txt    # Coverage metrics per gene
```

**Coverage BED format:**
```
chromosome  start      end        count_sample1  count_sample2
chr17       41196312   41196312   45             38
chr17       41196313   41196313   48             42
```

### buildAnnotation() - Annotate Low Coverage

**Purpose:** Identify variants at low-coverage positions

**Basic usage:**
```r
buildAnnotation(
  sample_data = "results/output/Mon_Nov_25_2024.bed",
  target_sample = "patient1",
  coverage_threshold = 20,
  genome = "hg38",
  output_formatted = "patient1_annotated.xlsx"
)
```

**Parameters:**

| Parameter | Required | Description |
|-----------|----------|-------------|
| `sample_data` | Yes | Path to coverage BED from buildInput() |
| `target_sample` | Yes | Sample name (column in BED file) |
| `coverage_threshold` | Yes | Coverage cutoff (e.g., 20) |
| `genome` | Yes | Reference genome (hg19/hg38) |
| `output_intersect` | No | Tab-separated output file |
| `output_formatted` | No | Formatted Excel file |

**Output files:**
```
patient1_annotated.xlsx    # Excel with conditional formatting
patient1_lowcov.tsv        # Tab-separated (if specified)
```

### Batch Processing Multiple Samples

```r
# Step 1: Generate coverage for all samples
buildInput(
  geneList = "panel_genes.txt",
  sampleList = "cohort_samples.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",
  MAPQ.min = 20,
  base.quality = 20,
  outDir = "./cohort_results"
)

# Step 2: Get output file path
coverage_file <- list.files("cohort_results/output", 
                            pattern = "\\.bed$", 
                            full.names = TRUE)[1]

# Step 3: Annotate each sample
samples <- c("patient1", "patient2", "patient3", "patient4", "patient5")

for (sample in samples) {
  cat("Processing:", sample, "\n")
  
  buildAnnotation(
    sample_data = coverage_file,
    target_sample = sample,
    coverage_threshold = 20,
    genome = "hg38",
    output_formatted = paste0("annotations/", sample, "_annotated.xlsx")
  )
}

cat("Cohort processing complete!\n")
```

### Parallel Processing (Advanced)

```r
library(parallel)

# Detect cores
n_cores <- detectCores() - 1

# Process samples in parallel
mclapply(samples, function(sample) {
  buildAnnotation(
    sample_data = coverage_file,
    target_sample = sample,
    coverage_threshold = 20,
    genome = "hg38",
    output_formatted = paste0("annotations/", sample, "_annotated.xlsx")
  )
}, mc.cores = n_cores)
```

---

## Input Files

### Gene List

**Format:** Plain text, one HGNC gene symbol per line

**Example (genes.txt):**
```
BRCA1
BRCA2
TP53
PTEN
APC
MLH1
```

**Requirements:**
- Official HGNC symbols (e.g., `TP53` not `p53`)
- One symbol per line
- No header
- Plain text (.txt extension)

### Target BED

**Format:** Tab-separated, 4 columns

**Example (targets.bed):**
```
chr17	41196312	41277500	BRCA1
chr13	32889617	32973809	BRCA2
chr17	7571720	7590868	TP53
chr10	89622870	89731687	PTEN
```

**Columns:**
1. Chromosome (chr1, chr2, ..., chrX, chrM or 1, 2, ..., X, M)
2. Start position (0-based or 1-based, must be consistent)
3. End position
4. Region name (e.g., gene symbol, amplicon ID)

**Use when:**
- Custom capture designs
- Specific exon subsets
- Non-gene regions
- Amplicon panels

### Sample List

**Format:** Plain text, one absolute file path per line

**Example (bams.list):**
```
/data/sequencing/patient_001.bam
/data/sequencing/patient_002.bam
/data/sequencing/patient_003.bam
```

**Example (coverage.list):**
```
/data/coverage/patient_001_coverage.bed
/data/coverage/patient_002_coverage.bed
/data/coverage/patient_003_coverage.bed
```

**Requirements:**
- Absolute paths (not relative)
- One file per line
- No header
- .list extension (recommended)

### BAM Files

**Requirements:**
- Indexed (`.bam.bai` file must exist in same directory)
- Aligned to hg19 or hg38
- Chromosome naming consistent with reference

**Quality filters:**
- `MAPQ.min`: Minimum mapping quality (default: 20)
  - Higher = more stringent (30 recommended for clinical)
- `base.quality`: Minimum base quality (default: 20)
  - Phred score threshold

### BED Coverage Files

**Format:** Tab-separated, at least 4 columns

**Example:**
```
chr1	100	101	25
chr1	101	102	30
chr1	102	103	18
```

**Columns:**
1. Chromosome
2. Start position
3. End position
4. Coverage depth

**Coordinate systems:**

**0-based (BED standard):**
- Intervals: `[start, end)` half-open
- Example: `chr1 100 101` covers position 101 only
- Used by: bedtools, mosdepth, UCSC

**1-based (IGV/VCF standard):**
- Intervals: `[start, end]` closed
- Example: `chr1 101 101` covers position 101 only
- Used by: samtools depth, IGV, VCF

**How to generate:**

```bash
# From BAM with bedtools (0-based)
bedtools genomecov -ibam sample.bam -bg > sample_coverage.bed

# From BAM with samtools (1-based)
samtools depth -a sample.bam > sample_depth.txt

# From BAM with mosdepth (0-based, fast)
mosdepth -n -t 4 sample_prefix sample.bam
```

---

## Output Files

### Coverage BED

**Location:** `outDir/output/DATE.bed`

**Format:**
```
chromosome  start      end        count_sample1  count_sample2
chr17       41196312   41196312   45             38
chr17       41196313   41196313   48             42
chr17       41196314   41196314   52             45
```

**Features:**
- One row per genomic position
- 1-based coordinates (regardless of input format)
- Column naming: `count_FILENAME`
- Can be loaded directly into interactive app

**Use for:**
- Visual inspection in app
- Further filtering/annotation
- Reanalysis with different parameters

### Statistical Summary

**Location:** `outDir/output/DATE_statistical_summary.txt`

**Format:**
```
SYMBOL  sample          Total_bases  Mean_coverage  Median_coverage  percentage_bases_under_20x
BRCA1   sample_patient1  80000       45.3           42               5.2
BRCA2   sample_patient1  84000       38.7           35               8.1
TP53    sample_patient1  20000       52.1           50               2.3
```

**Columns:**
- `SYMBOL`: Gene name
- `sample`: Sample identifier
- `Total_bases`: Total bases in target regions
- `Mean_coverage`: Weighted mean coverage
- `Median_coverage`: Median across all positions
- `number_of_intervals_under_20x`: Count of low-coverage intervals
- `bases_under_20x`: Total bases below 20x
- `percentage_bases_under_20x`: Percentage below threshold

**Use for:**
- Quality control reports
- Coverage metrics tracking
- Identifying problematic genes

### Annotated Excel File

**From:** `buildAnnotation()` or app download

**Features:**
- Conditional formatting (colors preserved)
- All dbNSFP annotation columns
- Filterable/sortable in Excel
- Ready for clinical review

**Color coding:**
- Red cells: Pathogenic predictions, ClinVar entries, rare variants
- Green cells: Benign predictions, common variants
- Yellow rows: High-impact variants (start/end columns highlighted)
- Blue cells: OMIM disease genes

**Use for:**
- Clinical variant review
- Prioritization of Sanger sequencing
- Reporting to clinicians

---

## Advanced Features

### Maximum Credible Allele Frequency (maxAF)

**Purpose:** Calculate disease-specific AF threshold for variant filtering

**Access:** "Calculate AF by allele frequency app" tab in interactive mode

**Requirements:**
- Must use "Filter by: Gene name" mode
- Gene-specific analysis only

**Parameters:**
- **Prevalence:** Disease prevalence (1 in X people)
- **Inheritance:** Monoallelic (dominant) or biallelic (recessive)
- **Allelic heterogeneity:** Proportion of cases due to this gene
- **Genetic heterogeneity:** Proportion of genetic vs environmental cases
- **Penetrance:** Probability carriers develop disease

**Example:**
- Disease: Hypertrophic cardiomyopathy
- Prevalence: 1 in 500
- Gene: MYH7 (accounts for 30% of cases)
- Inheritance: Autosomal dominant
- Penetrance: 50%

→ maxAF ≈ 0.0003 (0.03%)

**Output:**
- Filters variants with AF < maxAF
- Downloadable Excel file
- Useful for rare disease genetics

**Reference:** http://cardiodb.org/allelefrequencyapp/

### Binomial Probability Calculator

**Purpose:** Calculate detection probability for low-frequency variants

**Access:** "Binomial distribution" tab in interactive mode

**Use cases:**
- Somatic variants with low allele fraction
- Mosaic variants
- Subclonal tumor variants

**Inputs:**
- **Genomic position:** Position from low-coverage table
- **Allele fraction:** Expected variant reads fraction (e.g., 0.05)
- **Variant reads:** Minimum reads required for calling

**Outputs:**
- Probability distribution plot
- 95% confidence interval
- Interpretation (threshold appropriate or not)

**Example:**
- Coverage: 100x
- Expected AF: 5% (somatic variant)
- Your threshold: 10 reads

→ P(X ≥ 10) = 5% (likely to miss this variant!)

### Gene Coverage Plots

**Purpose:** Publication-quality visualization of coverage along gene

**Access:** "Gene coverage" tab in interactive mode

**Requirements:**
- Filter by: Gene name
- Clicked "Lookup UCSC Gene"
- Low coverage calculated

**Features:**
- Chromosome ideogram
- Gene structure (exons/introns)
- Coverage histogram (blue = adequate, red = low)
- Genomic coordinates

**Export:**
- Right-click → "Save image as..."
- PNG format

**Use for:**
- Publications
- Clinical reports
- Sanger sequencing planning

---

## Troubleshooting

### Installation Issues

**"Cannot install devtools"**
```r
# Try from CRAN mirror
install.packages("devtools", repos = "https://cloud.r-project.org/")
```

**"Package dependencies failed"**
```r
# Install Bioconductor packages first
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install(c("GenomicRanges", "Rsamtools", "Gviz"))

# Then retry uncoverappLib install
devtools::install_github("annaballestrazzi/uncoverappLib")
```

**"Java not found"**
- Install Java Runtime Environment
- Windows: Download from https://www.java.com
- Mac: `brew install openjdk`
- Linux: `sudo apt install default-jre`

### Annotation Setup Issues

**"Cannot download annotation files"**
```r
# Check internet connection
# Try manual download from Zenodo
browseURL("https://zenodo.org/records/17524340")

# Place files in:
# Linux/Mac: ~/.cache/uncoverappLib/
# Windows: %APPDATA%/uncoverappLib/
```

**"Annotation files corrupted"**
```r
# Remove cache and re-download
unlink("~/.cache/uncoverappLib", recursive = TRUE)
setup_uncoverapp()
```

### Processing Issues

**"No valid genes found"**
- ✅ Use official HGNC symbols: `TP53` not `p53`
- ✅ Check genome version (hg19 vs hg38)
- ✅ Verify gene exists in UCSC database

**"Cannot read BAM file"**
- ✅ Ensure BAM is indexed (`.bai` file exists)
- ✅ Use absolute paths in `.list` file
- ✅ Check file permissions
- ✅ Test: `samtools view file.bam | head`

**"No chromosome overlap"**
- ✅ Check chromosome naming:
  - BAM has "chr1" → `chromosome_notation = "chr"`
  - BAM has "1" → `chromosome_notation = "number"`
- ✅ Verify BAM has mapped reads: `samtools idxstats file.bam`

**"Processing too slow"**
- ✅ Use BED coverage files (10-50x faster)
- ✅ Increase quality filters (`MAPQ.min = 30`)
- ✅ Process in batches
- ✅ Use command-line instead of app interface

**"Genes have NO coverage"**
- ✅ Check gene regions overlap with sequencing targets
- ✅ Verify BAM chromosomes match reference
- ✅ May indicate wrong genome build

### App Issues

**"App doesn't open"**
```r
# Try specifying browser
options(browser = "firefox")  # or "chrome", "safari"
uncoverappLib.run()
```

**"Gene not found"**
- ✅ Click "Lookup UCSC Gene" first
- ✅ Check genome selection (hg19 vs hg38)
- ✅ Try gene aliases (e.g., "PARK2" vs "PRKN")

**"No low coverage found"**
- ✅ Good news! Your data has adequate coverage
- ✅ Try lowering threshold to verify
- ✅ Check sample name is correct

**"Plot doesn't appear"**
- ✅ Complete ALL steps: Lookup Gene → Calculate Coverage → Generate Plot
- ✅ Wait 1-2 minutes (large genes take time)
- ✅ Only works in "Gene name" filter mode

**"maxAF shows no results"**
- ✅ MUST use "Filter by: Gene name"
- ✅ Click "Calculate Annotations" first
- ✅ Not compatible with "All chromosomes" mode

---

## FAQ

**Q: Can I use uncoverappLib offline?**  
A: Yes! After running `setup_uncoverapp()` once, the app works completely offline.

**Q: What genome builds are supported?**  
A: hg19 (GRCh37) and hg38 (GRCh38)

**Q: Can I analyze whole exome/genome data?**  
A: Yes, but use BED coverage files for speed. Processing WES BAMs directly is very slow.

**Q: How long does processing take?**  
A: BAM: 5-15 min per sample. BED: 30 sec - 2 min per sample.

**Q: Can I process multiple samples together?**  
A: Yes! List all files in `.list` file. They're processed sequentially, output in single BED.

**Q: What's the difference between buildInput() and the app?**  
A: `buildInput()` = command-line batch processing. App = interactive exploration.

**Q: Can I customize annotation columns?**  
A: Currently no. All dbNSFP columns are included automatically.

**Q: How do I cite uncoverappLib?**  
A: See README.md citation section.

**Q: Is my data uploaded anywhere?**  
A: No! Everything runs locally on your computer. No data leaves your machine.

**Q: Can I use custom annotation databases?**  
A: Currently only dbNSFP is supported. Custom databases not yet available.

**Q: What if my BAM uses "1" instead of "chr1"?**  
A: Set `chromosome_notation = "number"` in `buildInput()`.

**Q: Can I analyze mitochondrial genes?**  
A: Yes! chrM is supported (or M if using `chromosome_notation = "number"`).

**Q: How do I get help?**  
A: Open an issue on GitHub: https://github.com/annaballestrazzi/uncoverappLib/issues

---

**Version:** 1.10.1  
**Last Updated:** February 2025

**Made with ❤️ for the clinical genomics community**
