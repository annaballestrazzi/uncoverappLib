---

## 🚀 Quick Reference

**Most Common Use Cases:**

1. **First-time processing from BAM:**
```r
   buildInput(geneList = "genes.txt", sampleList = "bams.list", 
              type_coverage = "bam", MAPQ.min = 20, base.quality = 20)
```

2. **Fast processing with pre-computed coverage:**
```r
   buildInput(geneList = "genes.txt", sampleList = "coverage.list", 
              type_coverage = "bed", input_coord_system = "0-based")
```

3. **Using target BED regions:**
```r
   buildInput(geneList = "targets.bed", type_input = "target", ...)
```

---
## Output Format

The resulting uncoverappLib input file is a **tab-separated BED file** containing:

**Required columns:**
- `chromosome`: Chromosome identifier (chr1, chr2, ..., chrX, chrM or 1, 2, ..., X, M)
- `start`: Start position (1-based after processing)
- `end`: End position (1-based)
- `count_sample1`, `count_sample2`, ...: Coverage depth for each sample

**Example output:**
```
chromosome    start      end        count_patient_001    count_patient_002
chr17         41196312   41196312   45                   38
chr17         41196313   41196313   48                   42
chr17         41196314   41196314   52                   45
```

The file contains coverage information for all target genes/regions across all samples listed in the `.list` file.

---

## Processing Time

Processing time depends on:
- **File type**: BED coverage files are ~10-50x faster than BAM files
- **Number of genes**: More genes = longer processing
- **Number of samples**: Each sample is processed sequentially
- **File size**: Larger BAM/BED files take longer

**Estimated times (1 sample, 100 genes):**
- BAM file processing: 5-15 minutes
- BED coverage file: 30 seconds - 2 minutes

**For large analyses (>1000 genes) or cohorts (>20 samples):**
We recommend using the `buildInput()` function in R console before launching the app, then loading the output BED file in the **Coverage Analysis** page.

---

## Input Parameters

### Required Parameters

**Reference Genome:**
- Select `hg19` or `hg38` to match your data

**Chromosome Notation:**
- `chr`: Chromosomes named chr1, chr2, ..., chrX, chrM
- `number`: Chromosomes named 1, 2, ..., X, M

**Load input filtering file:**
Upload one of:
- **Gene list** (`.txt`): HGNC gene symbols, one per line
```
  BRCA1
  BRCA2
  TP53
```
- **Target BED** (`.bed`): Genomic regions with 4 columns
```
  chr17    41196312    41277500    BRCA1
  chr13    32889617    32973809    BRCA2
```

**Choose the type of your input file:**
- `List of gene names`: For gene list input
- `Target BED`: For BED regions input

**Load genomic data file:**
Upload `.list` file with absolute paths to coverage files (one per line):

For BAM files:
```
/data/patient_001.bam
/data/patient_002.bam
/data/patient_003.bam
```

For BED coverage files:
```
/data/patient_001_coverage.bed
/data/patient_002_coverage.bed
/data/patient_003_coverage.bed
```

**File format:**
- `BAM file`: Raw sequencing data (requires pileup computation)
- `BED coverage file`: Pre-computed coverage (faster processing)

---

## Coverage File Type: BAM vs BED

### BAM Files (`type_coverage = "bam"`)

**When to use:**
- You have raw sequencing data
- You want full control over quality filters
- First-time processing of samples

**Required parameters:**
- `MAPQ.min`: Minimum mapping quality (default: 1, recommended: 20-30)
- `base.quality`: Minimum base quality (default: 1, recommended: 20)

**Processing:**
- Computes per-base coverage using `samtools pileup`
- Applies quality filters
- Generates coverage BED output
- **Slower but more flexible**

**Requirements:**
- BAM files must be indexed (`.bai` files present)
- `samtools` must be installed and in PATH

---

### BED Coverage Files (`type_coverage = "bed"`)

**When to use:**
- You already have coverage computed
- You want faster processing
- Re-analyzing with different gene sets

**Required parameters:**
- `input_coord_system`: Coordinate system of your BED files
  - `"0-based"`: Standard BED format (most common)
  - `"1-based"`: IGV/VCF format

**Processing:**
- Reads pre-computed coverage directly
- Trims to target regions
- Merges multiple samples
- **Much faster**

**Supported tools:**
- `bedtools genomecov -bg` (0-based)
- `samtools depth` (1-based, add header)
- `GATK DepthOfCoverage` (with formatting)
- `mosdepth` (0-based)

**BED coverage format:**
```
chr1    100    101    25
chr1    101    102    30
chr1    102    103    18
```
*Columns: chromosome, start, end, coverage_depth*

---

## Quality Filters (BAM files only)

**Minimum Mapping Quality (MAPQ):**
- Default: 1 (accept all mapped reads)
- Recommended: 20-30 (high-quality alignments)
- Higher values = more stringent filtering
- Filters out multi-mapped and low-confidence alignments

**Minimum Base Quality:**
- Default: 1 (accept all bases)
- Recommended: 20 (Phred score)
- Higher values = more confident base calls
- Filters out low-quality sequencing errors

**These parameters are ignored when using BED coverage files.**

---

## Coordinate Systems

**Critical:** Specify the correct coordinate system for BED coverage input.

### 0-based (BED standard)
- **Intervals:** `[start, end)` half-open
- **Example:** `chr1  100  101` covers position 101 only
- **Used by:** 
  - BED files
  - bedtools genomecov
  - UCSC Genome Browser
  - BAM files (internal storage)

### 1-based (IGV, VCF standard)
- **Intervals:** `[start, end]` closed
- **Example:** `chr1  101  101` covers position 101 only
- **Used by:**
  - VCF files
  - IGV browser
  - SAM format
  - samtools depth output
  - GATK tools

**Important:**
- BAM files are handled automatically (always converted to 1-based for output)
- For BED coverage input, you **must** specify `input_coord_system`
- The app converts all data to 1-based internally for consistency

---

## Output Files

After running `buildInput()`, two files are generated:

### 1. Coverage BED File

**Location:** `output/output/DATE.bed`

**Format:**
```
chromosome    start      end        count_patient_001    count_patient_002
chr15         89859516   89859516   68                   70
chr15         89859517   89859517   70                   72
chr15         89859518   89859518   73                   75
```

**Column naming:**
- `count_FILENAME`: Coverage for each sample
- Sample names derived from input filenames (e.g., `patient_001.bam` → `count_patient_001`)

**Coordinates:**
- All positions are **1-based** (regardless of input format)
- One row per genomic position
- Suitable for direct loading into the interactive app

---

### 2. Statistical Summary

**Location:** `output/output/DATE_statistical_summary.txt`

**Format:**
```
SYMBOL  sample              Total_bases  Mean_coverage  Median_coverage  percentage_bases_under_20x
BRCA1   sample_patient_001  80000        45.3           42               5.2
BRCA2   sample_patient_001  84000        38.7           35               8.1
TP53    sample_patient_001  20000        52.1           50               2.3
```

**Metrics per gene/region:**
- `Total_bases`: Total bases in target regions
- `Mean_coverage`: **Weighted** mean (accounts for interval lengths)
- `Median_coverage`: Median across all positions
- `number_of_intervals_under_20x`: Count of low-coverage intervals
- `bases_under_20x`: Total bases below 20x threshold
- `percentage_bases_under_20x`: Percentage of bases below threshold

**Download:** Click `Statistical_Summary` button in the app.

**Optional annotations:** If you provide an `annotation_file` (e.g., OMIM), disease/inheritance information is merged with coverage statistics.

---

## Using Pre-computed BED Coverage Files

### Generating BED Coverage

**From BAM with bedtools:**
```bash
bedtools genomecov -ibam sample.bam -bg > sample_coverage.bed
# Output is 0-based
```

**From BAM with samtools:**
```bash
samtools depth -a sample.bam > sample_depth.txt
# Output is 1-based, 3 columns: chr, pos, depth
```

**From BAM with mosdepth:**
```bash
mosdepth -n -t 4 sample_prefix sample.bam
# Output is 0-based BED
```

### Command-Line Usage

**With BAM files:**
```r
library(uncoverappLib)

buildInput(
  geneList = "genes.txt",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",       # <- BAM processing
  MAPQ.min = 20,                # <- Quality filters
  base.quality = 20,
  outDir = "./results"
)
```

**With BED coverage files:**
```r
buildInput(
  geneList = "genes.txt",
  sampleList = "coverage_beds.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",            # <- BED coverage
  input_coord_system = "0-based",   # <- Specify coordinate system!
  outDir = "./results"
)
```

**Advantages of BED coverage:**
- Much faster processing (10-50x speedup)
- Reproducible (same input = same output)
- Can reuse coverage computed by other pipelines
- No need for `samtools` or BAM index files

---

## Troubleshooting

### "No valid genes found"
- Check gene names are **official HGNC symbols** (e.g., `TP53` not `p53`)
- Verify correct reference genome (hg19 vs hg38)
- Check for typos in gene names
- Log file: `output/output/preprocessing_log.txt` lists unrecognized genes

### "Cannot read BAM file"
- Ensure BAM files are **indexed** (`.bai` files in same directory)
- Verify absolute paths in `.list` file
- Check file permissions (`chmod +r file.bam`)
- Test: `samtools view file.bam | head`

### "No chromosome overlap"
- Check chromosome naming matches between BAM and gene annotations
- Set `chromosome_notation` parameter correctly:
  - BAM has `chr1` → use `chromosome_notation = "chr"`
  - BAM has `1` → use `chromosome_notation = "number"`
- Verify BAM has mapped reads: `samtools idxstats file.bam`

### "Processing too slow"
- **Use BED coverage files instead of BAM** (fastest solution)
- Increase quality filters (`MAPQ.min = 30`, `base.quality = 30`)
- Process genes in batches (split gene list)
- Use command-line `buildInput()` instead of app interface

### "Genes have NO coverage"
- Check that gene regions overlap with sequencing targets
- Verify BAM chromosomes match reference genome
- List of genes without coverage: `output/output/genes_without_coverage.txt`
- May indicate: wrong genome build, poor capture, or annotation issues

### "Wrong coordinate system"
- BED coverage positions off by 1? Check `input_coord_system` parameter
- Compare output positions with IGV browser
- When in doubt: use `input_coord_system = "0-based"` for BED files

---

## Next Steps

After generating the coverage BED file:

1. **Go to Coverage Analysis page**
2. **Load BED file** using `Select input file` button
3. **Configure parameters:**
   - Reference genome (hg19/hg38)
   - Coverage threshold (e.g., 20x)
   - Sample name (e.g., `count_patient_001`)
   - Filter mode (gene/chromosome/region/all)
4. **Calculate Low Coverage:** Click button to identify gaps
5. **Calculate Annotations:** Click button to retrieve variant annotations from dbNSFP
6. **Explore results in multiple tabs:**
   - **"bed file" tab**: View raw input data
   - **"Low-coverage positions" tab**: See all positions below threshold
   - **"UCSC gene" tab**: Verify gene coordinates (gene mode only)
   - **"Gene coverage" tab**: Generate Gviz plot (gene mode only)
   - **"Annotations on low-coverage positions" tab**: View variant annotations with pathogenicity predictions
7. **Use analysis tools:**
   - **maxAF calculator** (tab "Calculate AF by allele frequency app"): Works with any filter mode (gene, chromosome, region, or all chromosomes)
   - **Binomial probability calculator** (tab "Binomial distribution"): Calculate probability of detecting variants at low coverage positions
   - **Generate gene coverage plot**: Create visual representation of coverage along gene structure (gene mode only)
8. **Download results:** Export annotated Excel files with conditional formatting

For detailed examples, see:
```r
vignette("uncoverapp-tutorial", package = "uncoverappLib")
```

---

## 📊 Gene Coverage Visualization

### Generate Gene Coverage Plot

The **Gene Coverage Plot** button creates an interactive Gviz visualization showing coverage tracks for a specific gene.

**When available:**
- Only in **"Gene name"** filter mode
- After verifying gene coordinates with "Lookup UCSC Gene"
- After calculating low coverage regions
- Requires valid gene in selected genome

**What it shows:**
1. **Chromosome ideogram** (top): Gene location on chromosome
2. **Gene model track**: Exons, introns, transcript structure
3. **Coverage track (blue)**: High coverage regions (above threshold)
4. **Coverage track (red)**: Low coverage regions (below threshold)
5. **Genome axis**: Genomic coordinates

---

### Complete Workflow for Gene Coverage Plot

**🔴 IMPORTANT: Follow ALL steps in order**

#### Step 1: Configure Gene Filter
1. Go to **Coverage Analysis** page
2. Load your BED file using **"Select input file"** button
3. Click **"load input file"** button
4. Set **"Filter By"** = **"Gene name"**
5. Enter gene symbol in **"Gene name"** field (e.g., `BRCA1`, `POLG`, `TP53`)

#### Step 2: Verify Gene Coordinates
6. Click **"Lookup UCSC Gene"** button (🔍 search icon)
7. Switch to **"UCSC gene"** tab
8. **Verify** the gene exists and coordinates are correct
   - Check chromosome, start, end positions
   - Confirm correct genome (hg19/hg38)
   - If table is empty → gene not found, check spelling/genome

#### Step 3: Calculate Coverage
9. Return to main view (first tab)
10. Set your **"Coverage threshold"** (e.g., 20)
11. Enter **"Sample"** name exactly as in column header (e.g., `count_patient_001`)
12. Click **"Calculate Low Coverage Regions"** button (green)
13. Wait for processing → switches to **"Low-coverage positions"** tab

#### Step 4: Generate Plot
14. Navigate to **"Gene coverage"** tab
15. Click **"Generate Gene Coverage Plot"** button (blue)
16. **Wait 1-2 minutes** for plot generation
17. Gviz plot appears showing coverage along gene structure

---

### Example Workflow (POLG gene)
```
Step-by-step example with real data:

1. Load file: example_coverage.bed
2. Click "load input file"
3. Filter By = "Gene name"
4. Gene name = "POLG"
5. Reference Genome = "hg19"
6. Click "Lookup UCSC Gene" → Verify POLG found
7. Go to "UCSC gene" tab → See coordinates chr15:89859516-89876985
8. Return to first tab
9. Coverage threshold = "20"
10. Sample = "count_example_POLG.bam"
11. Click "Calculate Low Coverage Regions"
12. Wait , "Low-coverage positions" tab shows results
13. Go to "Gene coverage" tab
14. Click "Generate Gene Coverage Plot"
15. Wait , Plot appears with blue/red coverage tracks
```

---

### Interpreting the Plot

**Coverage Tracks:**
- **Blue regions**: Coverage **above** threshold (adequate coverage)
  - Good for variant calling
  - No validation needed
- **Red regions**: Coverage **below** threshold (potential gaps)
  - May miss variants
  - Consider Sanger sequencing validation
  - Check if these are critical exons

**Gene Structure:**
- **Thick boxes**: Exons (coding regions)
- **Thin lines**: Introns (non-coding)
- **Direction**: Arrows show strand orientation
- **Transcript**: Usually canonical/longest transcript shown

**Common Patterns:**
- **First exon low**: Common due to GC content or capture design
- **Last exon low**: May be 3' UTR (check if coding)
- **Random gaps**: Technical issues, poor capture, repeats
- **Entire gene low**: Wrong sample? Check sample name

**Use Cases:**
- Identify which specific exons have insufficient coverage
- Prioritize regions for Sanger sequencing validation
- Assess systematic coverage issues across gene
- Generate publication-quality figures
- Compare coverage patterns between samples

---

### Troubleshooting Plot Generation

**"No plot appears after clicking button"**
- Did you complete **ALL steps** including "Lookup UCSC Gene"?
- Check gene name spelling (case-sensitive for some systems)
- Verify gene exists in selected genome (hg19 vs hg38)
- Check "Low-coverage positions" tab has data
- Wait longer (large genes take 1-2 minutes)

**"Gene not found" error**
- Click "Lookup UCSC Gene" first to verify gene exists
- Check correct genome selected (hg19/hg38)
- Try gene aliases (e.g., "PARK2" vs "PRKN")
- Verify gene is in UCSC database for selected genome

**"Plot is empty or all red"**
- Check sample name is correct (must match column header exactly)
- Verify coverage threshold is reasonable (try threshold = 1)
- Check if gene has any coverage data in your file
- May indicate poor sequencing of this gene

**"Plot generation is very slow"**
- Large genes (>100kb) take longer
- First plot generation loads packages (slower)
- Subsequent plots are faster
- Consider checking smaller genes first

**"Cannot export plot"**
- Right-click on plot → "Save image as..."
- Save as PNG format
- If right-click doesn't work, use screenshot tool
- Plot size is fixed (800x400 pixels)

---

### Plot Limitations

**Cannot use plot with:**
- "Chromosome" filter mode (use gene mode only)
- "Region coordinates" filter mode (use gene mode only)
- "All chromosomes" mode (use gene mode only)
- Multiple genes simultaneously (one gene at a time)

**Limitations:**
- Only shows canonical/longest transcript
- Cannot zoom or pan (automatic gene boundaries)
- Cannot customize colors or tracks
- Large genes (>100kb) may be slow
- Requires sufficient coverage data in gene region

**Works best with:**
- Single gene analysis
- Genes < 100kb
- Adequate coverage depth (mean >10x)
- Validated gene name (via "Lookup UCSC Gene")

---

### Tips and Best Practices

**Before generating plot:**
1.  **Always** click "Lookup UCSC Gene" first
2.  Verify gene coordinates in "UCSC gene" tab
3.  Ensure "Low-coverage positions" calculation completed
4.  Double-check sample name matches exactly

**For best results:**
- Use official HGNC gene symbols (e.g., `BRCA1` not `Brca1`)
- Start with small genes to test workflow
- Generate plots after confirming gene has coverage
- Export plots immediately after generation (they don't save automatically)

**Export options:**
- Right-click plot → "Save image as..." (PNG recommended)
- Use for publications, reports, presentations
- Include in Sanger sequencing request forms
- Share with clinical team to discuss low-coverage regions

---

### Understanding the Tabs

After loading your coverage file, the interface shows multiple tabs with different information:

#### 1. **bed file**
- Shows the **raw input data** (first 10,000 rows displayed)
- Columns: chromosome, start, end, coverage for each sample
- Use to verify:
  - Data loaded correctly
  - Sample names are correct (important for plot generation!)
  - Coverage values look reasonable
  - Chromosome notation matches expectations

#### 2. **Low-coverage positions**
- Appears after clicking **"Calculate Low Coverage Regions"**
- Shows **all genomic positions** below your coverage threshold
- Filtered by your selection (gene/chromosome/region/all)
- **Empty table = no gaps found** (indicates good coverage!)
- Contains: chromosome, start, end, coverage value
- **Required** before generating gene plot

#### 3. **UCSC gene**
- Shows **gene coordinates** from UCSC database
- Appears after clicking **"Lookup UCSC Gene"** button
- **CRITICAL** : You must click "Lookup UCSC Gene" and verify results here before generating plot
- Contains: chromosome, txStart, txEnd, name2 (gene symbol), strand
- Use to verify:
  - Gene exists in selected genome (hg19/hg38)
  - Correct gene coordinates
  - Transcript information
  - Gene boundaries match expectations

**Example output:**
```
name          chrom   strand  txStart    txEnd      name2
uc010whl.2    chr15   +       89859516   89876985   POLG
```

#### 4. **Gene coverage** (Gene mode only)
- Contains the **"Generate Gene Coverage Plot"** button
- Only active after completing Steps 1-3 (see workflow above)
- Shows **Gviz visualization** after clicking button
- Displays:
  - Gene model with exons/introns
  - High coverage regions (blue)
  - Low coverage regions (red)
  - Genomic coordinates
- Can be exported as image (right-click)

#### 5. **Annotations on low-coverage positions**
- Appears after clicking **"Calculate Annotations on Low Coverage"**
- Shows **variant annotations** from dbNSFP database for low coverage positions
- Contains critical columns:
  - **ClinVar**: Pathogenicity classification
  - **CADD_PHED**: Deleteriousness score (>20 = likely deleterious)
  - **MutationAssessor**: Functional impact (H = High, M = Medium)
  - **AF_gnomAD**: Population allele frequency
  - **dbSNP**: Variant identifier
  - **GENENAME**: Gene symbol
  - **HGVSc/HGVSp**: Variant nomenclature
- **Color coding**:
  - Red background = potentially pathogenic
  - Green background = likely benign
  - Yellow highlight = important variants (high impact + pathogenic + rare)
- **Download**: Excel file with conditional formatting preserved

---

### Lookup UCSC Gene Button (REQUIRED for Plot)

The **"Lookup UCSC Gene"** button is **mandatory** before generating gene coverage plots.

**Why it's required:**
- Verifies gene exists in selected genome
- Retrieves exact gene coordinates
- Loads transcript information needed for plot
- Prevents errors during plot generation

**How to use:**
1. Select **"Filter By = Gene name"**
2. Enter gene symbol in **"Gene name"** field (e.g., "TP53")
3. Click **"Lookup UCSC Gene"** button (🔍 icon)
4. **Wait** for query to complete
5. Switch to **"UCSC gene"** tab
6. **Verify** results:
   - Gene found → Proceed with coverage calculation
   - Empty table → Gene not found, check spelling/genome

**Output columns:**
- **name**: Transcript ID (e.g., uc001abc.1)
- **chrom**: Chromosome (e.g., chr17)
- **strand**: Orientation (+ or -)
- **txStart**: Transcript start position
- **txEnd**: Transcript end position
- **cdsStart**: Coding sequence start
- **cdsEnd**: Coding sequence end
- **exonCount**: Number of exons
- **name2**: Gene symbol (e.g., TP53)

**What to check:**
-  **name2** matches your input gene
-  **chrom** is expected chromosome
-  **txStart/txEnd** are reasonable coordinates
-  At least one transcript returned
-  Multiple transcripts = normal (longest is used for plot)

**Troubleshooting:**
- **Empty table** → Gene not found
  - Check spelling: "TP53" not "p53", "BRCA1" not "Brca1"
  - Verify genome selection (hg19/hg38)
  - Try gene alias (e.g., "PARK2" vs "PRKN")
  - Gene may not exist in UCSC database for selected genome
- **Wrong gene** → Homonym exists
  - Check chromosome matches expectation
  - Verify coordinates are reasonable
- **Unexpected coordinates** → Wrong genome build
  - Switch between hg19/hg38 and retry

**Example (BRCA1 in hg19):**
```
Input: Gene = "BRCA1", Genome = "hg19"

Click "Lookup UCSC Gene" → "UCSC gene" tab shows:

name          chrom   strand  txStart    txEnd      name2
uc010ctr.3    chr17   -       41196312   41277500   BRCA1
uc002ict.4    chr17   -       41196312   41322420   BRCA1

Gene found on chr17
Multiple transcripts (normal)
Ready to calculate coverage
```

---

## Advanced Options

### Adding Custom Annotations

You can merge additional gene-level annotations (e.g., OMIM, clinical databases) with coverage statistics:
```r
buildInput(
  geneList = "genes.txt",
  sampleList = "samples.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",
  outDir = "./results",
  annotation_file = "omim_annotations.txt"  # <- Custom annotations
)
```

**annotation_file format (tab-separated):**
```
SYMBOL    Disease                  Inheritance    OMIM_ID
BRCA1     Breast cancer           AD             113705
BRCA2     Breast cancer           AD             600185
TP53      Li-Fraumeni syndrome    AD             191170
```

**Requirements:**
- Must contain `SYMBOL` column matching HGNC gene symbols
- Additional columns are merged with statistical summary
- Tab-separated (`.txt` or `.tsv`)

### Batch Processing Multiple Gene Panels

Process different gene panels in parallel:
```r
# Cardiology panel
buildInput(
  geneList = "cardio_genes.txt",
  sampleList = "patients.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",
  outDir = "./cardio_results"
)

# Oncology panel
buildInput(
  geneList = "onco_genes.txt",
  sampleList = "patients.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",
  outDir = "./onco_results"
)
```

### Processing Whole Exome/Genome Data

For WES/WGS, use pre-computed BED coverage to speed up processing:
```bash
# Step 1: Generate coverage with bedtools (do once)
for bam in *.bam; do
  bedtools genomecov -ibam $bam -bg > ${bam%.bam}_coverage.bed
done

# Step 2: Create sample list
ls *_coverage.bed | xargs -I {} readlink -f {} > coverage.list
```
```r
# Step 3: Process in R (very fast)
buildInput(
  geneList = "all_genes.txt",      # Large gene list
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = "./wes_results"
)
```

---

## See Also

- **Interactive Analysis:** See vignette `vignette("interactive-mode")`
- **Batch Processing:** See vignette `vignette("batch-processing")`
- **API Reference:** `?buildInput`, `?buildAnnotation`
- **GitHub Issues:** Report bugs at [github.com/annaballestrazzi/uncoverappLib/issues](https://github.com/annaballestrazzi/uncoverappLib/issues)