## Processing Time

The processing time depends on:
- Size of BAM/BED files
- Number of genes to investigate
- Number of samples

**For large analyses (>1000 genes) or cohorts (>20 samples)**, we recommend using the `buildInput()` function 
in R console before launching the App, then loading the output file in the 
**Coverage Analysis** page.

---

## Input Parameters

To prepare the input file, the following **parameters** must be specified in the sidebar:

### Required Parameters

- **Reference Genome**: Choose reference genome (hg19 or hg38)

- **Chromosome Notation**: Chromosome notation in BAM files
  - `number`: Chromosomes as 1, 2, ..., X, M
  - `chr`: Chromosomes as chr1, chr2, ..., chrX, chrM

- **Load input filtering file**: Upload file containing:
  - **Gene list**: `.txt` file with HGNC gene names (one per row)
  - **OR Target BED**: BED file with 4 columns (chr, start, end, gene/amplicon name)

- **Choose the type of your input file**: 
  - `List of gene names`: For gene list input
  - `Target BED`: For BED file input

- **Load genomic data file**: Upload `.list` file containing absolute paths to:
  - **BAM files** (one per row), OR
  - **BED coverage files** (one per row)

- **File format**: Select type of genomic data
  - `BAM file`: For BAM files (pileup will be computed)
  - `BED coverage file`: For pre-computed BED coverage files (faster)

### Quality Filters (for BAM files)

- **Minimum Mapping Quality (MAPQ)**: Default value 1
  - Reads with MAPQ below this threshold are excluded
  - Higher values (e.g., 20) = more stringent filtering

- **Minimum Base Quality (QUAL)**: Default value 1
  - Bases with quality below this threshold are excluded
  - Higher values (e.g., 20) = more stringent filtering

---

## Coordinate Systems

**Important:** Specify the coordinate system of your input files:

### 0-based (BED standard)
- Intervals are `[start, end)` (half-open)
- **Example:** `chr1  100  101` covers base 101 only
- **Used by:** BED files, BAM files (internally)
- **Select:** `input_coord_system = "0-based"` for BED coverage input

### 1-based (IGV, VCF standard)
- Intervals are `[start, end]` (closed)
- **Example:** `chr1  101  101` covers base 101 only
- **Used by:** VCF files, IGV browser, SAM format
- **Select:** `input_coord_system = "1-based"` for most visualizations

**The app handles BAM files automatically using 1-based pileup positions.**

For BED coverage files, ensure you specify the correct coordinate system when running 
`buildInput()` from command line.

---

## Output Files

### 1. Coverage BED File

Location: `output/output/DATE.bed`

**Format:**
```
chromosome    start    end    count_sample1    count_sample2    count_sample3
chr15         89859516 89859516 68            70               65
chr15         89859517 89859517 70            72               68
chr15         89859518 89859518 73            75               71
```

**Columns:**
- `chromosome`: Chromosome identifier (with 'chr' prefix)
- `start`: Start position (1-based after processing)
- `end`: End position
- `count_sampleX`: Coverage depth for each sample

### 2. Statistical Summary

Location: `output/output/DATE_statistical_summary.txt`

**Format:**
```
SYMBOL  sample      Total_bases  Mean_coverage  Median_coverage  percentage_bases_under_20x
BRCA1   sample_1    80000        45.3           42               5.2
BRCA2   sample_1    84000        38.7           35               8.1
TP53    sample_1    20000        52.1           50               2.3
```

**Columns:**
- `SYMBOL`: Gene symbol
- `sample`: Sample identifier
- `Total_bases`: Total number of bases in gene target regions
- `Mean_coverage`: Weighted mean coverage across all positions
- `Median_coverage`: Median coverage value
- `number_of_intervals_under_20x`: Count of intervals below 20x
- `bases_under_20x`: Total bases with coverage below 20x
- `percentage_bases_under_20x`: Percentage of bases below 20x threshold

**Download:** Click `Statistical_Summary` button to download the report.

The report provides coverage metrics per gene (`List of genes name` input) or 
per amplicon (`Target BED` input).

---

## Using Pre-computed Coverage Files

If you have already computed coverage using other tools, you can use BED coverage files:

### Supported Tools

- **bedtools genomecov**: `-bg` option produces compatible output
- **samtools depth**: Produces compatible 3-column output (add gene annotation)
- **GATK DepthOfCoverage**: Can produce compatible output with formatting

### BED Coverage Format

**Required format:**
```
chr1    100    101    25
chr1    101    102    30
chr1    102    103    18
```

**Columns:**
- Column 1: Chromosome (with or without 'chr' prefix)
- Column 2: Start position
- Column 3: End position
- Column 4: Coverage depth

**Important - Coordinate System:** 

Specify the coordinate system used by your BED files:

- **0-based (default)**: Standard BED format
  - Most common: `bedtools genomecov -bg`, `samtools depth`
  - Intervals: `[start, end)` half-open
  
- **1-based**: Some tools use 1-based output
  - Less common for BED files
  - Intervals: `[start, end]` closed

The app will convert coordinates internally to 1-based for visualization. 
Select the correct input system in the UI or when calling `buildInput()`.

When running `buildInput()` with BED coverage:
```r
buildInput(
  geneList = "genes.txt",
  bamList = "coverage_files.list",
  genome = "hg19",
  type_bam = "chr",
  type_input = "genes",
  type_coverage = "bed",              # ← BED coverage
  input_coord_system = "0-based",     # ← Important!
  outDir = "./output"
)
```

**Advantages of using BED coverage:**
- ✅ Much faster processing (no pileup computation needed)
- ✅ Can use coverage computed by other pipelines
- ✅ Reproducible if same BED files used

---

## Troubleshooting

### "No valid genes found"
- Check that gene names are **official HGNC symbols**
- Verify the correct reference genome is selected (hg19 vs hg38)
- Check for typos in gene names

### "Cannot read BAM file"
- Ensure BAM files are **indexed** (`.bai` files present)
- Verify absolute paths in `.list` file are correct
- Check file permissions

### "Processing too slow"
- Consider using **BED coverage files** instead of BAM
- Reduce number of genes per batch
- Use `buildInput()` function from R console
- Increase MAPQ and base quality filters to reduce positions

### "Unrecognized chromosome names"
- Ensure chromosome notation matches BAM file format
- Check `Chromosome Notation` parameter setting
- BAM files should have consistent naming (all 'chr' prefix or all numeric)

---

## Command-Line Alternative

For large-scale analysis or automation, use the `buildInput()` function directly:

```r
library(uncoverappLib)

# For BAM files
buildInput(
  geneList = "cancer_panel.txt",
  bamList = "patient_bams.list",
  genome = "hg19",
  type_bam = "chr",
  type_input = "genes",
  type_coverage = "bam",
  outDir = "./results",
  MAPQ.min = 20,
  base.quality = 20
)

# For BED coverage files (faster)
buildInput(
  geneList = "cancer_panel.txt",
  bamList = "coverage_beds.list",
  genome = "hg19",
  type_bam = "chr",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = "./results"
)
```

See the package vignette for detailed examples:
```r
vignette("batch-processing", package = "uncoverappLib")
```

---

## Next Steps

After generating the input file:

1. Go to **Coverage Analysis** page
2. Load the generated BED file using `Select input file`
3. Configure analysis parameters:
   - Reference genome
   - Coverage threshold
   - Sample name
   - Filtering options (gene/chromosome/region)
4. Click `Calculate Low Coverage Regions` to identify gaps
5. Click `Calculate Annotations` to retrieve variant annotations
6. Explore results and download annotated tables

For **maxAF analysis**, ensure you filter by a specific gene using "Filter by: Gene name".