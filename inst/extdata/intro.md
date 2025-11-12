### An interactive graphical web application for clinical assessment of sequence coverage at base-pair level

This is a web application for clinical assessment of sequence coverage. 
unCOVERApp allows:

- to display interactive plots showing sequence gene coverage down to base-pair
resolution and functional/clinical annotations of sequence 
positions within coverage gaps (`Coverage Analysis` page).

- to calculate the [maximum credible population allele frequency](http://cardiodb.org/allelefrequencyapp/) (AF) to be applied as AF 
filtering threshold tailored to the model of the disease-of-interest 
instead of a general AF cut-off (e.g. 1% or 0.1%) 
(`Calculate AF by allele frequency app` page).
  - **Note:** maxAF filtering requires gene-specific analysis. Use "Filter by: Gene name" in Coverage Analysis.

- to calculate the 95% probability of the binomial distribution to observe at 
least N variant-supporting reads (N is the number of successes) based on a 
user-defined allele fraction that is expected for the variant 
(which is the probability of success). Especially useful to obtain the 
range of variant-supporting reads that is most likely to occur at a given 
depth of coverage (DoC) (which is the number of trials) for somatic variants
with low allele fraction
(`Binomial distribution` page). 

---

## Documentation 

All unCOVERApp functionalities are based on the availability of a bed file 
containing tab-separated specification of genomic coordinates (chromosome, 
start position, end position), coverage value and reference:alternate 
allele counts read for each position.

### Input Files

In the **Processing and Statistical Summary** page, users prepare the bed file providing 
input files containing a list of genes and a list of BAM files, respectively: 

#### Gene List File

A text file, with `.txt` extension, containing HGNC official gene name(s) one per 
row and to be uploaded to `Load input file` box. An example file is
included in extdata of uncoverappLib package. 

Users can also use a **target BED file** with 4 columns: chromosome, start, end and gene/amplicon name.

**Example gene list:**

```
DNAJC8
GNB1
PEX10
RPL22
```

#### BAM/BED File List

A text file, with `.list` extension, containing absolute paths to BAM or BED coverage files
(one per row) and to be uploaded to `Load genomic data file` box.

In the output file, sample 1, 2, 3... correspond to the samples in the list file in row 1, 2, 3...

**Example BAM list (Unix/macOS):**

```
/home/user/bam/sample1.bam
/home/user/bam/sample2.bam
/home/user/bam/sample3.bam
```

**Example BED coverage list:**

```
/home/user/coverage/sample1_coverage.bed
/home/user/coverage/sample2_coverage.bed
/home/user/coverage/sample3_coverage.bed
```

**Important:** BED coverage files should use **0-based coordinates** (BED standard).

While all inputs are loading, a progress bar appears during processing phase. 
unCOVERApp input file generation fails if incorrect gene names are specified. 
An unrecognized gene name(s) table is displayed if such a case occurs.

### Output Format

**Example output BED file from Preprocessing:**

```
chr15   89859516        89859516        68
chr15   89859517        89859517        70
chr15   89859518        89859518        73
chr15   89859519        89859519        73
chr15   89859520        89859520        74
chr15   89859521        89859521        75
```

### Processing Time

The preprocessing time depends on the size of the BAM file(s) and on the number 
of genes to investigate. 

**For large analyses (>50 genes)**, we recommend using the `buildInput()` function 
in R console before launching the app, or using the batch mode (see below).

### Alternative Tools

Alternatively, other tools can generate compatible input files:
- [bedtools](https://bedtools.readthedocs.io/en/latest/#)
- [samtools depth](http://www.htslib.org/doc/samtools-depth.html)
- [GATK DepthOfCoverage](https://gatk.broadinstitute.org/hc/en-us)

In this case, users can load the file directly on the **Coverage Analysis** page 
in `Select input file` box.

---

## Coverage Analysis

Once the bed file is ready, users can move to **Coverage Analysis** page and push 
`Load prepared input file` button.

### Input Parameters

To assess sequence coverage, the following **input parameters** must be 
specified in the sidebar:

- **Reference Genome**: hg19 or hg38
- **Coverage threshold**: Required minimum coverage depth (default: 20x)
- **Sample**: Sample name to analyze (matches column name in coverage file)
- **Coordinate system**: 0-based (BED standard) or 1-based (IGV, VCF)

### Filtering Options

**Filter by:**
- **Gene name**: Enter HGNC official gene symbol (e.g., BRCA1)
  - Push `Apply` button to load gene coordinates
  - **Required for maxAF analysis**
- **Chromosome**: Select specific chromosome (chr1-chr22, chrX, chrY)
- **All chromosomes**: Analyze entire genome
- **Region coordinates**: Specify genomic region (e.g., chr7:37453745-45627643)

### Outputs

unCOVERApp generates the following **outputs**:

- **bed file**: Unfiltered BED file
- **Low coverage positions**: Filtered dataset showing positions below threshold
- **UCSC gene**: Gene coordinates table
- **Gene coverage**: Interactive plot displaying:
  - Chromosome ideogram
  - Genomic location
  - Gene annotations from *Ensembl*
  - Transcripts annotation from *UCSC*
  - Related table shows the number of uncovered positions in each exon
- **Annotations on low-coverage positions**: dbNSFP annotation table
  - Click `Calculate Annotations` button to retrieve annotations
  - Click `download` button to save as formatted Excel file
  - Cells are colored according to pre-specified thresholds for:
    - Allele frequency (AF_gnomAD)
    - CADD score
    - M-CAP prediction
    - SIFT prediction
    - Polyphen2 prediction
    - ClinVar classification
    - OMIM ID, HGVSp and HGVSc notations

---

## Calculate Maximum Credible Allele Frequency

In the **Calculate AF by allele frequency app** page, users can set 
allele frequency cut-offs based on specific assumptions about the genetic 
architecture of the disease.

**Important:** This analysis requires **gene-specific filtering**. Use "Filter by: Gene name" 
in the Coverage Analysis page before accessing this tab.

Users may click on the `download` button and save the resulting table as formatted 
Excel file.

**Parameters:**
- **Prevalence**: Disease prevalence (1 in X people)
- **Inheritance**: Monoallelic or biallelic
- **Allelic heterogeneity**: Proportion of cases due to variants in this gene
- **Genetic heterogeneity**: Proportion of cases due to genetic factors
- **Penetrance**: Probability that carriers develop disease

---

## Binomial Distribution

The **Binomial distribution** page calculates the 95% binomial probability 
distribution of variant-supporting reads at a specific genomic position.

**Required inputs:**
- **Genomic position**: Position to analyze (from low-coverage positions)
- **Allele fraction**: Expected fraction of variant reads (probability of success)
- **Variant reads**: Minimum number of variant reads required to support variant calling

The tool provides:
- Probability distribution plot
- Cumulative distribution function
- 95% confidence interval
- Interpretation of variant detection probability

---

## Batch Mode (Command-Line)

For automated pipelines and processing multiple samples, uncoverappLib can be used 
from the command line with the following functions:

### Generate Coverage Data

```r
library(uncoverappLib)

buildInput(
  geneList = "genes.txt",           # Gene list or BED file
  bamList = "samples.list",          # List of BAM/BED files
  genome = "hg19",                   # hg19 or hg38
  type_bam = "chr",                  # Chromosome notation
  type_input = "genes",              # "genes" or "target"
  type_coverage = "bam",             # "bam" or "bed"
  outDir = "./output",
  MAPQ.min = 20,                     # Minimum mapping quality
  base.quality = 20,                 # Minimum base quality
  input_coord_system = "1-based"    # Coordinate system
)

# Output:
# - output/output/DATE.bed (coverage data)
# - output/output/DATE_statistical_summary.txt
```

### Annotate Low-Coverage Positions

```r
annotate_all_lowcov(
  sample_data = "coverage.bed",
  target_sample = "sample1",
  coverage_threshold = 20,
  genome = "hg19",
  output_intersect = "annotated.tsv",
  output_formatted = "annotated.xlsx"
)

# Output:
# - annotated.tsv (tab-separated data)
# - annotated.xlsx (formatted Excel with conditional coloring)
```

### Batch Processing Multiple Samples

```r
samples <- c("sample1", "sample2", "sample3")
coverage_file <- "output/output/DATE.bed"

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

### Documentation

For detailed batch processing examples and tutorials, see the package vignette:

```r
# View vignette
vignette("batch-processing", package = "uncoverappLib")

# Or after installation
browseVignettes("uncoverappLib")
```

---

## Links and Resources

- **GitHub Repository:** https://github.com/annaballestrazzi/uncoverappLib
- **Installation:** 
  ```r
  devtools::install_github("annaballestrazzi/uncoverappLib", 
                           build_vignettes = TRUE)
  ```
- **Issues/Support:** https://github.com/annaballestrazzi/uncoverappLib/issues
- **BioRxiv Paper:** https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1
- **Zenodo (Annotation Files):** https://zenodo.org/records/17524340

---

## Citation

If you use uncoverappLib in your research, please cite:

Iovino E, Pippucci T, et al. (2020). unCOVERApp: a web application for clinical 
assessment and annotation of coverage gaps in target genes. *bioRxiv*. 
doi: 10.1101/2020.02.10.939769
