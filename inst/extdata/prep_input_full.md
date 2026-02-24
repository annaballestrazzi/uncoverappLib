## Processing and Statistical Summary

In this section you can prepare the BED coverage file required for all downstream analysis in uncoverappLib. Upload your target regions and sequencing data, set parameters, and click **Process Coverage Files**.

The output BED file can be loaded directly into the **Coverage Analysis** tab for interactive exploration, or used as input to `buildAnnotation()` for batch annotation.

---

## Required inputs

### 1. Input filtering file

Upload one of the following:

**Gene list** (`.txt`) — one official HGNC gene symbol per line, no header. The app retrieves genomic coordinates from the UCSC database automatically.
```
BRCA1
BRCA2
TP53
POLG
```

**Target BED** (`.bed`) — four tab-separated columns: chromosome, start, end, region name. Use this for custom capture designs or amplicon panels.
```
chr17    41196312    41277500    BRCA1
chr13    32889617    32973809    BRCA2
chr15    89859516    89876985    POLG
```

Select the corresponding option: **List of gene names** or **Target BED**.

---

### 2. Genomic data file

Upload a `.list` file — a plain text file with one **absolute path** per line, one file per sample.

For BAM files:
```
/absolute/path/to/patient_001.bam
/absolute/path/to/patient_002.bam
```

For BED coverage files:
```
/absolute/path/to/patient_001_coverage.bed
/absolute/path/to/patient_002_coverage.bed
```

---

### 3. Reference genome and chromosome notation

Select `hg19` or `hg38` to match your data. Select `chr` if your files use `chr1`-style chromosome names, `number` if they use `1`-style.

---

### 4. Coverage file type

**BAM file** — uncoverappLib computes per-base coverage using `Rsamtools::pileup()`. Additional parameters required:
- **Minimum mapping quality (MAPQ)** — default 1, recommended 20–30
- **Minimum base quality** — default 1, recommended 20

BAM files must be indexed: a `.bam.bai` file must exist in the same directory.

**BED coverage file** — pre-computed per-base or per-interval coverage from bedtools, mosdepth, or samtools. Much faster (10–50×). Additional parameter required:
- **Coordinate system** — `0-based` for bedtools/mosdepth output, `1-based` for samtools depth output

| Tool | Coordinate system |
|------|------------------|
| `bedtools genomecov -bg` | 0-based |
| `mosdepth` | 0-based |
| `samtools depth` | 1-based |

---

## Output files

After processing, two files are saved in `outDir/output/`:

### Coverage BED (`DATE.bed`)

Tab-separated, one row per genomic position, one `count_<sample>` column per sample. All coordinates are 1-based.

```
chromosome    start      end        SYMBOL    count_patient_001    count_patient_002
chr17         41196312   41196312   BRCA1     45                   38
chr15         89859516   89859516   POLG      12                   67
```

This file can be loaded directly into the **Coverage Analysis** tab.

### Statistical summary (`DATE_statistical_summary.txt`)

Tab-separated, one row per gene per sample. Includes width-weighted mean and median coverage, bases below threshold, and OMIM disease annotation for genes in the NDD gene database.

```
SYMBOL  sample              Mean_coverage  Median_coverage  bases_under_20x  percentage_under_20x
BRCA1   count_patient_001   45.3           42               1250             5.2
POLG    count_patient_001   12.1           11               4320             68.4
```

Click **Statistical Summary** to download this file.

---

## Optional: filter and download the coverage table

After processing completes, a filter panel appears above the coverage table:

- **Coverage threshold** — positions where all selected samples are above this value are hidden; values above the threshold are displayed as `>N`
- **Filter on samples** — select one or more samples to include
- **Apply Filter** — applies the threshold and sample filter to the table view
- **Download XLSX** — exports the data: if *Apply Filter* was clicked, exports only the filtered subset; otherwise exports the full coverage table

---

## Processing time

| Input type | Estimated time (1 sample, 100 genes) |
|------------|--------------------------------------|
| BAM file | 5–15 minutes |
| BED coverage file | 30 seconds – 2 minutes |

For large analyses (>1000 genes) or cohorts (>20 samples), use `buildInput()` from the R console and load the output BED into the Coverage Analysis tab.

---

## Troubleshooting

**"No valid genes found"** — use official HGNC symbols (`TP53` not `p53`). Check the log file at `outDir/output/preprocessing_log.txt` for unrecognised genes.

**"Cannot read BAM file"** — ensure `.bai` index is in the same directory. Use absolute paths in the `.list` file. Check permissions with `chmod +r file.bam`.

**"No chromosome overlap"** — chromosome naming mismatch. Set `chromosome_notation = "chr"` for `chr1`-style, `"number"` for `1`-style.

**"Positions off by 1"** — wrong coordinate system for BED input. Use `0-based` for bedtools/mosdepth, `1-based` for samtools depth.

**"Genes have no coverage"** — check that gene regions overlap your sequencing targets. Unmatched genes are listed in `outDir/output/genes_without_coverage.txt`.

**"Processing too slow"** — switch to BED coverage files (10–50× faster), or use `buildInput()` from the R console.
