# Quick Start Guide

**Choose your scenario, get the exact command.**

---

## 📊 Pick Your Input Type

| # | You Have | Jump To |
|---|----------|---------|
| 1 | Single BED coverage + Gene list | [Scenario 1](#scenario-1-single-bed-coverage--gene-list) |
| 2 | Multiple BED coverage + Gene list | [Scenario 2](#scenario-2-multiple-bed-coverage--gene-list) |
| 3 | Single BED coverage + Target BED | [Scenario 3](#scenario-3-single-bed-coverage--target-bed) |
| 4 | Multiple BED coverage + Target BED | [Scenario 4](#scenario-4-multiple-bed-coverage--target-bed) |
| 5 | Single BAM + Gene list | [Scenario 5](#scenario-5-single-bam--gene-list) |
| 6 | Multiple BAMs + Gene list | [Scenario 6](#scenario-6-multiple-bams--gene-list) |
| 7 | Single BAM + Target BED | [Scenario 7](#scenario-7-single-bam--target-bed) |
| 8 | Multiple BAMs + Target BED | [Scenario 8](#scenario-8-multiple-bams--target-bed) |

---

## Scenario 1: Single BED Coverage + Gene List

**Your files:**
```
sample1_coverage.bed    → Pre-computed coverage
genes.txt               → HGNC gene symbols
```

**Command:**
```r
library(uncoverappLib)

buildInput(
  geneList = "genes.txt",
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based",    # ⚠️ Specify coordinate system
  outDir = "./results"
)
```

**coverage.list:**
```
/full/path/to/sample1_coverage.bed
```

**genes.txt:**
```
BRCA1
BRCA2
TP53
```

⏱️ **Time:** ~30 seconds  
💡 **Best for:** Quick single-sample analysis

---

## Scenario 2: Multiple BED Coverage + Gene List

**Your files:**
```
patient1.bed    → Coverage for patient 1
patient2.bed    → Coverage for patient 2
patient3.bed    → Coverage for patient 3
genes.txt       → Gene panel
```

**Command:**
```r
buildInput(
  geneList = "genes.txt",
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = "./results"
)
```

**coverage.list:**
```
/full/path/to/patient1.bed
/full/path/to/patient2.bed
/full/path/to/patient3.bed
```

⏱️ **Time:** ~1-2 minutes  
💡 **Best for:** Cohort analysis, multiple patients

**Output columns:**
```
chromosome  start  end  count_patient1  count_patient2  count_patient3
```

---

## Scenario 3: Single BED Coverage + Target BED

**Your files:**
```
sample1_coverage.bed    → Pre-computed coverage
targets.bed             → Custom genomic regions
```

**Command:**
```r
buildInput(
  geneList = "targets.bed",          # ⚠️ Target BED, not gene list
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "target",             # ⚠️ "target" not "genes"
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = "./results"
)
```

**targets.bed:**
```
chr17	41196312	41277500	BRCA1
chr13	32889617	32973809	BRCA2
chr17	7571720	7590868	TP53
```

⏱️ **Time:** ~30 seconds  
💡 **Best for:** Custom regions, amplicon panels

---

## Scenario 4: Multiple BED Coverage + Target BED

**Your files:**
```
patient1.bed    → Coverage files
patient2.bed
targets.bed     → Custom regions
```

**Command:**
```r
buildInput(
  geneList = "targets.bed",
  sampleList = "coverage.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "target",
  type_coverage = "bed",
  input_coord_system = "0-based",
  outDir = "./results"
)
```

**coverage.list:**
```
/full/path/to/patient1.bed
/full/path/to/patient2.bed
```

⏱️ **Time:** ~1-2 minutes  
💡 **Best for:** Custom panel on multiple samples

---

## Scenario 5: Single BAM + Gene List

**Your files:**
```
sample1.bam         → Raw sequencing data
sample1.bam.bai     → BAM index (required!)
genes.txt           → Gene panel
```

**Command:**
```r
buildInput(
  geneList = "genes.txt",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",             # ⚠️ BAM processing
  MAPQ.min = 20,                     # ⚠️ Required for BAM
  base.quality = 20,                 # ⚠️ Required for BAM
  outDir = "./results"
)
```

**bams.list:**
```
/full/path/to/sample1.bam
```

⏱️ **Time:** ~5-15 minutes (depends on BAM size)  
💡 **Best for:** First-time processing from raw data

---

## Scenario 6: Multiple BAMs + Gene List

**Your files:**
```
patient1.bam        → BAM files
patient2.bam
patient3.bam
*.bam.bai           → Index files (required!)
genes.txt           → Gene panel
```

**Command:**
```r
buildInput(
  geneList = "genes.txt",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "genes",
  type_coverage = "bam",
  MAPQ.min = 20,
  base.quality = 20,
  outDir = "./results"
)
```

**bams.list:**
```
/full/path/to/patient1.bam
/full/path/to/patient2.bam
/full/path/to/patient3.bam
```

⏱️ **Time:** ~5-15 minutes per BAM (sequential)  
💡 **Best for:** Cohort from raw sequencing data

---

## Scenario 7: Single BAM + Target BED

**Your files:**
```
sample1.bam         → Raw data
sample1.bam.bai     → Index
targets.bed         → Custom regions
```

**Command:**
```r
buildInput(
  geneList = "targets.bed",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "target",             # ⚠️ Target mode
  type_coverage = "bam",
  MAPQ.min = 20,
  base.quality = 20,
  outDir = "./results"
)
```

**bams.list:**
```
/full/path/to/sample1.bam
```

⏱️ **Time:** ~5-15 minutes  
💡 **Best for:** Custom amplicon panel from BAM

---

## Scenario 8: Multiple BAMs + Target BED

**Your files:**
```
patient1.bam        → Multiple BAMs
patient2.bam
*.bam.bai           → Indices
targets.bed         → Custom regions
```

**Command:**
```r
buildInput(
  geneList = "targets.bed",
  sampleList = "bams.list",
  genome = "hg38",
  chromosome_notation = "chr",
  type_input = "target",
  type_coverage = "bam",
  MAPQ.min = 20,
  base.quality = 20,
  outDir = "./results"
)
```

**bams.list:**
```
/full/path/to/patient1.bam
/full/path/to/patient2.bam
```

⏱️ **Time:** ~5-15 minutes per BAM  
💡 **Best for:** Custom panel cohort from raw data

---

## 🎯 Decision Tree

```
┌─ Do you have BAM or BED coverage?
│
├─ BAM files
│  ├─ Gene list?
│  │  └─ type_input = "genes", type_coverage = "bam"
│  │     MAPQ.min + base.quality required
│  │
│  └─ Target BED?
│     └─ type_input = "target", type_coverage = "bam"
│        MAPQ.min + base.quality required
│
└─ BED coverage files
   ├─ Gene list?
   │  └─ type_input = "genes", type_coverage = "bed"
   │     input_coord_system required
   │
   └─ Target BED?
      └─ type_input = "target", type_coverage = "bed"
         input_coord_system required

Multiple samples? → List all paths in .list file (one per line)
```

---

## 📋 Parameter Cheat Sheet

### Input Type

| You Have | Set Parameter |
|----------|---------------|
| Gene list (genes.txt) | `type_input = "genes"` |
| Target BED (regions) | `type_input = "target"` |

### Coverage Type

| You Have | Set Parameters |
|----------|----------------|
| BAM files | `type_coverage = "bam"`<br>`MAPQ.min = 20`<br>`base.quality = 20` |
| BED coverage | `type_coverage = "bed"`<br>`input_coord_system = "0-based"` |

### Coordinate System (BED only)

| Your Tool | Set Value |
|-----------|-----------|
| bedtools genomecov | `"0-based"` |
| mosdepth | `"0-based"` |
| samtools depth | `"1-based"` |
| GATK DepthOfCoverage | Check output format |

---

## ⚡ Speed Comparison

| Input Type | Processing Time | Use When |
|------------|-----------------|----------|
| **BED coverage** | ⚡ Fast (30 sec - 2 min) | You have pre-computed coverage |
| **BAM files** | 🐢 Slow (5-15 min per sample) | First-time processing, need quality control |

**Pro tip:** Generate BED coverage once with `bedtools`, reuse for multiple analyses.

---

## 🔧 After Processing

Once `buildInput()` completes:

### Option 1: Interactive App
```r
uncoverappLib.run()
```
Then:
1. Load the output BED file
2. Select sample
3. Set coverage threshold
4. Filter by gene/chromosome/region
5. Calculate annotations
6. Export Excel files

### Option 2: Batch Annotation
```r
buildAnnotation(
  sample_data = "results/output/DATE.bed",
  target_sample = "patient1",
  coverage_threshold = 20,
  genome = "hg38",
  output_formatted = "patient1_annotated.xlsx"
)
```

---

## 💡 Pro Tips

**Speed up processing:**
- Use BED coverage files (10-50x faster than BAM)
- Increase quality filters (`MAPQ.min = 30`)
- Process genes in batches

**For large cohorts:**
- Pre-compute coverage with `bedtools genomecov -bg`
- Use `type_coverage = "bed"` for all samples
- Process samples in parallel (see USER_GUIDE.md)

**File organization:**
```
project/
├── data/
│   ├── genes.txt
│   ├── targets.bed
│   ├── coverage.list
│   └── coverage/
│       ├── patient1.bed
│       └── patient2.bed
├── results/
│   └── output/
│       └── DATE.bed
└── annotations/
    ├── patient1_annotated.xlsx
    └── patient2_annotated.xlsx
```

---

## 🐛 Common Issues

**"No valid genes found"**
- Use official HGNC symbols: `TP53` not `p53`
- Check genome version (hg19 vs hg38)

**"Cannot open file"**
- Use absolute paths in .list files
- Check file permissions

**"No chromosome overlap"**
- BAM has "chr1" → use `chromosome_notation = "chr"`
- BAM has "1" → use `chromosome_notation = "number"`

**"Processing too slow"**
- Switch to BED coverage files
- Increase quality filters
- Process in batches

---

**Next:** See [USER_GUIDE.md](USER_GUIDE.md) for complete documentation

**Questions?** Open an issue: https://github.com/annaballestrazzi/uncoverappLib/issues
