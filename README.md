# uncoverappLib

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17524340.svg)](https://doi.org/10.5281/zenodo.17524340)

**Clinical-grade coverage analysis. Interactive. Offline-capable. Ready for production.**

Identify coverage gaps and annotate variants at base-pair resolution. Built for clinical diagnostics, optimized for research workflows.

[Read the paper](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1) • [Quick Start](#-quick-start) • [Documentation](#-documentation)

---

## ✨ What You Get

**📊 Interactive Coverage Analysis**  
Visual inspection of sequencing coverage gaps across genes, chromosomes, or custom regions. Supports both BAM files and pre-computed BED coverage.

**🧬 Clinical Annotations**  
Automatic integration with dbNSFP v4.0 for pathogenicity predictions (CADD, SIFT, MutationAssessor), population frequencies (gnomAD), and ClinVar classifications.

**📈 Statistical Tools**  
Built-in calculators for maximum credible allele frequency (rare disease genetics) and binomial probability (somatic variant detection).

**💾 Offline-First Design**  
After one-time setup, runs completely without internet. Perfect for clinical labs with network restrictions.

---

## 🚀 Quick Start

### Installation

```r
# Install from GitHub
install.packages("devtools")
devtools::install_github("annaballestrazzi/uncoverappLib")
```

### One-Time Setup

```r
library(uncoverappLib)

# Download annotation files (~1GB, requires internet)
setup_uncoverapp()

# Verify setup completed
check_annotations()
```

### Launch App

```r
# Start interactive interface
uncoverappLib.run()
```

The app opens in your default browser. No server required—everything runs locally.

---

## 💡 Common Workflows

### Single Sample Analysis
**You have:** 1 coverage file + gene list

```r
buildInput(
  geneList = "genes.txt",
  genome = "hg38",
  chromosome_notation = "chr",
  sampleList = "coverage.list",
  outDir = "./results",
  type_input = "genes",
  type_coverage = "bed",
  input_coord_system = "0-based"
)
```

### Cohort Processing
**You have:** Multiple BAM files + gene panel

```r
buildInput(
  geneList = "panel_genes.txt",
  genome = "hg38",
  chromosome_notation = "chr",
  sampleList = "cohort_bams.list",
  outDir = "./cohort_results",
  type_input = "genes",
  type_coverage = "bam",
  MAPQ.min = 20,
  base.quality = 20
)
```

### Batch Annotation
**You have:** Coverage data, need variant annotations

```r
buildAnnotation(
  sample_data = "results/coverage.bed",
  target_sample = "patient_001",
  coverage_threshold = 20,
  genome = "hg38",
  output_formatted = "patient_001_annotated.xlsx"
)
```

**→** See [QUICK_START.md](QUICK_START.md) for all scenarios

---

## 📋 What You Need

### Required
- **R** ≥ 4.0.0
- **Java Runtime** (for Excel export)
- **2 GB disk space** (annotation cache)
- **4 GB RAM** minimum

### Input Files
- **Coverage data:** BAM files OR BED coverage files
- **Gene targets:** Gene list (HGNC symbols) OR custom BED regions
- **Reference:** hg19 or hg38

### Supported Formats
```
Gene list:        genes.txt (one symbol per line)
Target regions:   targets.bed (chr, start, end, name)
BAM files:        *.bam (with .bai index)
BED coverage:     *.bed (chr, start, end, depth)
Sample list:      *.list (absolute paths, one per line)
```

---

## 🎯 Key Features

### Interactive Mode
Launch the Shiny app for real-time exploration:
- Filter by gene, chromosome, region, or genome-wide
- Generate publication-quality coverage plots
- Calculate binomial probabilities for low-frequency variants
- Export formatted Excel files with color-coded annotations

### Batch Mode
Process large cohorts from the command line:
- Parallel processing for multiple samples
- Reproducible workflows with R scripts
- Integration with existing pipelines
- Automated quality reports

### Clinical Integration
- **OMIM gene highlighting** in annotation tables
- **ClinVar** pathogenicity classifications
- **ACMG criteria support** via prediction scores
- **Offline operation** for air-gapped systems

**Excel Color Coding:**
| Color | Meaning | Example |
|-------|---------|---------|
| 🔴 Red | Pathogenic predictions, ClinVar entries, rare variants | AF < 0.01, CADD > 20 |
| 🟢 Green | Benign predictions, common variants | AF > 0.01, SIFT Tolerated |
| 🟡 Yellow | High-impact variants (rare + pathogenic + high impact) | H/M + ClinVar + AF < 0.01 |
| 🔵 Blue | OMIM disease-associated genes | Gene in OMIM database |

---

## 📖 Documentation

**Quick References:**
- [QUICK_START.md](QUICK_START.md) - Choose your scenario, get the command
- [USER_GUIDE.md](USER_GUIDE.md) - Complete documentation with examples

**In-App Help:**
```r
# View package vignettes
browseVignettes("uncoverappLib")

# Function documentation
?buildInput
?buildAnnotation
```

**Online Resources:**
- **Paper:** https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1
- **Issues:** https://github.com/annaballestrazzi/uncoverappLib/issues
- **Zenodo:** https://zenodo.org/records/17524340

---

## 🔒 Offline Operation

**Network Requirements:**
- ✅ **Setup phase:** Internet required (one-time, ~10 min)
- ✅ **Daily use:** Fully offline after setup
- ✅ **Annotations:** Cached locally (~1GB)
- ✅ **Updates:** Manual (via GitHub)

**Perfect for:**
- Clinical diagnostic labs
- Secure research environments
- Teaching/training settings
- Field work on laptops

---

## 🐛 Troubleshooting

**"No annotation files found"**
```r
setup_uncoverapp()  # Re-download
check_annotations() # Verify
```

**"Gene not found"**
- Use official HGNC symbols (e.g., `TP53` not `p53`)
- Check genome version matches data (hg19 vs hg38)

**"No coverage data"**
- Verify chromosome naming ("chr1" vs "1")
- Check coordinate system (0-based vs 1-based for BED)
- Ensure BAM files are indexed (.bai present)

**"App is slow"**
- Use BED coverage files instead of BAM (10-50x faster)
- Process large gene lists in batch mode
- Increase quality filters to reduce data volume

**More help:** See [USER_GUIDE.md](USER_GUIDE.md) troubleshooting section

---

## 📊 Performance

| Task | BAM Input | BED Input |
|------|-----------|-----------|
| Single sample, 10 genes | 2-5 min | 10-30 sec |
| Single sample, 100 genes | 10-20 min | 1-2 min |
| 10 samples, 50 genes | 1-2 hours | 5-10 min |

**Recommendation:** Pre-compute coverage with `bedtools` or `mosdepth` for faster processing.

---

## 📜 Citation

If you use uncoverappLib in your research, please cite:

```
Ballestraz et al. (2020). unCOVERApp: a web application for 
clinical assessment of sequence coverage at the base-pair level.
bioRxiv. doi: 10.1101/2020.02.10.939769
```

---

## 🏗️ Version History

**v1.10.1** (February 2025)
- ✅ Offline-first design (system fonts, no external dependencies)
- ✅ Fixed gene coordinate retrieval for alternative contigs
- ✅ Enhanced UI responsiveness and layout
- ✅ Improved error messages and user feedback
- ✅ Better handling of edge cases in coverage calculation

**v1.10.0** (November 2024)
- ✅ BED coverage file support
- ✅ Coordinate system selection (0-based/1-based)
- ✅ Enhanced maxAF calculator
- ✅ Updated parameter names (`sampleList`, `chromosome_notation`)

---

## 📄 License

GPL-3 | See LICENSE file

---

## 👥 Contributors

**Original concept and first version:** Emanuela Iovino  
**Current development:** Anna Ballestrazzi

**Affiliation:** University of Bologna

**Acknowledgments:**
- dbNSFP database (variant annotations)
- UCSC Genome Browser (gene coordinates)
- gnomAD project (population frequencies)
- ClinVar (clinical variant database)