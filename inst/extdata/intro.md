# UNCOVERAPP

## Interactive web application for clinical assessment of sequence coverage at base-pair level

This is a web application for clinical assessment of sequence coverage in diagnostic genomics.

---

## What You Can Do

### 📊 Coverage Analysis
Display interactive plots showing gene coverage down to base-pair resolution with functional and clinical annotations of positions within coverage gaps.

**Access:** "Coverage Analysis" page

**Features:**
- Visual inspection of coverage gaps
- Gene-level plots with exon structure
- Filter by gene, chromosome, or custom regions
- Identify positions requiring validation

---

### 🧬 Maximum Credible Allele Frequency (maxAF)
Calculate the maximum credible population allele frequency tailored to your disease model instead of using a general AF cut-off (e.g., 1% or 0.1%).

**Access:** "Calculate AF by allele frequency app" page

**Use for:**
- Rare disease genetics
- Gene-specific filtering
- Personalized AF thresholds based on inheritance pattern and penetrance

**Reference:** [maxAF Calculator](http://cardiodb.org/allelefrequencyapp/)

---

### 📈 Binomial Probability Calculator
Calculate the 95% probability of observing at least N variant-supporting reads based on expected allele fraction and sequencing depth.

**Access:** "Binomial distribution" page

**Especially useful for:**
- Somatic variants with low allele fraction
- Mosaic variants
- Determining detection probability at given coverage depth
- Planning validation experiments

---

## Input Requirements

All uncoverappLib functionalities require a BED file containing:

**Required columns:**
- Chromosome
- Start position (1-based)
- End position (1-based)
- Coverage value per sample

**Optional columns:**
- Reference and alternate allele counts (for BAM pileup data)

**How to generate:**
1. **Use the "Processing and Statistical Summary" page:**
   - Upload gene list + BAM/BED files
   - Click "Process Coverage Files"
   - Download generated BED file

2. **Or generate externally:**
   - Use `buildInput()` in R console
   - Use bedtools, mosdepth, or samtools
   - Load pre-computed coverage file directly

---

## Getting Started

### 1. Prepare Your Data

**Option A: Process from BAM/BED files**
- Go to "Processing and Statistical Summary" page
- Upload your gene list and sample files
- Click "Process Coverage Files"

**Option B: Load existing coverage file**
- Go to "Coverage Analysis" page
- Click "Select input file"
- Choose your coverage BED file

### 2. Analyze Coverage

**Step 1:** Select sample and set coverage threshold (e.g., 20x)

**Step 2:** Choose filter mode:
- **Gene name:** Analyze specific gene (most common)
- **Chromosome:** Analyze entire chromosome
- **Region:** Analyze specific coordinates
- **All chromosomes:** Genome-wide analysis

**Step 3:** Click "Calculate Low Coverage Regions"

**Step 4:** Click "Calculate Annotations" (optional)

### 3. Export Results

- Download annotated Excel files with color-coded pathogenicity predictions
- Export gene coverage plots for publications
- Save statistical summaries for quality control

---

## Key Features

✅ **Offline operation** - Works without internet after initial setup  
✅ **Multiple filtering modes** - Gene, chromosome, region, or genome-wide  
✅ **Clinical annotations** - ClinVar, CADD, gnomAD, OMIM integration  
✅ **Publication-ready plots** - Export high-quality gene coverage visualizations  
✅ **Batch processing** - Process multiple samples via R console  
✅ **Excel export** - Color-coded annotations for clinical review  

---

## Color Coding in Excel Files

| Color | Meaning |
|-------|---------|
| 🔴 Red | Pathogenic predictions, ClinVar pathogenic, rare variants (AF < 0.01) |
| 🟢 Green | Benign predictions, common variants (AF > 0.01) |
| 🟡 Yellow | High-impact variants (pathogenic + rare + high functional impact) |
| 🔵 Blue | OMIM disease-associated genes |

---

## Need Help?

**Documentation:**
- See "QUICK_START.md" for scenario-based commands
- See "USER_GUIDE.md" for complete documentation

**Function help:**
```r
?buildInput
?buildAnnotation
```

**Report issues:**
[GitHub Issues](https://github.com/annaballestrazzi/uncoverappLib/issues)

---

## Citation

If you use uncoverappLib in your research, please cite:

Ballestraz et al. (2020). unCOVERApp: a web application for clinical assessment of sequence coverage at the base-pair level. *bioRxiv*. doi: 10.1101/2020.02.10.939769

---

**Version:** 1.10.1  
**Last Updated:** February 2025

**Made with ❤️ for the clinical genomics community**