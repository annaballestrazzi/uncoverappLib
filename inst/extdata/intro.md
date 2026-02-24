# uncoverappLib

**Interactive coverage analysis and variant annotation for targeted clinical sequencing**

---

Reliable variant calling depends on reliable coverage. A variant that falls in a poorly sequenced region is a variant that may never be reported — not because it is absent, but because it was never observed with sufficient confidence to be called. uncoverappLib makes those gaps visible, annotated, and ready for clinical review.

---

## What this application does

Starting from a coverage file — generated here or imported from your pipeline — uncoverappLib identifies every genomic position below a sequencing depth threshold you define. For each of those positions it queries a local annotation database and returns all known variants: their population frequency in gnomAD, their pathogenicity in ClinVar, and functional impact predictions from CADD, MutationAssessor, SIFT, PolyPhen-2, and M-CAP.

The result is a prioritised, colour-coded table of positions where a variant could have been missed — annotated with the evidence needed to decide whether Sanger validation or re-sequencing is warranted.

---

## Three analysis tools

**Coverage Analysis** — the core of the application. Load a coverage file, filter by gene, chromosome, or custom region, set your depth threshold, and identify gaps. Generate a Gviz plot of the full gene with exons and coverage tracks. Annotate all low-coverage positions and export to Excel with conditional formatting.

**Maximum Credible Allele Frequency** — instead of a generic allele frequency cut-off, compute a gene-specific threshold based on the inheritance pattern, disease prevalence, allelic heterogeneity, and penetrance of the condition you are investigating. Variants below this threshold are flagged as potentially relevant.

**Binomial Probability** — given the observed coverage at a position and the expected allele fraction of a variant of interest, compute the probability of detecting at least N supporting reads. A direct answer to the question: *is this depth enough to reliably call this variant?*

---

## Before you begin

You need a **coverage file**: a tab-separated BED file with columns for chromosome, start, end, and per-sample coverage depth.

You can generate it directly here: go to the **Processing and Statistical Summary** tab, upload a gene list (or a target BED specifying custom chromosomal regions) and a list of BAM or pre-computed coverage files, and click *Process Coverage Files*. After processing, an optional filter panel lets you subset by sample and threshold and download the result as XLSX.

You can also generate the coverage file from the R console using `buildInput()` and load it here for interactive exploration — this is the recommended approach for large gene panels or many samples.

If this is your first time, the **Quick Start** and **User Guide** tabs in the menu walk you through every use case step by step.

---

*Everything runs locally on your computer. No sequencing data, no coverage files, and no patient information are transmitted to any external server.*

*Iovino E, Pippucci T, Ballestrazzi A (2020) · [doi:10.1101/2020.02.10.939769](https://www.biorxiv.org/content/10.1101/2020.02.10.939769v1)*
