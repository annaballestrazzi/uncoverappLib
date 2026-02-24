# Quick Start — uncoverappLib

---

## Install

### Step 1 — Install uncoverappLib from GitHub

`BiocManager` resolves both Bioconductor and CRAN dependencies in a single step:

```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("annaballestrazzi/uncoverappLib")
```

All dependencies — Bioconductor packages (`Gviz`, `Homo.sapiens`, `org.Hs.eg.db`, `TxDb.Hsapiens.UCSC.hg19.knownGene`, `TxDb.Hsapiens.UCSC.hg38.knownGene`, `EnsDb.Hsapiens.v75`, `EnsDb.Hsapiens.v86`, `OrganismDbi`, `Rsamtools`, `GenomicRanges`, `AnnotationDbi`, `GenomeInfoDb`, `IRanges`, `S4Vectors`) and CRAN packages (`shiny`, `shinyjs`, `shinyBS`, `shinyWidgets`, `shinycssloaders`, `DT`, `openxlsx`, `condformat`, `stringr`, `rappdirs`, `rlist`, `processx`, `dplyr`, `tidyselect`, `waiter`, `markdown`) are installed automatically.

### Step 2 — Download annotation files (one-time, ~1 GB)

```r
library(uncoverappLib)
setup_uncoverapp()    # downloads dbNSFP/VEP annotation files — 10–20 min
check_annotations()   # verifies everything is in place
```

After this step the app works fully offline. You never need to repeat this unless you want to update the annotation files.

---

## Launch the app

```r
library(uncoverappLib)
uncoverappLib.run()
```

The app opens in your browser. No data leaves your machine.

---

**For all usage scenarios — both interactive and standalone — see the [User Guide](USER_GUIDE.md).**