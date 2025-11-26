## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----coverage-----------------------------------------------------------------
# library(uncoverappLib)
# 
# buildInput(
#   geneList = "genes.txt",
#   bamList = "samples.list",
#   genome = "hg19",
#   type_bam = "chr",
#   type_input = "genes",
#   type_coverage = "bam",
#   outDir = "./output"
# )
# 
# # Output: ./output/output/DATE.bed

## ----annotate-----------------------------------------------------------------
# annotate_all_lowcov(
#   sample_data = "./output/output/Mon_Nov_11_2024.bed",
#   target_sample = "sample1",
#   coverage_threshold = 20,
#   genome = "hg19",
#   output_formatted = "sample1_annotated.xlsx"
# )

## ----example1-----------------------------------------------------------------
# # Input files
# cat("BRCA1\nBRCA2\nTP53\nATM\n", file = "panel.txt")
# cat("/data/patient001.bam\n", file = "patient001.list")
# 
# # Step 1: Coverage
# buildInput(
#   geneList = "panel.txt",
#   bamList = "patient001.list",
#   genome = "hg19",
#   type_bam = "chr",
#   type_input = "genes",
#   type_coverage = "bam",
#   outDir = "./patient001",
#   MAPQ.min = 20,
#   base.quality = 20
# )
# 
# # Step 2: Annotate
# coverage_file <- list.files(
#   "./patient001/output",
#   pattern = "\\.bed$",
#   full.names = TRUE
# )[1]
# 
# annotate_all_lowcov(
#   sample_data = coverage_file,
#   target_sample = "patient001",
#   coverage_threshold = 20,
#   genome = "hg19",
#   output_formatted = "patient001_lowcov.xlsx"
# )

## ----example2-----------------------------------------------------------------
# # Input
# samples <- c("patient001", "patient002", "patient003")
# 
# # Create BAM list
# bam_files <- paste0("/data/", samples, ".bam")
# writeLines(bam_files, "cohort.list")
# 
# # Step 1: Coverage (all samples together)
# buildInput(
#   geneList = "panel.txt",
#   bamList = "cohort.list",
#   genome = "hg19",
#   type_bam = "chr",
#   type_input = "genes",
#   type_coverage = "bam",
#   outDir = "./cohort"
# )
# 
# # Step 2: Annotate each sample
# coverage_file <- list.files("./cohort/output", pattern = "\\.bed$", full.names = TRUE)[1]
# 
# for (sample in samples) {
#   cat("Processing", sample, "...\n")
# 
#   annotate_all_lowcov(
#     sample_data = coverage_file,
#     target_sample = sample,
#     coverage_threshold = 20,
#     genome = "hg19",
#     output_formatted = paste0(sample, "_annotated.xlsx")
#   )
# }

## ----example3-----------------------------------------------------------------
# # If you have pre-computed coverage files
# # (e.g., from `samtools depth`)
# 
# # bedcov.list:
# # /data/sample1_coverage.bed
# # /data/sample2_coverage.bed
# 
# buildInput(
#   geneList = "panel.txt",
#   bamList = "bedcov.list",
#   genome = "hg19",
#   type_bam = "chr",
#   type_input = "genes",
#   type_coverage = "bed",          # ← BED instead of BAM
#   input_coord_system = "0-based", # ← BED uses 0-based
#   outDir = "./output"
# )

## ----parallel-----------------------------------------------------------------
# library(parallel)
# 
# # Annotate samples in parallel
# mclapply(samples, function(s) {
#   annotate_all_lowcov(
#     sample_data = coverage_file,
#     target_sample = s,
#     coverage_threshold = 20
#   )
# }, mc.cores = 4)

## ----custom_ann---------------------------------------------------------------
# # Use local annotation database
# annotate_all_lowcov(
#   sample_data = "coverage.bed",
#   target_sample = "sample1",
#   annotation_file = "/custom/sorted_hg19.bed.gz"
# )

## ----omim---------------------------------------------------------------------
# # Add OMIM data during buildInput
# buildInput(
#   geneList = "genes.txt",
#   bamList = "samples.list",
#   genome = "hg19",
#   type_bam = "chr",
#   type_input = "genes",
#   type_coverage = "bam",
#   outDir = "./output",
#   annotation_file = "omim_phenotypes.txt"  # ← Tab-separated with SYMBOL column
# )

## ----download_ann-------------------------------------------------------------
# # Download annotation databases
# getAnnotationFiles()

## ----sessioninfo, eval=TRUE---------------------------------------------------
sessionInfo()

