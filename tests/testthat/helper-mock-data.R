# ==============================================================================
# HELPER FUNCTIONS AND MOCK DATA FOR TESTING
# ==============================================================================

# ==============================================================================
# MOCK DATA GENERATORS
# ==============================================================================

#' Create mock coverage file
#' @param n_rows Number of rows
#' @param coverage_range Range of coverage values
#' @param sample_name Name of sample
create_mock_coverage <- function(n_rows = 100, 
                                  coverage_range = c(5, 200),
                                  sample_name = "test_sample") {
  data.frame(
    chromosome = paste0("chr", sample(1:22, n_rows, replace = TRUE)),
    start = sample(1000:100000, n_rows),
    end = sample(1000:100000, n_rows),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(
      end = start + sample(10:100, n_rows, replace = TRUE),
      !!paste0("count_", sample_name) := sample(
        coverage_range[1]:coverage_range[2], 
        n_rows, 
        replace = TRUE
      )
    ) %>%
    dplyr::arrange(chromosome, start)
}

#' Create mock gene list file
create_mock_gene_list <- function(genes = c("BRCA1", "TP53", "EGFR")) {
  tmpfile <- tempfile(fileext = ".txt")
  writeLines(genes, tmpfile)
  tmpfile
}

#' Create mock BAM list file
create_mock_bam_list <- function(n_samples = 2) {
  tmpfile <- tempfile(fileext = ".list")
  # Create fake BAM paths
  bam_paths <- paste0("/fake/path/sample", seq_len(n_samples), ".bam")
  writeLines(bam_paths, tmpfile)
  tmpfile
}

#' Create mock annotation file (BED format)
#' @param n_variants Number of variants
create_mock_annotation <- function(n_variants = 50) {
  data.frame(
    chromosome = paste0("chr", sample(1:22, n_variants, replace = TRUE)),
    start = sample(1000:100000, n_variants),
    end = sample(1000:100000, n_variants),
    REF = sample(c("A", "T", "C", "G"), n_variants, replace = TRUE),
    ALT = sample(c("A", "T", "C", "G"), n_variants, replace = TRUE),
    dbsnp = paste0("rs", sample(1000:9999, n_variants)),
    GENENAME = sample(c("BRCA1", "TP53", "EGFR", "KRAS"), n_variants, replace = TRUE),
    PROTEIN_ensembl = paste0("ENSP", sample(10000:99999, n_variants)),
    field9 = ".",
    MutationAssessor = sample(c("H", "M", "L", "N"), n_variants, replace = TRUE),
    SIFT = sample(c("D", "T"), n_variants, replace = TRUE),
    Polyphen2 = sample(c("D", "P", "B"), n_variants, replace = TRUE),
    M_CAP = sample(c("D", "T"), n_variants, replace = TRUE),
    CADD_PHED = runif(n_variants, 0, 30),
    AF_gnomAD = runif(n_variants, 0, 0.5),
    ClinVar = sample(c(".", "Pathogenic", "Likely_pathogenic"), n_variants, replace = TRUE),
    clinvar_MedGen_id = paste0("C", sample(100000:999999, n_variants)),
    HGVSc_VEP = paste0("c.", sample(1:1000, n_variants), "A>G"),
    HGVSp_VEP = paste0("p.Arg", sample(1:500, n_variants), "His"),
    stringsAsFactors = FALSE
  ) %>%
    dplyr::mutate(end = start + 1) %>%
    dplyr::arrange(chromosome, start)
}

# ==============================================================================
# MOCK FILE SYSTEM FUNCTIONS
# ==============================================================================

#' Write mock coverage file to temp location
#' @return Path to temp file
write_mock_coverage_file <- function(n_rows = 100, 
                                      coverage_range = c(5, 200),
                                      sample_name = "test_sample") {
  df <- create_mock_coverage(n_rows, coverage_range, sample_name)
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  tmpfile
}

#' Write mock annotation file (compressed + indexed)
#' @return List with paths to .bed.gz and .bed.gz.tbi
write_mock_annotation_file <- function(n_variants = 50) {
  df <- create_mock_annotation(n_variants)
  
  # Write to temp BED file
  tmpbed <- tempfile(fileext = ".bed")
  write.table(df, tmpbed, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  # Compress with bgzip (if available)
  tmpgz <- paste0(tmpbed, ".gz")
  
  # For testing, we'll just use gzip (bgzip may not be available in test env)
  # In real usage, annotation files are properly bgzipped and tabix-indexed
  R.utils::gzip(tmpbed, tmpgz, remove = TRUE, overwrite = TRUE)
  
  # Create fake index file
  tmpidx <- paste0(tmpgz, ".tbi")
  file.create(tmpidx)
  
  list(bed = tmpgz, tbi = tmpidx)
}

# ==============================================================================
# VALIDATION HELPERS
# ==============================================================================

#' Check if file is valid TSV with expected columns
check_tsv_format <- function(filepath, expected_cols) {
  if (!file.exists(filepath)) {
    return(list(valid = FALSE, reason = "File does not exist"))
  }
  
  df <- tryCatch({
    read.table(filepath, header = TRUE, sep = "\t", stringsAsFactors = FALSE, nrows = 10)
  }, error = function(e) {
    return(list(valid = FALSE, reason = paste("Read error:", e$message)))
  })
  
  if (!all(expected_cols %in% colnames(df))) {
    missing <- setdiff(expected_cols, colnames(df))
    return(list(valid = FALSE, reason = paste("Missing columns:", paste(missing, collapse = ", "))))
  }
  
  list(valid = TRUE)
}

#' Check if Excel file has expected structure
check_excel_format <- function(filepath, expected_sheet = "Low Coverage Variants") {
  if (!file.exists(filepath)) {
    return(list(valid = FALSE, reason = "File does not exist"))
  }
  
  sheets <- tryCatch({
    openxlsx::getSheetNames(filepath)
  }, error = function(e) {
    return(list(valid = FALSE, reason = paste("Cannot read Excel:", e$message)))
  })
  
  if (!expected_sheet %in% sheets) {
    return(list(valid = FALSE, reason = paste("Sheet", expected_sheet, "not found")))
  }
  
  list(valid = TRUE, sheets = sheets)
}

# ==============================================================================
# CLEANUP HELPERS
# ==============================================================================

#' Clean up temporary test files
cleanup_test_files <- function(file_list) {
  for (f in file_list) {
    if (file.exists(f)) {
      unlink(f, force = TRUE)
    }
  }
}