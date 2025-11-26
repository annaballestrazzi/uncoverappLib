# ==============================================================================
# HELPER FUNCTIONS FOR TESTING
# ==============================================================================

# Skip conditions
skip_if_no_annotation <- function() {
  cache_dir <- rappdirs::user_cache_dir("uncoverapp")
  hg19_file <- file.path(cache_dir, "sorted_hg19.bed.gz")
  
  if (!file.exists(hg19_file)) {
    skip("Annotation files not available. Run setup_uncoverapp() first.")
  }
}

skip_if_offline <- function() {
  if (!curl::has_internet()) {
    skip("No internet connection")
  }
}

# ==============================================================================
# MOCK DATA GENERATORS
# ==============================================================================

#' Create mock coverage data
create_mock_coverage <- function(n_rows = 100, 
                                  coverage_range = c(5, 200),
                                  sample_name = "test_sample") {
  set.seed(123)  # Reproducible
  
  data.frame(
    chromosome = paste0("chr", sample(1:22, n_rows, replace = TRUE)),
    start = sort(sample(10000000:20000000, n_rows)),
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

#' Write mock coverage to file
write_mock_coverage_file <- function(n_rows = 100, 
                                      coverage_range = c(5, 200),
                                      sample_name = "test_sample") {
  df <- create_mock_coverage(n_rows, coverage_range, sample_name)
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = TRUE)
  tmpfile
}

#' Create mock gene list
create_mock_gene_list <- function(genes = c("BRCA1", "TP53", "EGFR")) {
  tmpfile <- tempfile(fileext = ".txt")
  writeLines(genes, tmpfile)
  tmpfile
}

#' Create mock sample list
create_mock_sample_list <- function(n_samples = 2, type = "bam") {
  tmpfile <- tempfile(fileext = ".list")
  ext <- if (type == "bam") ".bam" else ".bed"
  sample_paths <- paste0(tempdir(), "/sample", seq_len(n_samples), ext)
  writeLines(sample_paths, tmpfile)
  tmpfile
}

#' Create mock target BED
create_mock_target_bed <- function() {
  tmpfile <- tempfile(fileext = ".bed")
  df <- data.frame(
    chr = c("chr17", "chr17", "chr13"),
    start = c(41196312, 7571720, 32889617),
    end = c(41277500, 7590868, 32973809),
    gene = c("BRCA1", "TP53", "BRCA2")
  )
  write.table(df, tmpfile, sep = "\t", quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  tmpfile
}

# ==============================================================================
# CLEANUP
# ==============================================================================

cleanup_test_files <- function(file_list) {
  for (f in file_list) {
    if (file.exists(f)) {
      unlink(f, force = TRUE, recursive = TRUE)
    }
  }
  invisible(NULL)
}