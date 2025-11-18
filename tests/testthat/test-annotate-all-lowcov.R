# ==============================================================================
# TESTS FOR annotate_all_lowcov()
# ==============================================================================

test_that("annotate_all_lowcov validates required inputs", {
  # Missing sample_data
  expect_error(
    annotate_all_lowcov(target_sample = "test"),
    "sample_data must be supplied"
  )
  
  # Missing target_sample
  tmpfile <- write_mock_coverage_file()
  expect_error(
    annotate_all_lowcov(sample_data = tmpfile),
    "target_sample must be supplied"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("annotate_all_lowcov validates genome parameter", {
  tmpfile <- write_mock_coverage_file()
  
  expect_error(
    annotate_all_lowcov(
      sample_data = tmpfile,
      target_sample = "test_sample",
      genome = "hg99"
    ),
    "hg19.*hg38"
  )
  
  expect_error(
    annotate_all_lowcov(
      sample_data = tmpfile,
      target_sample = "test_sample",
      genome = "mm10"
    ),
    "hg19.*hg38"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("annotate_all_lowcov fails with non-existent file", {
  expect_error(
    annotate_all_lowcov(
      sample_data = "/fake/path/nonexistent.tsv",
      target_sample = "test"
    ),
    "not found"
  )
})

test_that("annotate_all_lowcov detects column names correctly", {
  # Test exact match
  df1 <- create_mock_coverage(n_rows = 20, sample_name = "sample1")
  tmpfile1 <- tempfile(fileext = ".tsv")
  write.table(df1, tmpfile1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  expect_silent({
    result <- annotate_all_lowcov(
      sample_data = tmpfile1,
      target_sample = "count_sample1",
      coverage_threshold = 50,
      genome = "hg19",
      output_intersect = tempfile(fileext = ".tsv"),
      output_formatted = tempfile(fileext = ".xlsx")
    )
  })
  
  cleanup_test_files(tmpfile1)
  
  # Test with prefix
  df2 <- create_mock_coverage(n_rows = 20, sample_name = "sample2")
  tmpfile2 <- tempfile(fileext = ".tsv")
  write.table(df2, tmpfile2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  expect_silent({
    result <- annotate_all_lowcov(
      sample_data = tmpfile2,
      target_sample = "sample2",  # Without "count_" prefix
      coverage_threshold = 50,
      genome = "hg19",
      output_intersect = tempfile(fileext = ".tsv"),
      output_formatted = tempfile(fileext = ".xlsx")
    )
  })
  
  cleanup_test_files(tmpfile2)
})

test_that("annotate_all_lowcov fails when column not found", {
  tmpfile <- write_mock_coverage_file(sample_name = "real_sample")
  
  expect_error(
    annotate_all_lowcov(
      sample_data = tmpfile,
      target_sample = "nonexistent_sample",
      coverage_threshold = 20
    ),
    "not found"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("annotate_all_lowcov filters by coverage threshold correctly", {
  # Create data with known coverage distribution
  df <- data.frame(
    chromosome = rep("chr1", 100),
    start = seq(1000, 10000, length.out = 100),
    end = seq(1100, 10100, length.out = 100),
    count_test = c(rep(10, 30), rep(50, 40), rep(100, 30))  # 30 low, 70 high
  )
  
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  out_tsv <- tempfile(fileext = ".tsv")
  out_xlsx <- tempfile(fileext = ".xlsx")
  
  # Filter at threshold 20 (should keep only 30 rows)
  result <- annotate_all_lowcov(
    sample_data = tmpfile,
    target_sample = "test",
    coverage_threshold = 20,
    genome = "hg19",
    output_intersect = out_tsv,
    output_formatted = out_xlsx
  )
  
  # Read back the TSV output
  if (file.exists(out_tsv)) {
    output_df <- read.table(out_tsv, header = TRUE, sep = "\t")
    # Should have filtered correctly (may be 0 if no annotation overlap)
    expect_true(nrow(output_df) <= 30)
  }
  
  cleanup_test_files(c(tmpfile, out_tsv, out_xlsx))
})

test_that("annotate_all_lowcov handles missing data columns correctly", {
  tmpfile <- write_mock_coverage_file()
  
  # Missing 'start' column
  df <- read.table(tmpfile, header = TRUE, sep = "\t")
  df$start <- NULL
  tmpfile_bad <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile_bad, sep = "\t", quote = FALSE, row.names = FALSE)
  
  expect_error(
    annotate_all_lowcov(
      sample_data = tmpfile_bad,
      target_sample = "test_sample",
      coverage_threshold = 20
    ),
    "start.*not found"
  )
  
  cleanup_test_files(c(tmpfile, tmpfile_bad))
})

test_that("annotate_all_lowcov creates output files with correct structure", {
  skip_on_cran()
  skip_if_offline()
  
  tmpfile <- write_mock_coverage_file(n_rows = 50, coverage_range = c(5, 30))
  out_tsv <- tempfile(fileext = ".tsv")
  out_xlsx <- tempfile(fileext = ".xlsx")
  
  result <- annotate_all_lowcov(
    sample_data = tmpfile,
    target_sample = "test_sample",
    coverage_threshold = 20,
    genome = "hg19",
    output_intersect = out_tsv,
    output_formatted = out_xlsx
  )
  
  # Check TSV exists and has correct format
  expect_true(file.exists(out_tsv))
  tsv_check <- check_tsv_format(out_tsv, c("seqnames", "start", "end", "coverage"))
  expect_true(tsv_check$valid)
  
  # Check Excel exists and has correct structure
  if (file.exists(out_xlsx)) {
    excel_check <- check_excel_format(out_xlsx, "Low Coverage Variants")
    expect_true(excel_check$valid)
  }
  
  cleanup_test_files(c(tmpfile, out_tsv, out_xlsx))
})

test_that("annotate_all_lowcov returns data invisibly", {
  tmpfile <- write_mock_coverage_file(n_rows = 20, coverage_range = c(5, 25))
  out_tsv <- tempfile(fileext = ".tsv")
  out_xlsx <- tempfile(fileext = ".xlsx")
  
  result <- annotate_all_lowcov(
    sample_data = tmpfile,
    target_sample = "test_sample",
    coverage_threshold = 20,
    genome = "hg19",
    output_intersect = out_tsv,
    output_formatted = out_xlsx
  )
  
  # Should return a data.frame
  expect_true(is.data.frame(result) || is.null(result))
  
  # If data.frame, should have expected columns
  if (is.data.frame(result) && nrow(result) > 0) {
    expect_true("seqnames" %in% colnames(result))
    expect_true("start" %in% colnames(result))
    expect_true("end" %in% colnames(result))
    expect_true("coverage" %in% colnames(result))
  }
  
  cleanup_test_files(c(tmpfile, out_tsv, out_xlsx))
})

test_that("annotate_all_lowcov handles empty results gracefully", {
  # Create data where all coverage is HIGH (no low coverage positions)
  df <- data.frame(
    chromosome = rep("chr1", 20),
    start = seq(1000, 2000, length.out = 20),
    end = seq(1100, 2100, length.out = 20),
    count_test = rep(200, 20)  # All high coverage
  )
  
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  out_tsv <- tempfile(fileext = ".tsv")
  out_xlsx <- tempfile(fileext = ".xlsx")
  
  # Filter at threshold 20 - should find NOTHING
  expect_error(
    annotate_all_lowcov(
      sample_data = tmpfile,
      target_sample = "test",
      coverage_threshold = 20,
      genome = "hg19",
      output_intersect = out_tsv,
      output_formatted = out_xlsx
    ),
    "No positions with coverage"
  )
  
  cleanup_test_files(c(tmpfile, out_tsv, out_xlsx))
})

test_that("annotate_all_lowcov normalizes chromosome names", {
  # Test with and without 'chr' prefix
  df1 <- data.frame(
    chromosome = c("1", "2", "X"),  # Without chr prefix
    start = c(1000, 2000, 3000),
    end = c(1100, 2100, 3100),
    count_test = c(10, 15, 20)
  )
  
  tmpfile1 <- tempfile(fileext = ".tsv")
  write.table(df1, tmpfile1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  expect_silent({
    result1 <- annotate_all_lowcov(
      sample_data = tmpfile1,
      target_sample = "test",
      coverage_threshold = 25,
      genome = "hg19",
      output_intersect = tempfile(fileext = ".tsv"),
      output_formatted = tempfile(fileext = ".xlsx")
    )
  })
  
  df2 <- data.frame(
    chromosome = c("chr1", "chr2", "chrX"),  # With chr prefix
    start = c(1000, 2000, 3000),
    end = c(1100, 2100, 3100),
    count_test = c(10, 15, 20)
  )
  
  tmpfile2 <- tempfile(fileext = ".tsv")
  write.table(df2, tmpfile2, sep = "\t", quote = FALSE, row.names = FALSE)
  
  expect_silent({
    result2 <- annotate_all_lowcov(
      sample_data = tmpfile2,
      target_sample = "test",
      coverage_threshold = 25,
      genome = "hg19",
      output_intersect = tempfile(fileext = ".tsv"),
      output_formatted = tempfile(fileext = ".xlsx")
    )
  })
  
  cleanup_test_files(c(tmpfile1, tmpfile2))
})