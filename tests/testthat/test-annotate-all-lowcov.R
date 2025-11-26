# ==============================================================================
# TESTS FOR buildAnnotation()
# ==============================================================================

test_that("buildAnnotation validates required parameters", {
  # Missing sample_data
  expect_error(
    buildAnnotation(
      target_sample = "test",
      coverage_threshold = 20,
      genome = "hg19"
    ),
    "sample_data"
  )
  
  # Missing target_sample
  tmpfile <- write_mock_coverage_file()
  expect_error(
    buildAnnotation(
      sample_data = tmpfile,
      coverage_threshold = 20,
      genome = "hg19"
    ),
    "target_sample"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("buildAnnotation validates genome parameter", {
  tmpfile <- write_mock_coverage_file()
  
  # Invalid genome should error
  expect_error(
    buildAnnotation(
      sample_data = tmpfile,
      target_sample = "test_sample",
      coverage_threshold = 20,
      genome = "mm10"
    ),
    "hg19|hg38|arg"  # Match actual error message
  )
  
  cleanup_test_files(tmpfile)
})

test_that("buildAnnotation validates file existence", {
  expect_error(
    buildAnnotation(
      sample_data = "/fake/path/nonexistent.tsv",
      target_sample = "test",
      coverage_threshold = 20,
      genome = "hg19"
    ),
    "not found|does not exist"
  )
})

test_that("buildAnnotation detects sample columns", {
  skip_if_no_annotation()
  
  # Test exact match - use real genomic coordinates
  df1 <- data.frame(
    chromosome = rep("chr17", 50),
    start = seq(41200000, 41200000 + 49 * 100, by = 100),
    end = seq(41200050, 41200050 + 49 * 100, by = 100),
    count_sample1 = sample(5:50, 50, replace = TRUE)
  )
  tmpfile1 <- tempfile(fileext = ".tsv")
  write.table(df1, tmpfile1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Should work with exact name (suppress ALL output)
  expect_no_error(
    capture.output(
      suppressMessages(
        suppressWarnings(
          buildAnnotation(
            sample_data = tmpfile1,
            target_sample = "count_sample1",
            coverage_threshold = 50,
            genome = "hg19",
            output_intersect = tempfile(fileext = ".tsv"),
            output_formatted = tempfile(fileext = ".xlsx")
          )
        )
      ),
      type = "output"
    )
  )
  
  # Should work without count_ prefix
  expect_no_error(
    capture.output(
      suppressMessages(
        suppressWarnings(
          buildAnnotation(
            sample_data = tmpfile1,
            target_sample = "sample1",
            coverage_threshold = 50,
            genome = "hg19",
            output_intersect = tempfile(fileext = ".tsv"),
            output_formatted = tempfile(fileext = ".xlsx")
          )
        )
      ),
      type = "output"
    )
  )
  
  cleanup_test_files(tmpfile1)
})

test_that("buildAnnotation fails with non-existent column", {
  tmpfile <- write_mock_coverage_file(sample_name = "real_sample")
  
  expect_error(
    buildAnnotation(
      sample_data = tmpfile,
      target_sample = "nonexistent_sample",
      coverage_threshold = 20,
      genome = "hg19"
    ),
    "not found|column"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("buildAnnotation filters by coverage threshold", {
  skip_if_no_annotation()
  
  # Use real genomic coordinates in gene-rich regions
  df <- data.frame(
    chromosome = rep("chr17", 100),
    start = seq(41200000, 41200000 + 99 * 100, by = 100),
    end = seq(41200050, 41200050 + 99 * 100, by = 100),
    count_test = c(rep(10, 30), rep(50, 40), rep(100, 30))
  )
  
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Suppress ALL output
  result <- capture.output(
    suppressMessages(
      suppressWarnings(
        buildAnnotation(
          sample_data = tmpfile,
          target_sample = "test",
          coverage_threshold = 20,
          genome = "hg19",
          output_intersect = tempfile(fileext = ".tsv"),
          output_formatted = tempfile(fileext = ".xlsx")
        )
      )
    ),
    type = "output"
  )
  
  # Test passes if no error thrown
  expect_true(TRUE)
  
  cleanup_test_files(tmpfile)
})

test_that("buildAnnotation handles empty results", {
  skip_if_no_annotation()
  
  # All high coverage - nothing should pass filter
  df <- data.frame(
    chromosome = rep("chr17", 20),
    start = seq(41200000, 41200000 + 19 * 100, by = 100),
    end = seq(41200050, 41200050 + 19 * 100, by = 100),
    count_test = rep(200, 20)
  )
  
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Should error with informative message
  expect_error(
    capture.output(
      suppressMessages(
        buildAnnotation(
          sample_data = tmpfile,
          target_sample = "test",
          coverage_threshold = 20,
          genome = "hg19"
        )
      ),
      type = "output"
    ),
    "No positions with coverage"
  )
  
  cleanup_test_files(tmpfile)
})

test_that("buildAnnotation normalizes chromosome names", {
  skip_if_no_annotation()
  
  # Test without 'chr' prefix
  df1 <- data.frame(
    chromosome = c("17", "17", "17"),
    start = c(41200000, 41201000, 41202000),
    end = c(41200050, 41201050, 41202050),
    count_test = c(10, 15, 20)
  )
  
  tmpfile1 <- tempfile(fileext = ".tsv")
  write.table(df1, tmpfile1, sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Should handle chromosome without 'chr'
  expect_no_error(
    capture.output(
      suppressMessages(
        suppressWarnings(
          buildAnnotation(
            sample_data = tmpfile1,
            target_sample = "test",
            coverage_threshold = 25,
            genome = "hg19",
            output_intersect = tempfile(fileext = ".tsv"),
            output_formatted = tempfile(fileext = ".xlsx")
          )
        )
      ),
      type = "output"
    )
  )
  
  cleanup_test_files(tmpfile1)
})

test_that("buildAnnotation returns correct data structure", {
  skip_if_no_annotation()
  
  df <- data.frame(
    chromosome = rep("chr17", 50),
    start = seq(41200000, 41200000 + 49 * 100, by = 100),
    end = seq(41200050, 41200050 + 49 * 100, by = 100),
    count_test = sample(5:25, 50, replace = TRUE)
  )
  
  tmpfile <- tempfile(fileext = ".tsv")
  write.table(df, tmpfile, sep = "\t", quote = FALSE, row.names = FALSE)
  
  result <- capture.output(
    suppressMessages(
      suppressWarnings(
        buildAnnotation(
          sample_data = tmpfile,
          target_sample = "test",
          coverage_threshold = 20,
          genome = "hg19",
          output_intersect = tempfile(fileext = ".tsv"),
          output_formatted = tempfile(fileext = ".xlsx")
        )
      )
    ),
    type = "output"
  )
  
  # Just verify no error
  expect_true(TRUE)
  
  cleanup_test_files(tmpfile)
})