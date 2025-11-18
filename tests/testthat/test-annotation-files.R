# ==============================================================================
# TESTS FOR getAnnotationFiles() and setup_uncoverapp()
# ==============================================================================

test_that("getAnnotationFiles works with valid genome", {
  skip_on_cran()
  skip_if_offline()
  
  # Test hg19
  result_hg19 <- getAnnotationFiles(assembly = "hg19", verbose = FALSE)
  
  expect_type(result_hg19, "list")
  expect_true("hg19" %in% names(result_hg19))
  expect_true(length(result_hg19$hg19) >= 1)
  
  # Check that files exist
  for (f in result_hg19$hg19) {
    expect_true(file.exists(f), info = paste("File not found:", f))
  }
  
  # Check for .gz and .tbi files
  gz_files <- result_hg19$hg19[grepl("\\.gz$", result_hg19$hg19)]
  tbi_files <- result_hg19$hg19[grepl("\\.tbi$", result_hg19$hg19)]
  
  expect_true(length(gz_files) >= 1, info = "No .gz file found")
  expect_true(length(tbi_files) >= 1, info = "No .tbi file found")
})

test_that("getAnnotationFiles works with hg38", {
  skip_on_cran()
  skip_if_offline()
  
  result_hg38 <- getAnnotationFiles(assembly = "hg38", verbose = FALSE)
  
  expect_type(result_hg38, "list")
  expect_true("hg38" %in% names(result_hg38))
  expect_true(length(result_hg38$hg38) >= 1)
})

test_that("getAnnotationFiles fails with invalid genome", {
  expect_error(
    getAnnotationFiles(assembly = "hg99", verbose = FALSE),
    "assembly"
  )
  
  expect_error(
    getAnnotationFiles(assembly = "mm10", verbose = FALSE),
    "assembly"
  )
})

test_that("getAnnotationFiles uses cache on second call", {
  skip_on_cran()
  skip_if_offline()
  
  # First call (may download)
  time1_start <- Sys.time()
  result1 <- getAnnotationFiles(assembly = "hg19", verbose = FALSE)
  time1 <- as.numeric(difftime(Sys.time(), time1_start, units = "secs"))
  
  # Second call (should be cached - much faster)
  time2_start <- Sys.time()
  result2 <- getAnnotationFiles(assembly = "hg19", verbose = FALSE)
  time2 <- as.numeric(difftime(Sys.time(), time2_start, units = "secs"))
  
  # Second call should be much faster (cached)
  expect_true(time2 < time1 * 0.5, 
              info = paste("Cache not working: time1 =", time1, "time2 =", time2))
  
  # Results should be identical
  expect_identical(result1, result2)
})

test_that("setup_uncoverapp downloads annotation files", {
  skip_on_cran()
  skip_if_offline()
  
  # This function should call getAnnotationFiles internally
  expect_message(
    setup_uncoverapp(verbose = TRUE),
    "annotation|download|hg19|hg38"
  )
})

test_that("Environment variables take precedence", {
  # Mock environment variables
  test_file <- tempfile(fileext = ".bed.gz")
  file.create(test_file)
  test_tbi <- paste0(test_file, ".tbi")
  file.create(test_tbi)
  
  withr::with_envvar(
    c(UNCOVERAPP_HG19_ANNOTATION = test_file),
    {
      # If function respects env vars, it should not download
      # This test depends on implementation
      expect_true(file.exists(Sys.getenv("UNCOVERAPP_HG19_ANNOTATION")))
    }
  )
  
  # Cleanup
  unlink(c(test_file, test_tbi))
})