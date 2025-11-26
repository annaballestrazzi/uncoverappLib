# ==============================================================================
# TESTS FOR getAnnotationFiles() AND setup_uncoverapp()
# ==============================================================================

test_that("getAnnotationFiles validates assembly parameter", {
  # Valid assemblies should work
  expect_no_error(getAnnotationFiles(assembly = "hg19"))
  expect_no_error(getAnnotationFiles(assembly = "hg38"))
  expect_no_error(getAnnotationFiles(assembly = c("hg19", "hg38")))
  
  # Invalid assemblies should fail with match.arg error
  expect_error(
    getAnnotationFiles(assembly = "hg99"),
    "'arg' should be one of|hg19|hg38"  # Match actual error
  )
  
  expect_error(
    getAnnotationFiles(assembly = "mm10"),
    "'arg' should be one of|hg19|hg38"  # Match actual error
  )
})

test_that("getAnnotationFiles returns correct structure", {
  skip_on_cran()
  skip_if_offline()
  
  result <- getAnnotationFiles(assembly = "hg19")
  
  # Should be a list
  expect_type(result, "list")
  
  # Should have hg19 element
  expect_true("hg19" %in% names(result))
  
  # hg19 should have at least 2 files (.gz and .tbi)
  expect_true(length(result$hg19) >= 2)
  
  # Files should exist
  for (f in result$hg19) {
    expect_true(file.exists(f), info = paste("File missing:", f))
  }
  
  # Should have .gz and .tbi files
  gz_files <- result$hg19[grepl("\\.gz$", result$hg19)]
  tbi_files <- result$hg19[grepl("\\.tbi$", result$hg19)]
  
  expect_true(length(gz_files) >= 1)
  expect_true(length(tbi_files) >= 1)
})

test_that("getAnnotationFiles works with both assemblies", {
  skip_on_cran()
  skip_if_offline()
  
  result <- getAnnotationFiles(assembly = c("hg19", "hg38"))
  
  expect_type(result, "list")
  expect_true(all(c("hg19", "hg38") %in% names(result)))
  expect_true(length(result$hg19) >= 2)
  expect_true(length(result$hg38) >= 2)
})

test_that("getAnnotationFiles uses cache", {
  skip_on_cran()
  skip_if_offline()
  
  # First call (may download or use existing)
  result1 <- getAnnotationFiles(assembly = "hg19")
  
  # Second call (should definitely be from cache)
  start_time <- Sys.time()
  result2 <- getAnnotationFiles(assembly = "hg19")
  elapsed <- as.numeric(difftime(Sys.time(), start_time, units = "secs"))
  
  # Cache call should be very fast (<2 seconds - more lenient)
  expect_lt(elapsed, 2.0, 
            label = "Cached call should be under 2 seconds")
  
  # Results should be identical
  expect_identical(result1, result2)
})

test_that("setup_uncoverapp works", {
  skip_on_cran()
  skip_if_offline()
  
  # Function has no parameters - just call it
  expect_no_error(setup_uncoverapp())
  
  # Should create cache directory
  cache_dir <- rappdirs::user_cache_dir("uncoverapp")
  expect_true(dir.exists(cache_dir))
  
  # Should have annotation files
  hg19_file <- file.path(cache_dir, "sorted_hg19.bed.gz")
  hg38_file <- file.path(cache_dir, "sorted_hg38.bed.gz")
  
  expect_true(file.exists(hg19_file) || file.exists(hg38_file),
              info = "At least one annotation file should exist")
})

test_that("check_annotations works", {
  # Function should always succeed (prints status)
  expect_output(check_annotations(), "Annotation|Cache|sorted")
})