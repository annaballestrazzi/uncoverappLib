# tests/testthat/test-getAnnotationFiles.R
# Test essenziali per la funzione getAnnotationFiles()

test_that("getAnnotationFiles validates assembly parameter", {
  # Assembly invalido
  expect_error(
    getAnnotationFiles(assembly = "hg99"),
    "hg19|hg38"
  )
  
  # Assembly vuoto (matcha il vero messaggio di errore)
  expect_error(
    getAnnotationFiles(assembly = ""),
    "should be one of"
  )
})


test_that("getAnnotationFiles returns valid structure for hg19", {
  skip_on_cran()
  skip_if_offline()
  
  result <- getAnnotationFiles(assembly = "hg19", verbose = FALSE)
  
  expect_type(result, "list")
  expect_true("hg19" %in% names(result))
  expect_type(result$hg19, "character")
  expect_true(nchar(result$hg19[1]) > 0)
})


test_that("getAnnotationFiles returns valid structure for hg38", {
  skip_on_cran()
  skip_if_offline()
  
  result <- getAnnotationFiles(assembly = "hg38", verbose = FALSE)
  
  expect_type(result, "list")
  expect_true("hg38" %in% names(result))
  expect_type(result$hg38, "character")
})


skip_if_offline <- function() {
  has_internet <- tryCatch({
    con <- url("https://www.google.com", open = "rb", timeout = 2)
    close(con)
    TRUE
  }, error = function(e) FALSE)
  
  if (!has_internet) {
    skip("No internet connection")
  }
}