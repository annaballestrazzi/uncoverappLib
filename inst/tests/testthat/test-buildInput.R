# tests/testthat/test-buildInput.R
# Test essenziali per la funzione buildInput()

test_that("buildInput rejects invalid genome", {
  temp_genes <- tempfile(fileext = ".txt")
  temp_bams <- tempfile(fileext = ".list")
  writeLines("BRCA1", temp_genes)
  writeLines("/fake/path.bam", temp_bams)
  
  expect_error(
    buildInput(
      geneList = temp_genes,
      bamList = temp_bams,
      genome = "hg99",
      type_bam = "chr",
      type_input = "genes",
      outDir = tempdir()
    )
  )
  
  unlink(c(temp_genes, temp_bams))
})


test_that("buildInput rejects invalid MAPQ", {
  temp_genes <- tempfile(fileext = ".txt")
  temp_bams <- tempfile(fileext = ".list")
  writeLines("BRCA1", temp_genes)
  writeLines("/fake/path.bam", temp_bams)
  
  expect_error(
    buildInput(
      geneList = temp_genes,
      bamList = temp_bams,
      genome = "hg19",
      type_bam = "chr",
      type_input = "genes",
      outDir = tempdir(),
      MAPQ.min = -5
    ),
    "MAPQ.min must be greater than 0"
  )
  
  unlink(c(temp_genes, temp_bams))
})


test_that("buildInput validates type_coverage", {
  temp_genes <- tempfile(fileext = ".txt")
  temp_bams <- tempfile(fileext = ".list")
  writeLines("BRCA1", temp_genes)
  writeLines("/fake/path.bam", temp_bams)
  
  expect_error(
    buildInput(
      geneList = temp_genes,
      bamList = temp_bams,
      genome = "hg19",
      type_bam = "chr",
      type_input = "genes",
      type_coverage = "invalid",
      outDir = tempdir()
    ),
    "type_coverage must be 'bam' or 'bed'"
  )
  
  unlink(c(temp_genes, temp_bams))
})