# ==============================================================================
# TESTS FOR buildInput()
# ==============================================================================

test_that("buildInput validates required inputs", {
  # Missing gene list
  expect_error(
    buildInput(
      bam_list = create_mock_bam_list(),
      genome = "hg19"
    ),
    "gene.*required|must.*supply"
  )
  
  # Missing BAM list
  expect_error(
    buildInput(
      gene_list = create_mock_gene_list(),
      genome = "hg19"
    ),
    "bam.*required|must.*supply"
  )
})

test_that("buildInput validates genome parameter", {
  gene_file <- create_mock_gene_list()
  bam_file <- create_mock_bam_list()
  
  expect_error(
    buildInput(
      gene_list = gene_file,
      bam_list = bam_file,
      genome = "mm10"
    ),
    "hg19|hg38"
  )
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput validates file types correctly", {
  gene_file <- create_mock_gene_list()
  bam_file <- create_mock_bam_list()
  
  # Invalid type_coverage
  expect_error(
    buildInput(
      gene_list = gene_file,
      bam_list = bam_file,
      genome = "hg19",
      type_coverage = "invalid"
    ),
    "bam|bed"
  )
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput fails with non-existent gene file", {
  bam_file <- create_mock_bam_list()
  
  expect_error(
    buildInput(
      gene_list = "/fake/path/genes.txt",
      bam_list = bam_file,
      genome = "hg19"
    ),
    "not found|does not exist"
  )
  
  cleanup_test_files(bam_file)
})

test_that("buildInput fails with non-existent BAM file", {
  gene_file <- create_mock_gene_list()
  
  expect_error(
    buildInput(
      gene_list = gene_file,
      bam_list = "/fake/path/bams.list",
      genome = "hg19"
    ),
    "not found|does not exist"
  )
  
  cleanup_test_files(gene_file)
})

test_that("buildInput handles invalid gene names gracefully", {
  # Create gene list with invalid genes
  invalid_genes <- c("FAKEGENE1", "NOTREALGENE", "INVALIDGENE")
  gene_file <- create_mock_gene_list(invalid_genes)
  bam_file <- create_mock_bam_list()
  
  # Should warn about invalid genes
  expect_warning(
    buildInput(
      gene_list = gene_file,
      bam_list = bam_file,
      genome = "hg19",
      type_coverage = "bed"
    ),
    "not found|invalid|skip"
  )
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput handles mix of valid and invalid genes", {
  # Mix of real and fake genes
  mixed_genes <- c("BRCA1", "FAKEGENE", "TP53")
  gene_file <- create_mock_gene_list(mixed_genes)
  bam_file <- create_mock_bam_list()
  
  # Should warn but continue with valid genes
  expect_warning(
    buildInput(
      gene_list = gene_file,
      bam_list = bam_file,
      genome = "hg19",
      type_coverage = "bed"
    ),
    "FAKEGENE"
  )
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput respects MAPQ and base quality filters", {
  skip("Requires real BAM files")
  
  gene_file <- create_mock_gene_list(c("BRCA1"))
  bam_file <- create_mock_bam_list()
  
  # Test different MAPQ values
  result1 <- buildInput(
    gene_list = gene_file,
    bam_list = bam_file,
    genome = "hg19",
    type_coverage = "bam",
    mapq = 10
  )
  
  result2 <- buildInput(
    gene_list = gene_file,
    bam_list = bam_file,
    genome = "hg19",
    type_coverage = "bam",
    mapq = 30
  )
  
  # Higher MAPQ should result in different coverage
  expect_false(identical(result1, result2))
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput handles coordinate system conversion", {
  skip("Requires real files")
  
  gene_file <- create_mock_gene_list(c("BRCA1"))
  bam_file <- create_mock_bam_list()
  
  # 0-based (BED standard)
  result_0based <- buildInput(
    gene_list = gene_file,
    bam_list = bam_file,
    genome = "hg19",
    type_coverage = "bed",
    coordinate_system = "0-based"
  )
  
  # 1-based (IGV, VCF)
  result_1based <- buildInput(
    gene_list = gene_file,
    bam_list = bam_file,
    genome = "hg19",
    type_coverage = "bed",
    coordinate_system = "1-based"
  )
  
  # Start positions should differ by 1
  expect_equal(
    result_0based$start + 1,
    result_1based$start
  )
  
  cleanup_test_files(c(gene_file, bam_file))
})

test_that("buildInput creates statistical summary", {
  skip("Requires real files")
  
  gene_file <- create_mock_gene_list(c("BRCA1", "TP53"))
  bam_file <- create_mock_bam_list()
  
  result <- buildInput(
    gene_list = gene_file,
    bam_list = bam_file,
    genome = "hg19",
    type_coverage = "bed"
  )
  
  # Should include statistical summary columns
  expected_stat_cols <- c(
    "SYMBOL",
    "Total_bases",
    "Mean_coverage",
    "Median_coverage",
    "bases_under_20x",
    "percentage_bases_under_20x"
  )
  
  # Check if any stat columns present (exact structure depends on implementation)
  expect_true(any(expected_stat_cols %in% colnames(result)))
  
  cleanup_test_files(c(gene_file, bam_file))
})