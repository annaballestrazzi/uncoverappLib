utils::globalVariables(c(
  "ENTREZID", "ALIAS", "SYMBOL", "chr", "count",
  "seqnames.1", "start.1", "end.1", "width.1", "strand.1", "value"
))
#' Build input file for coverage analysis
#' @importFrom stats median
#' @importFrom utils download.file head read.csv
#' 
#' @description
#' Function to build input file for coverage analysis with optimized processing
#' and support for both BAM and BED coverage files.
#'
#' @details
#' **Coverage File Types:**
#' 
#' - **BAM files** (`type_coverage = "bam"`): 
#'   - Requires `MAPQ.min` and `base.quality` parameters
#'   - Automatically uses 1-based coordinates
#'   - Generates per-base coverage via samtools pileup
#'   
#' - **BED coverage files** (`type_coverage = "bed"`):
#'   - Requires `input_coord_system` parameter ("0-based" or "1-based")
#'   - `MAPQ.min` and `base.quality` are ignored
#'   - Faster processing for pre-computed coverage
#'   - Expected format: chromosome, start, end, coverage_depth
#'
#' @return Two files: 
#' 1. A .bed file containing tab-separated genomic coordinates (chromosome, start, end)
#'    with coverage values for each position
#' 2. A .txt file with statistical summary of coverage per gene and sample
#'
#' @param geneList A text file (.txt) containing HGNC official gene name(s) one per row,
#'   or a BED file with genomic coordinates (chr, start, end, SYMBOL)
#' @param genome (char) Reference genome: "hg19" or "hg38"
#' @param chromosome_notation (char) Chromosome notation: "chr" (chr1,chr2,...,chrX,chrM) or 
#'   "number" (1,2,...,X,M)
#' @param sampleList A text file (.list) containing absolute paths to coverage files 
#'   (BAM or BED) one per row
#' @param outDir (char) Directory where output files will be stored
#' @param type_input (char) Type of input target: "genes" (gene list) or "target" (BED file)
#' @param type_coverage (char) Type of coverage files: "bam" or "bed" (default: "bam")
#' @param MAPQ.min (integer) **[BAM only]** Minimum mapping quality for alignments (default: 1). 
#'   Ignored when `type_coverage = "bed"`
#' @param base.quality (integer) **[BAM only]** Minimum base quality for nucleotides (default: 1).
#'   Ignored when `type_coverage = "bed"`
#' @param input_coord_system (char) **[BED only]** Coordinate system of input BED coverage files: 
#'   "0-based" (BED standard) or "1-based" (IGV/VCF). Default: "1-based". 
#'   Ignored when `type_coverage = "bam"` (BAM is always 1-based)
#' @param annotation_file (char) Optional path to annotation file (tab-separated 
#'   with SYMBOL column for adding gene-level annotations)
#'
#' @export
#' @import Rsamtools
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import GenomicRanges
#' @import OrganismDbi
#' @import org.Hs.eg.db
#' @importFrom utils write.table read.table
#' @importFrom S4Vectors queryHits subjectHits
#' @importFrom dplyr filter group_by summarise mutate select inner_join full_join
#' @importFrom tidyselect all_of
#'
#' @examples
#' # Show function signature
#' args(buildInput)
#' 
#' \dontrun{
#' # Process BAM files
#' buildInput(
#'   geneList = "genes.txt",
#'   sampleList = "samples.list",
#'   genome = "hg38",
#'   chromosome_notation = "chr",
#'   type_input = "genes",
#'   type_coverage = "bam",
#'   outDir = tempdir()
#' )
#' }
#' \dontrun{
#' # Example 1: Process BAM files
#' gene_list <- system.file("extdata", "mygene.txt", package = "uncoverappLib")
#' bam_example <- system.file("extdata", "example_POLG.bam", package = "uncoverappLib")
#' cat(bam_example, file = "bam.list", sep = "\n")
#' 
#' buildInput(
#'   geneList = gene_list, 
#'   genome = "hg19", 
#'   chromosome_notation = "chr",
#'   sampleList = "bam.list",
#'   type_input = "genes", 
#'   outDir = tempdir(),
#'   type_coverage = "bam",
#'   MAPQ.min = 20,
#'   base.quality = 20
#' )
#' 
#' # Example 2: Process BED coverage files
#' buildInput(
#'   geneList = gene_list,
#'   genome = "hg38",
#'   chromosome_notation = "chr",
#'   sampleList = "coverage.list",
#'   type_input = "genes",
#'   outDir = tempdir(),
#'   type_coverage = "bed",
#'   input_coord_system = "0-based"
#' )
#' }



buildInput <- function(geneList, 
                       genome, 
                       chromosome_notation, 
                       sampleList, 
                       outDir, 
                       type_input,
                       MAPQ.min = 1, 
                       base.quality = 1, 
                       type_coverage = "bam",
                       input_coord_system = "1-based",
                       annotation_file = NULL) {
                        script_start <- Sys.time()

  
  # ==============================================================================
  # INPUT VALIDATION
  # ==============================================================================
  
  if (missing(geneList)) stop("geneList must be supplied.\n")
  if (missing(sampleList)) stop("sampleList must be supplied.\n")
  if (MAPQ.min <= 0) stop("MAPQ.min must be greater than 0")
  if (base.quality <= 0) stop("base.quality must be greater than 0")
  if (!type_coverage %in% c("bam", "bed")) {
    stop("type_coverage must be 'bam' or 'bed'")
  }
  if (!input_coord_system %in% c("0-based", "1-based")) {
    stop("input_coord_system must be '0-based' or '1-based'")
  }
  
  # ==============================================================================
  # INITIALIZATION
  # ==============================================================================

  message("=== COVERAGE ANALYSIS STARTED ===")
  message(Sys.time())
  message("\n=== CONFIGURATION ===")
  message("Gene file: ", geneList)
  message("Reference genome: ", genome)
  message("Chromosome notation: ", chromosome_notation)
  message("Sample list: ", sampleList)
  message("Output directory: ", outDir)
  message("Input type: ", type_input)
  message("Coverage type: ", type_coverage)
  message("Coordinate system: ", input_coord_system)
  # Conditional parameters (only for BAM)
  if (type_coverage == "bam") {
    message("MAPQ minimum: ", MAPQ.min)
    message("Base quality minimum: ", base.quality)
  }
  # Optional annotation
  if (!is.null(annotation_file)) {
    message("Annotation file: ", annotation_file)
  }
  message("=============================\n")
  # ==============================================================================
  # STEP 1: LOAD INPUT FILES
  # ==============================================================================
  
  message("\n=== STEP 1: LOADING INPUT FILES ===")
  start_time <- Sys.time()
  
  # Load gene list or BED file
  if (type_input == "target") {
    message("Loading target BED file...")
    gene.List.r <- utils::read.table(geneList, stringsAsFactors = FALSE)
    gene.List <- gene.List.r[1:4]
    colnames(gene.List) <- c('chr', 'start', 'end', 'SYMBOL')
    message(paste("Loaded", nrow(gene.List), "regions from BED file"))
  } else {
    message("Loading gene list...")
    gene.List <- base::scan(geneList, character(), quote = "")
    message(paste("Loaded", length(gene.List), "genes"))
  }
  
  # Load coverage file list
  list_coverage <- base::scan(sampleList, character(), quote = "")
  message(paste("Loaded", length(list_coverage), "coverage files"))
  
  # Select genome
  if (genome == "hg19") {
    all_gene <- TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    all_gene <- TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  message(paste("Using genome:", genome))
  
  message(paste("Input loading took:", 
                round(difftime(Sys.time(), start_time, units="secs"), 2), "seconds\n"))
  
  # ==============================================================================
  # CREATE OUTPUT DIRECTORY EARLY (for log files)
  # ==============================================================================
  
  dir_users <- sprintf("%s/output/", outDir)
  myDir <- dir_users
  if (!file.exists(myDir)) {
    dir.create(myDir, recursive = TRUE, showWarnings = FALSE)
  }
  message(paste("Output directory:", myDir, "\n"))
  
  # ==============================================================================
  # STEP 2: GENE VALIDATION (for gene list input)
  # ==============================================================================
  
  validated_genes <- NULL
  
  if (type_input != "target") {
    message("=== STEP 2: GENE VALIDATION ===")
    start_time <- Sys.time()
    
    # Stage 1: Convert to ENTREZ IDs
    message("Stage 1: Converting gene names to ENTREZ IDs...")
    my_gene_name <- OrganismDbi::select(
      org.Hs.eg.db, 
      key = gene.List, 
      columns = "ENTREZID", 
      keytype = "ALIAS"
    ) %>%
      as.data.frame() %>%
      dplyr::filter(!is.na(ENTREZID)) %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()
    
    found_genes <- unique(my_gene_name$ALIAS)
    not_found <- setdiff(gene.List, found_genes)
    
    if (length(not_found) > 0) {
      warning(paste("Genes not found in org.Hs.eg.db:", 
                    paste(not_found, collapse=", ")))
      log_file <- file.path(myDir, 'preprocessing_log.txt')
      utils::write.table(not_found, file = log_file, 
                        quote = FALSE, row.names = FALSE, col.names = FALSE)
      message(paste("List of unrecognized genes written to:", log_file))
    }
    
    if (nrow(my_gene_name) == 0) {
      stop("ERROR: No valid genes found in org.Hs.eg.db.")
    }
    
    message(paste("Stage 1:", nrow(my_gene_name), "mappings from", 
                  length(found_genes), "genes"))
    
    # Stage 2: Validate against genome
    message("Stage 2: Validating against reference genome...")
    all_entrez_ids <- unique(my_gene_name$ENTREZID)
    genome_keys <- AnnotationDbi::keys(all_gene, "GENEID")
    valid_entrez_ids <- intersect(all_entrez_ids, genome_keys)
    invalid_entrez_ids <- setdiff(all_entrez_ids, genome_keys)
    
    if (length(invalid_entrez_ids) > 0) {
      invalid_genes <- my_gene_name %>% 
        dplyr::filter(ENTREZID %in% invalid_entrez_ids) %>%
        dplyr::pull(ALIAS) %>%
        unique()
      
      warning(paste("Genes not in reference genome:", 
                    paste(invalid_genes, collapse=", ")))
      log_file <- file.path(myDir, 'preprocessing_log1.txt')
      utils::write.table(invalid_genes, 
                        file = log_file, 
                        quote = FALSE, row.names = FALSE, col.names = FALSE)
      message(paste("List of genes not in genome written to:", log_file))
    }
    
    my_gene_name <- my_gene_name %>%
      dplyr::filter(ENTREZID %in% valid_entrez_ids)
    
    if (length(valid_entrez_ids) == 0) {
      stop("ERROR: No genes passed validation.")
    }
    
    message(paste("Stage 2:", length(valid_entrez_ids), "validated genes"))
    message(paste("Gene validation took:", 
                  round(difftime(Sys.time(), start_time, units="secs"), 2), 
                  "seconds\n"))
    
    validated_genes <- list(
      my_gene_name = my_gene_name,
      valid_entrez_ids = valid_entrez_ids
    )
  }
  
  # ==============================================================================
  # STEP 3: CREATE TARGET BED
  # ==============================================================================
  
  message("=== STEP 3: CREATING TARGET BED ===")
  start_time <- Sys.time()
  
  if (type_input == "target") {
    for_bed <- gene.List
  } else {
    # Extract exon coordinates
    my_gene_name <- validated_genes$my_gene_name
    ID <- validated_genes$valid_entrez_ids
    
    cds <- OrganismDbi::exonsBy(all_gene, by = "tx", use.names = TRUE)
    txid_map <- OrganismDbi::select(all_gene, ID, "TXNAME", "GENEID")
    
    # Extract coordinates for each gene (optimized)
    pre_list <- lapply(split(txid_map, txid_map$GENEID), function(gene_txs) {
      tx_names <- gene_txs$TXNAME
      coor <- as.data.frame(cds[names(cds) %in% tx_names])
      if (nrow(coor) == 0) return(NULL)
      coor$ENTREZID <- gene_txs$GENEID[1]
      return(coor)
    })
    
    pre_list <- pre_list[!sapply(pre_list, is.null)]
    
    if (length(pre_list) == 0) {
      stop("No valid transcript coordinates found!")
    }
    
    # Combine and format
    pre <- dplyr::bind_rows(pre_list) %>%
      dplyr::inner_join(my_gene_name, by = "ENTREZID", relationship = "many-to-many") %>%
      dplyr::mutate(
        start = start - 10,
        end = end + 10
      )
    
    for_bed <- pre %>%
      dplyr::transmute(
        chr = as.character(seqnames),
        start = start,
        end = end,
        SYMBOL = ALIAS
      ) %>%
      dplyr::distinct()
    
    # Remove alternative contigs
    for_bed <- for_bed[!grepl("_", for_bed$chr), ]
    
    # Apply notation
    if (chromosome_notation == "number") {
      for_bed <- for_bed %>%
        dplyr::mutate(chr = sub("^.{0,3}", "", chr, perl = TRUE))
    }
  }
  
  message(paste("Created BED with", nrow(for_bed), "intervals"))
  message(paste("BED creation took:", 
                round(difftime(Sys.time(), start_time, units="secs"), 2), 
                "seconds\n"))
  
  # ==============================================================================
  # STEP 4: PROCESS COVERAGE DATA
  # ==============================================================================
  
  message("=== STEP 4: PROCESSING COVERAGE ===")
  coverage_time <- Sys.time()
  
  message(paste("Input coordinate system:", input_coord_system))
  
  # Convert target BED to GRanges (1-based)
  target_bed_1based <- for_bed
  if (input_coord_system == "0-based") {
    target_bed_1based <- target_bed_1based %>% dplyr::mutate(start = start + 1)
  }
  
  for_grange <- GenomicRanges::makeGRangesFromDataFrame(
    target_bed_1based, 
    keep.extra.columns = TRUE
  )
  
if (type_coverage == "bam") {
    message("Processing BAM files with pileup...")
    
    # Check chromosome availability
    message("Checking BAM chromosome availability...")
    idx_stats <- system(paste("samtools idxstats", list_coverage[1]), intern = TRUE)
    idx_df <- read.table(text = idx_stats, header = FALSE, 
                         col.names = c("seqnames", "length", "mapped", "unmapped"),
                         stringsAsFactors = FALSE)
    chrs_with_data <- idx_df$seqnames[idx_df$mapped > 0]
    target_chrs <- as.character(unique(GenomicRanges::seqnames(for_grange)))
    available_chrs <- intersect(target_chrs, chrs_with_data)
    
    message(paste("Chromosomes in BAM:", nrow(idx_df)))
    message(paste("Chromosomes with data:", length(chrs_with_data)))
    message(paste("Target chromosomes:", length(target_chrs)))
    message(paste("Available overlap:", length(available_chrs)))
    
    if (length(available_chrs) == 0) {
      stop("No chromosome overlap between BAM and target regions!")
    }
    
    # Filter target regions
    for_grange_filtered <- for_grange[
      as.character(GenomicRanges::seqnames(for_grange)) %in% available_chrs
    ]
    
    param <- Rsamtools::ScanBamParam(which = for_grange_filtered)
    p_param <- Rsamtools::PileupParam(
      distinguish_nucleotides = FALSE,
      distinguish_strands = FALSE,
      min_mapq = MAPQ.min,
      min_base_quality = base.quality,
      min_nucleotide_depth = 1,
      max_depth = 150000
    )
    
    lst3 <- lapply(list_coverage, function(i) {
      message(paste("  Processing:", basename(i)))
      pileup_df <- Rsamtools::pileup(i, scanBamParam = param, pileupParam = p_param)
      pileup_df <- pileup_df[, -5]
      pileup_df <- pileup_df[!duplicated(pileup_df), ]
      
      pileup_df %>%
        dplyr::mutate(end = pos) %>%
        dplyr::group_by(seqnames, pos, end) %>%
        dplyr::summarise(count = sum(count), .groups = "drop")
    })
    
    pp <- Reduce(function(...) merge(..., by = c("seqnames", "pos", "end")), lst3)
    pp[is.na(pp)] <- 0
    colnames(pp)[1:3] <- c("seqnames", "start", "end")
    
    # Add SYMBOL
    message("Adding gene annotations to BAM coverage...")
    pp_gr <- GenomicRanges::makeGRangesFromDataFrame(pp, keep.extra.columns = TRUE)
    overlaps <- GenomicRanges::findOverlaps(pp_gr, for_grange_filtered, type = "any")
    pp_gr_matched <- pp_gr[queryHits(overlaps)]
    mcols(pp_gr_matched)$SYMBOL <- for_grange_filtered[subjectHits(overlaps)]$SYMBOL
    pp <- as.data.frame(pp_gr_matched) %>%
      dplyr::select(-width, -strand) %>%
      dplyr::distinct()
    
    message(paste("Added SYMBOL to", nrow(pp), "positions"))
    message("BAM pileup: per-base coverage (1-based positions)")
    
  } else if (type_coverage == "bed") {
    message("Processing BED coverage files...")
    
    df_list <- lapply(seq_along(list_coverage), function(idx) {
      file <- list_coverage[idx]
      message(paste("Sample", idx, ":", basename(file)))
      
      df <- utils::read.table(
        file, 
        header = FALSE,
        stringsAsFactors = FALSE,
        colClasses = c("character", "integer", "integer", "numeric")
      )
      colnames(df) <- c("seqnames", "start", "end", "count")
      
      message(paste("  Raw intervals:", nrow(df)))
      
      # Convert to 1-based
      if (input_coord_system == "0-based") {
        df <- df %>% dplyr::mutate(start = start + 1)
      }
      
      df$interval_id <- paste0("interval_", seq_len(nrow(df)))
      
      df_gr <- GenomicRanges::makeGRangesFromDataFrame(df, keep.extra.columns = TRUE)
      overlaps <- GenomicRanges::findOverlaps(df_gr, for_grange, type = "any")
      
      if (length(overlaps) == 0) {
        warning(paste("No overlaps found for", basename(file)))
        return(data.frame(
          seqnames = character(), 
          start = integer(), 
          end = integer(), 
          count = numeric(),
          stringsAsFactors = FALSE
        ))
      }
      
      df_gr_overlapping <- df_gr[queryHits(overlaps)]
      for_grange_overlapping <- for_grange[subjectHits(overlaps)]
      
      df_gr_trimmed <- GenomicRanges::pintersect(
        df_gr_overlapping,
        for_grange_overlapping,
        ignore.strand = TRUE
      )
      
      df_trimmed <- as.data.frame(df_gr_trimmed)
      df_trimmed <- df_trimmed %>%
        dplyr::select(seqnames, start, end, count) %>%
        dplyr::distinct()
      
      message(paste("  Final intervals:", nrow(df_trimmed)))
      return(df_trimmed)
    })
    
    df_list <- df_list[sapply(df_list, nrow) > 0]
    
    if (length(df_list) == 0) {
      stop("ERROR: No coverage data overlaps with target regions!")
    }
    
    if (length(df_list) == 1) {
      pp <- df_list[[1]]
    } else {
      pp <- df_list[[1]]
      for (i in 2:length(df_list)) {
        pp <- merge(pp, df_list[[i]], by = c("seqnames", "start", "end"), 
                    all = TRUE, suffixes = c("", paste0(".sample", i)))
      }
    }
    
    pp[is.na(pp)] <- 0
  }
  
  # Metadata (OUTSIDE if/else)
  attr(pp, "coordinate_system") <- "1-based"
  attr(pp, "trimmed_to_targets") <- TRUE
  attr(pp, "original_input_system") <- input_coord_system
  
  # Normalize chromosomes (OUTSIDE if/else)
  if (chromosome_notation == "chr") {
    pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
  } else {
    pp$seqnames <- sub("^chr", "", as.character(pp$seqnames))
  }
  
  # Check missing genes (OUTSIDE if/else)
  if ("SYMBOL" %in% colnames(pp)) {
    target_genes <- unique(for_grange$SYMBOL)
    covered_genes <- unique(pp$SYMBOL)
    missing_genes <- setdiff(target_genes, covered_genes)
    
    if (length(missing_genes) > 0) {
      warning(paste(length(missing_genes), "genes have NO coverage!"))
      message("First 20 missing genes:", paste(head(missing_genes, 20), collapse=", "))
      
      log_file <- file.path(myDir, 'genes_without_coverage.txt')
      write.table(missing_genes, file = log_file, 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      message(paste("Full list written to:", log_file))
    } else {
      message("SUCCESS: All target genes have coverage data!")
    }
  }
  
  message(paste("Final coverage dimensions:", nrow(pp), "intervals x", ncol(pp)-3, "samples"))
  message(paste("Coverage processing took:", 
                round(difftime(Sys.time(), coverage_time, units = "mins"), 2), "minutes\n"))

  
  # ==============================================================================
  # STEP 5: CALCULATE STATISTICS
  # ==============================================================================
  
  message("=== STEP 5: CALCULATING STATISTICS ===")
  start_time <- Sys.time()
  
  ppinp <- as.data.frame(pp)
  colnames(ppinp)[1:3] <- c("seqnames", "start", "end")
  
  # Sample names
  n <- length(colnames(ppinp)[-1:-3])
  samples <- tools::file_path_sans_ext(basename(list_coverage))
  
  if (length(samples) != n) {
    samples <- paste0("sample_", seq_len(n))
  }
  
  colnames(ppinp)[-1:-3] <- paste0("sample_", samples)[seq_len(n)]
  
  # Create GRanges
  for_range_pp <- GenomicRanges::makeGRangesFromDataFrame(
    ppinp,
    keep.extra.columns = TRUE
  )
  
  # Find overlaps
  tp <- GenomicRanges::findOverlaps(
    query = for_range_pp,
    subject = for_grange,
    type = "any",
    select = "all"
  )
  
  # Merge coverage with gene annotations
  sts_df <- data.frame(
    for_range_pp[queryHits(tp), ],
    for_grange[subjectHits(tp), ]
  )
  
  statistiche <- sts_df[!duplicated(sts_df[, c("start", "end", "SYMBOL")]), ]
  statistiche <- subset(
    statistiche,
    select = -c(width, strand, seqnames.1, start.1, end.1, width.1, strand.1)
  )
  
  colnames(statistiche)[1:3] <- c("chromosome", "start", "end")
  
  merge_g <- dplyr::full_join(target_bed_1based, statistiche, by = "SYMBOL",
                              relationship = "many-to-many")
  
  col_name <- colnames(merge_g)
  col.sub <- col_name[grepl("sample_", col_name)]
  merge_g[col.sub] <- sapply(merge_g[col.sub], as.numeric)
  
  # Calculate interval widths (1-based: end - start + 1)
  merge_g$width <- merge_g$end.x - merge_g$start.x + 1
  
  # Prepare per-sample data
  x <- lapply(col.sub, function(i) {
    x_df <- merge_g %>%
      dplyr::select("chromosome", "start.x", "end.x", "SYMBOL", "width", 
                   tidyselect::all_of(i))
    colnames(x_df)[6] <- "value"
    x_df$sample <- i
    x_df
  })
  
  # Calculate statistics with proper weighting
  statistical_operation <- function(df) {
    df <- df %>% dplyr::filter(!is.na(value), !is.na(width))
    
    if (nrow(df) == 0) {
      return(data.frame(
        SYMBOL = character(),
        sample = character(),
        Total_bases = numeric(),
        Mean_coverage = numeric(),
        Median_coverage = numeric(),
        number_of_intervals_under_20x = numeric(),
        bases_under_20x = numeric(),
        percentage_bases_under_20x = numeric(),
        stringsAsFactors = FALSE
      ))
    }
    
    df %>%
      dplyr::group_by(SYMBOL, sample) %>%
      dplyr::summarize(
        Total_bases = sum(width, na.rm = TRUE),
        Mean_coverage = sum(value * width, na.rm = TRUE) / sum(width, na.rm = TRUE),
        Median_coverage = median(value, na.rm = TRUE),
        number_of_intervals_under_20x = sum(value < 20, na.rm = TRUE),
        bases_under_20x = sum(width[value < 20], na.rm = TRUE),
        percentage_bases_under_20x = (sum(width[value < 20], na.rm = TRUE) / 
                                      sum(width, na.rm = TRUE)) * 100,
        .groups = "drop"
      ) %>%
      as.data.frame() %>%
      dplyr::mutate_if(is.numeric, round, 3)
  }
  
  stat_summ <- do.call(rbind, lapply(x, statistical_operation))
  
  message(paste("Statistics computation took:",
                round(difftime(Sys.time(), start_time, units = "secs"), 2),
                "seconds\n"))
  
  # ==============================================================================
  # STEP 6: ADD OPTIONAL ANNOTATIONS
  # ==============================================================================
  
  if (!is.null(annotation_file)) {
    message("=== STEP 6: ADDING ANNOTATIONS ===")
    start_time <- Sys.time()
    
    message(paste("Loading annotations from:", annotation_file))
    
    if (!file.exists(annotation_file)) {
      warning(paste("Annotation file not found:", annotation_file))
    } else {
      tryCatch({
        omim_gene <- utils::read.table(
          annotation_file, 
          header = TRUE, 
          sep = "\t", 
          stringsAsFactors = FALSE,
          quote = ""
        )
        
        if (!"SYMBOL" %in% colnames(omim_gene)) {
          warning("Annotation file must contain a 'SYMBOL' column. Skipping annotations.")
        } else {
          message(paste("Loaded", nrow(omim_gene), "annotations"))
          
          stat_summ <- merge(
            stat_summ, 
            omim_gene, 
            by = "SYMBOL", 
            all.x = TRUE
          )
          
          message("Annotations added successfully")
          message(paste("Annotation took:", 
                       round(difftime(Sys.time(), start_time, units = "secs"), 2), 
                       "seconds\n"))
        }
      }, error = function(e) {
        warning(paste("Error loading annotation file:", e$message))
      })
    }
  }
  
  # ==============================================================================
  # STEP 7: WRITE OUTPUTS
  # ==============================================================================
  
  message("=== STEP 7: WRITING OUTPUTS ===")
  
  # Directory already created at the beginning
  # Generate timestamp for filenames
  timestamp <- format(Sys.time(), "%a_%b_%d_%Y")
  
  # Rename columns before writing
  pp_output <- pp

  if ("SYMBOL" %in% colnames(pp_output)) {
  message("Removing SYMBOL from coverage output (for compute-annotation.R)")
  pp_output <- pp_output %>% dplyr::select(-SYMBOL)
  }

  # 1. Rename first column from "seqnames" to "chromosome"
  colnames(pp_output)[1] <- "chromosome"

  # 2. Rename count columns to include filename prefix
  sample_names <- tools::file_path_sans_ext(basename(list_coverage))
  count_cols <- colnames(pp_output)[-(1:3)]  # All columns except chr, start, end

  # Create new column names with "count_" prefix
  new_count_names <- paste0("count_", sample_names)

  # Apply new names to count columns
  colnames(pp_output)[4:ncol(pp_output)] <- new_count_names
  # Write coverage file
  coverage_file <- file.path(myDir, paste0(timestamp, '.bed'))
  utils::write.table(
    x = pp_output, 
    file = coverage_file,
    quote = FALSE, 
    sep = "\t", 
    eol = "\n", 
    row.names = FALSE,
    col.names = TRUE
  )
  message(paste("Coverage data written to:", coverage_file))
  
  # Write statistics file
  stats_file <- file.path(myDir, paste0(timestamp, '_statistical_summary.txt'))
  utils::write.table(
    x = stat_summ,
    file = stats_file,
    quote = FALSE, 
    sep = "\t", 
    eol = "\n", 
    row.names = FALSE,
    col.names = TRUE, 
    dec = "."
  )
  message(paste("Statistical summary written to:", stats_file))
  
  message("\n=== ANALYSIS COMPLETE ===")
  message(Sys.time())
  message("Output directory: ", myDir)
  cat("Total script execution:",round(difftime(Sys.time(), script_start, units = "secs"), 2),"seconds\n")

  # Return results invisibly
  invisible(list(
    coverage = pp,
    statistics = stat_summ,
    output_dir = myDir
  ))
}