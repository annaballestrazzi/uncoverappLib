
# ============================================================================
# REACTIVE BASE: Calcolo una sola volta
# ============================================================================

# ============================================================================
# compute-annotation.R - ANNOTAZIONE VARIANTI
# ============================================================================

# ============================================================================
# STEP 1: QUERY ANNOTATION DATABASE (calcola UNA SOLA VOLTA)
# ============================================================================

intBED <- eventReactive(input$calc_annotations, {
  cat("\n=== ANNOTATION: intBED (TRIGGERED BY BUTTON) ===\n")
  
  if (is.null(filtered_low_nucl())) {
    cat("ERROR: No filtered data\n")
    return(NULL)
  }
  
  bedA <- filtered_low_nucl()
  cat("Low coverage positions:", nrow(bedA), "\n")
  
  if (nrow(bedA) == 0) {
    cat("ERROR: No positions to annotate\n")
    return(NULL)
  }
  
  # FILTER BY USER REGION (optional)
  query <- input$query_Database
  
  if (!is.null(query) && query != "") {
    cat("Filtering by region:", query, "\n")
    
    query.regions <- tryCatch({
      read.table(text=gsub("[:-]+", " ", query, perl=TRUE),
                 header=FALSE, col.names = c("chr", "start", "end"))
    }, error = function(e) NULL)
    
    if (!is.null(query.regions)) {
      query_chr <- if (grepl("^chr", query.regions$chr)) {
        query.regions$chr
      } else {
        paste0("chr", query.regions$chr)
      }
      
      bedA <- bedA %>%
        dplyr::filter(
          chromosome == query_chr,
          end >= query.regions$start,
          start <= query.regions$end
        )
      
      if (nrow(bedA) == 0) {
        cat("WARNING: No data in region\n")
        return(NULL)
      }
    }
  }
  
  # SELECT ANNOTATION FILE
  env_hg19 <- Sys.getenv("UNCOVERAPP_HG19_ANNOTATION", unset = "")
  env_hg38 <- Sys.getenv("UNCOVERAPP_HG38_ANNOTATION", unset = "")
  
  file.name <- NULL
  
  # Priority 1: Environment variables (for production deployment)
  if (identical(input$UCSC_Genome, "hg19") && nzchar(env_hg19) && file.exists(env_hg19)) {
    file.name <- env_hg19
    cat("Using hg19 from environment variable\n")
  } else if (identical(input$UCSC_Genome, "hg38") && nzchar(env_hg38) && file.exists(env_hg38)) {
    file.name <- env_hg38
    cat("Using hg38 from environment variable\n")
  } else {
    # Priority 2: Download via getAnnotationFiles()
    cat("Environment variables not set, downloading annotation files...\n")
    
    # Call getAnnotationFiles with assembly parameter
    m <- uncoverappLib::getAnnotationFiles(assembly = input$UCSC_Genome, verbose = TRUE)
    
    # Extract files for requested assembly
    if (is.list(m) && input$UCSC_Genome %in% names(m)) {
      data_files <- m[[input$UCSC_Genome]]
      
      # Filter for .gz file (not .tbi)
      gz_files <- data_files[grepl("\\.gz$", data_files) & !grepl("\\.tbi$", data_files)]
      
      if (length(gz_files) > 0) {
        file.name <- gz_files[1]
        cat("Using annotation file:", file.name, "\n")
      } else {
        stop("ERROR: No .gz annotation file found for ", input$UCSC_Genome)
      }
    } else {
      stop("ERROR: Cannot retrieve annotation files for ", input$UCSC_Genome)
    }
  }
  
  if (is.null(file.name) || !file.exists(file.name)) {
    stop("ERROR: Annotation file not found: ", file.name)
  }
  
  cat("Final annotation file:", file.name, "\n")
  
  # QUERY TABIX
  cat("Querying Tabix...\n")
  
  bedA_for_tabix <- bedA %>%
    dplyr::mutate(chromosome = sub("^chr", "", chromosome))
  
  bedA_gr <- tryCatch({
    GenomicRanges::makeGRangesFromDataFrame(
      bedA_for_tabix,
      seqnames.field = "chromosome",
      start.field = "start",
      end.field = "end",
      keep.extra.columns = FALSE,
      ignore.strand = TRUE
    )
  }, error = function(e) {
    cat("ERROR creating GRanges:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(bedA_gr)) return(NULL)
  
  query_gr <- if (length(bedA_gr) > 1000) {
    cat("Optimizing: merging intervals...\n")
    reduced <- GenomicRanges::reduce(bedA_gr, min.gapwidth = 100)
    cat("Reduced from", length(bedA_gr), "to", length(reduced), "\n")
    reduced
  } else {
    bedA_gr
  }
  
  result <- try(Rsamtools::scanTabix(file.name, param = query_gr), silent = FALSE)
  
  if (inherits(result, "try-error")) {
    cat("ERROR in Tabix\n")
    return(NULL)
  }
  
  lengths_result <- sapply(result, length)
  cat("Variants found:", sum(lengths_result), "\n")
  
  if (sum(lengths_result) == 0) {
    cat("No variants\n")
    return(NULL)
  }
  
  dff <- lapply(result, function(elt) {
    if (length(elt) == 0) return(data.frame())
    read.csv(textConnection(elt), sep="\t", header=FALSE, stringsAsFactors = FALSE)
  })
  
  valid_dfs <- dff[sapply(dff, nrow) > 0]
  if (length(valid_dfs) == 0) return(NULL)
  
  bedB <- do.call(rbind, valid_dfs)
  cat("Combined annotation:", nrow(bedB), "variants\n")
  
  # PROCESS ANNOTATION
  ncols <- ncol(bedB)
  
  if (ncols == 19) {
    colnames(bedB) <- c('Chromo', 'start','end','REF','ALT',
                        'dbsnp','GENENAME', 'PROTEIN_ensembl',
                        'MutationAssessor','SIFT','Polyphen2',
                        'M_CAP','CADD_PHED','AF_gnomAD','ClinVar',
                        'clinvar_MedGen_id','clinvar_OMIM_id','HGVSc_VEP','HGVSp_VEP')
  } else {
    colnames(bedB) <- paste0("V", 1:ncols)
  }
  
  bedB$Chromosome <- paste0("chr", bedB[[1]])
  bedB <- bedB[, -1]
  bedB$Chromosome <- as.character(bedB$Chromosome)
  bedB$AF_gnomAD <- suppressWarnings(as.numeric(bedB$AF_gnomAD))
  bedB$CADD_PHED <- suppressWarnings(as.numeric(bedB$CADD_PHED))
  
  # OVERLAP
  cat("Computing overlap...\n")
  # Normalize chromosome naming: remove "chr" prefix from both datasets
  bedA_norm <- bedA %>%
    dplyr::mutate(chromosome = sub("^chr", "", chromosome))

  bedB_norm <- bedB %>%
    dplyr::mutate(Chromosome = sub("^chr", "", Chromosome))

  cat("bedA chromosomes:", unique(bedA_norm$chromosome), "\n")
  cat("bedB chromosomes:", unique(bedB_norm$Chromosome), "\n")
  bed1_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedA_norm, seqnames.field = "chromosome", ignore.strand = TRUE, keep.extra.columns = TRUE
  )
  
  bed2_gr <- GenomicRanges::makeGRangesFromDataFrame(
    bedB_norm, seqnames.field = "Chromosome", ignore.strand = TRUE, keep.extra.columns = TRUE
  )
  
  tp <- GenomicRanges::findOverlaps(query = bed2_gr, subject = bed1_gr, type = "any")
  
  cat("Overlaps:", length(tp), "\n")
  
  if (length(tp) == 0) return(data.frame())
  
  hits_bed2_idx <- S4Vectors::queryHits(tp)
  hits_bed1_idx <- S4Vectors::subjectHits(tp)
  
  bed1_hits <- bedA[hits_bed1_idx, ]
  bed2_hits <- bedB[hits_bed2_idx, ]
  
  intersect_df <- data.frame(
    seqnames = bed2_hits$Chromosome,
    start = bed2_hits$start,
    end = bed2_hits$end,
    coverage = bed1_hits$coverage,
    REF = bed2_hits$REF,
    ALT = bed2_hits$ALT,
    dbsnp = bed2_hits$dbsnp,
    GENENAME = bed2_hits$GENENAME,
    PROTEIN_ensembl = bed2_hits$PROTEIN_ensembl,
    MutationAssessor = bed2_hits$MutationAssessor,
    SIFT = bed2_hits$SIFT,
    Polyphen2 = bed2_hits$Polyphen2,
    M_CAP = bed2_hits$M_CAP,
    CADD_PHED = bed2_hits$CADD_PHED,
    AF_gnomAD = bed2_hits$AF_gnomAD,
    ClinVar = bed2_hits$ClinVar,
    clinvar_MedGen_id = bed2_hits$clinvar_MedGen_id,
    HGVSc_VEP = bed2_hits$HGVSc_VEP,
    HGVSp_VEP = bed2_hits$HGVSp_VEP,
    stringsAsFactors = FALSE
  )
  
  cat("Final result:", nrow(intersect_df), "variants\n")
  return(intersect_df)
}, ignoreNULL = TRUE, ignoreInit = TRUE)


# ============================================================================
# STEP 2: CLEAN DATAFRAME (per operazioni dplyr, export, ecc.)
# ============================================================================

annotated_variants_data <- reactive({
  cat("\n=== ANNOTATION: Preparing clean dataframe ===\n")
  
  data <- intBED()
  
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  # Aggiungi colonna helper per highlight
  data$highlight_important <- grepl("H|M", data$MutationAssessor) & 
                              data$ClinVar != "." & 
                              !is.na(data$AF_gnomAD) & 
                              data$AF_gnomAD < 0.5
  
  cat("Dataframe ready:", nrow(data), "variants\n")
  return(data)
})

# ============================================================================
# STEP 3: DT WIDGET (SOLO per visualizzazione)
# ============================================================================

condform_table <- reactive({
  cat("\n=== ANNOTATION: Building DT widget ===\n")
  
  progress <- shiny::Progress$new()
  on.exit(progress$close())
  progress$set(message = "Formatting table...", value = 0.5)
  
  data <- annotated_variants_data()
  
  if (is.null(data) || nrow(data) == 0) {
    return(NULL)
  }
  
  cat("Formatting", nrow(data), "rows\n")
  
  dt <- DT::datatable(
    data,
    options = list(
      pageLength = 25,
      scrollX = TRUE,
      scrollY = "600px",
      dom = 'Bfrtip',
      buttons = c('copy', 'csv', 'excel'),
      columnDefs = list(
        list(targets = which(colnames(data) == "highlight_important") - 1, visible = FALSE)
      )
    ),
    extensions = 'Buttons',
    rownames = FALSE
  ) %>%
    DT::formatStyle('ClinVar',
                    backgroundColor = DT::styleEqual(c('.'), c('lightgreen'), default = 'lightcoral')) %>%
    DT::formatStyle('CADD_PHED',
                    backgroundColor = DT::styleInterval(c(20), c('lightgreen', 'lightcoral'))) %>%
    DT::formatStyle('MutationAssessor',
                    backgroundColor = DT::styleEqual(c('H', 'M'), c('lightcoral', 'lightyellow'), default = 'lightgreen')) %>%
    DT::formatStyle('M_CAP',
                    backgroundColor = DT::styleEqual(c('D'), c('lightcoral'), default = 'lightgreen')) %>%
    DT::formatStyle('AF_gnomAD',
                    backgroundColor = DT::styleInterval(c(0.5), c('lightcoral', 'lightgreen'))) %>%
    DT::formatStyle(c('start', 'end'), 'highlight_important',
                    backgroundColor = DT::styleEqual(c(TRUE, FALSE), c('yellow', 'white'))) %>%
    DT::formatStyle(c('start', 'end'), 'highlight_important',
                    color = DT::styleEqual(c(TRUE, FALSE), c('red', 'green')))
  
  cat("DT widget created\n")
  return(dt)
})