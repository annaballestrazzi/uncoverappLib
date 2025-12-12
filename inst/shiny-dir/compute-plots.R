# ============================================================================
# compute-plots.R - GENE COVERAGE PLOTS
# ============================================================================

# ============================================================================
# DATA TRACKS (Coverage visualization)
# ============================================================================
dtrack1 <- reactive({
  grcoverage <- filtered_low()
  dt1 <- Gviz::DataTrack(
    range = grcoverage,
    type = "histogram", 
    fill.histogram = "red", 
    col.histogram = NA,
    genome = input$UCSC_Genome,
    name = "Seq. Depth",
    cex.title = 1.2 
  )
  
  cat("DataTrack LOW created OK\n")
  return(dt1)
})

dtrackHigh <- reactive({
  grcoverage_high <- filtered_high()
  dt2 <- Gviz::DataTrack(
    range = grcoverage_high,
    type = "histogram", 
    fill.histogram = "dodgerblue",
    col.histogram = NA, 
    genome = input$UCSC_Genome,
    name = "Seq. Depth",
    cex.title = 1.2 
  )
  
  cat("DataTrack HIGH created OK\n")
  return(dt2)
})

# ============================================================================
# CHROMOSOME IDEOGRAM TRACK
# ============================================================================

itrack <- reactive({
  cat("Creating IdeogramTrack for", input$UCSC_Genome, "chromosome", Chromosome(), "\n")
  Gviz::IdeogramTrack(genome=input$UCSC_Genome, chromosome=Chromosome())
})

# ============================================================================
# MAIN PLOT GENERATION (triggered by button)
# ============================================================================

p1 <- eventReactive(input$generate_gene_plot, {
  cat("\n======================================\n")
  cat("=== GENE PLOT GENERATION START ===\n")
  cat("======================================\n")
  
  # Show waiter FIRST
  show_uncoverapp_waiter(
    message = "Generating gene coverage plot...",
    detail = "Creating Gviz tracks - This may take a minute"
  )
  
  # Ensure waiter is hidden even if error occurs
  on.exit(waiter::waiter_hide(), add = TRUE)
  
  # =========================================
  # STEP 1: Validate inputs
  # =========================================
  cat("\nSTEP 1: Validating inputs...\n")
  
  if (is.null(coord())) {
    cat("ERROR: coord() is NULL\n")
    showNotification("Please lookup a gene first!", type="error", duration=5)
    return(NULL)
  }
  
  chr <- Chromosome()
  if (is.null(chr) || chr == "") {
    cat("ERROR: Chromosome is empty\n")
    showNotification("Invalid chromosome!", type="error", duration=5)
    return(NULL)
  }
  
  cat("Chromosome:", chr, "\n")
  
  # =========================================
  # STEP 2: Get gene coordinates
  # =========================================
  cat("\nSTEP 2: Getting gene coordinates...\n")
  
  coord_data <- coord()
  coord_chr <- coord_data[coord_data$seqnames == chr, ]
  
  if (nrow(coord_chr) == 0) {
    cat("ERROR: No coordinates for chromosome", chr, "\n")
    showNotification(paste("Gene not found on", chr), type="error", duration=5)
    return(NULL)
  }

  start_gene <- min(coord_chr$start)
  end_gene <- max(coord_chr$end)
  cat("Gene region:", start_gene, "-", end_gene, "\n")
  
  # =========================================
  # STEP 3: Create tracks (with error handling)
  # =========================================
  cat("\nSTEP 3: Creating visualization tracks...\n")
  
  # 3a. Ideogram
  cat("  Creating ideogram...\n")
  it <- tryCatch({
    itrack()
  }, error = function(e) {
    cat("ERROR creating ideogram:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(it)) {
    showNotification("Failed to create chromosome ideogram", type="error", duration=5)
    return(NULL)
  }
  
  # 3b. Genome axis
  cat("  Creating genome axis...\n")
  gtrack <- Gviz::GenomeAxisTrack()
  
  # 3c. Data tracks
  cat("  Creating low coverage track...\n")
  dt_low <- dtrack1()
  
  cat("  Creating high coverage track...\n")
  dt_high <- dtrackHigh()
  
  cat("  Creating overlay track...\n")
  ot <- Gviz::OverlayTrack(trackList = list(dt_low, dt_high))
  
  # 3d. Gene annotation
  cat("  Creating gene region track...\n")
  gr_ex_track <- tryCatch({
    Gviz::GeneRegionTrack(
      txdb(),
      chromosome = chr,
      start = start_gene, 
      end = end_gene,
      showId = FALSE,
      name = "Gene Annotation",
      cex.title = 1.2
    )
  }, error = function(e) {
    cat("ERROR creating gene track:", e$message, "\n")
    return(NULL)
  })
  cat("DEBUG: gr_ex_track@range has", length(gr_ex_track@range), "elements\n")
  cat("DEBUG: Names:", paste(names(gr_ex_track@range), collapse=", "), "\n")
  transcript_list <<- NULL
  if (length(gr_ex_track@range) > 0) {
    # Per GRanges, i metadati sono in mcols()
    transcript_metadata <- GenomicRanges::mcols(gr_ex_track@range)
    
    cat("DEBUG: Metadata columns:", paste(colnames(transcript_metadata), collapse=", "), "\n")
    
    # Prova varie colonne possibili
    if ("transcript" %in% colnames(transcript_metadata)) {
      transcript_ids <- as.character(transcript_metadata$transcript)
    } else if ("tx_name" %in% colnames(transcript_metadata)) {
      transcript_ids <- as.character(transcript_metadata$tx_name)
    } else if ("tx_id" %in% colnames(transcript_metadata)) {
      transcript_ids <- as.character(transcript_metadata$tx_id)
    } else {
      transcript_ids <- NULL
    }
    
    cat("DEBUG: Found transcript IDs:", paste(head(transcript_ids, 5), collapse=", "), "\n")
    cat("DEBUG: Total transcripts:", length(unique(transcript_ids)), "\n")
    
    # ✅ SALVA i transcript_ids in una variabile accessibile dopo plotTracks
    if (!is.null(transcript_ids)) {
      unique_transcripts <- unique(transcript_ids)
      transcript_list <<- unique_transcripts
    } else {
      unique_transcripts <- NULL
    }
  } else {
    unique_transcripts <- NULL
  }

  Gviz::displayPars(gr_ex_track) <- list(
    collapse = FALSE,
    collapseTranscripts = FALSE,
    showId = FALSE,
    cex = 1.0,
    fontsize = 11
  )
  
  if (is.null(gr_ex_track)) {
    showNotification("Failed to create gene annotation track", type="error", duration=5)
    return(NULL)
  }
  
  # =========================================
  # STEP 4: Calculate plot parameters
  # =========================================
  cat("\nSTEP 4: Calculating plot parameters...\n")
  
  # Get y-axis limits from data
  ylims <- tryCatch({
    all_data <- c(dt_low@data, dt_high@data)
    if (length(all_data) > 0) {
      grDevices::extendrange(range(all_data, na.rm = TRUE))
    } else {
      c(0, 100)  # default range
    }
  }, error = function(e) {
    cat("WARNING: Could not calculate ylims, using default\n")
    c(0, 100)
  })
  
  cat("Y-axis limits:", ylims, "\n")
  
  # Get threshold
  threshold <- as.numeric(input$coverage_co)
  cat("Coverage threshold:", threshold, "\n")
  
  Gviz::displayPars(it) <- list(
    background.title = "white",
    col.title = "black"
  )
  
  # =========================================
  # STEP 5: Generate plot
  # =========================================
  cat("\nSTEP 5: Rendering plot...\n")
  cat("This may take 30-60 seconds for large regions...\n")
  par(mar = c(5, 15, 4, 2) + 0.1)
  plot_output <- tryCatch({
    Gviz::plotTracks(
      list(it, gtrack, ot, gr_ex_track),
      from = start_gene,
      to = end_gene,
      innerMargin = 5,   
      chromosome = chr,
      reverseStrand = TRUE,
      ylim = ylims,
      type = "histogram",
      baseline = threshold,
      col.baseline = "black",
      lwd.baseline = 1,
      sizes = c(0.5, 1, 3.5, 3.5),
      # background.title = "gray50",  # <-- Grigio scuro
      # col.title = "white",          # <-- Scritte bianche
      col.axis = "gray50",
      cex.axis = 0.8,
      fontsize = 14,
      cex.title = 0.9
    )
  }, error = function(e) {
    cat("\n❌ ERROR rendering plot:", e$message, "\n")
    showNotification(
      paste("Plot generation failed:", e$message), 
      type="error", 
      duration=10
    )
    return(NULL)
  })
  
  cat("\nPlot rendered successfully!\n")
  cat("\n======================================\n")
  cat("=== GENE PLOT GENERATION COMPLETE ===\n")
  cat("======================================\n\n")
  
  # Return the plot output (plotTracks returns NULL but draws as side effect)
  invisible(plot_output)
  
}, ignoreNULL = TRUE, ignoreInit = TRUE)