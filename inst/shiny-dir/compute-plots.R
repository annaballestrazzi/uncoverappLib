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
  
  tryCatch({
    Gviz::IdeogramTrack(genome=input$UCSC_Genome, chromosome=Chromosome())
  }, error = function(e) {
    cat("WARNING: Could not create IdeogramTrack:", conditionMessage(e), "\n")
    cat("This is not critical - the plot will work without the ideogram\n")
    return(NULL)
  })
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
  
  # 3a. Ideogram (optional)
  cat("  Creating ideogram...\n")
  it <- tryCatch({
    itrack()
  }, error = function(e) {
    cat("WARNING: Ideogram track creation failed:", conditionMessage(e), "\n")
    return(NULL)
  })
  
  if (is.null(it)) {
    cat("  ⚠️ Ideogram track not available - continuing without it\n")
  } else {
    cat("  ✓ Ideogram created successfully\n")
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
  
  # # 3d. Gene annotation
  # cat("  Creating gene region track...\n")
  # gr_ex_track <- tryCatch({
  #   Gviz::GeneRegionTrack(
  #     txdb(),
  #     chromosome = chr,
  #     start = start_gene, 
  #     end = end_gene,
  #     showId = FALSE,
  #     name = "Gene Annotation",
  #     cex.title = 1.2
  #   )
  # }, error = function(e) {
  #   cat("ERROR creating gene track:", e$message, "\n")
  #   return(NULL)
  # })

  # 3d. Gene annotation
  cat("  Creating gene region track...\n")
  
  # Add 50kb margin to ensure we catch all transcripts
  # (transcript boundaries may not match exactly with gene boundaries)
  margin <- 50000
  plot_start <- max(1, start_gene - margin)
  plot_end <- end_gene + margin
  
  cat("  Querying wider region:", plot_start, "-", plot_end, 
      "(gene region +", margin/1000, "kb margin)\n")
  
  
  gr_ex_track <- tryCatch({
    # CRITICAL FIX: Query by gene ID instead of coordinates
    # GeneRegionTrack(txdb, chr, start, end) returns 0 features for SHANK3!
    
    # First, get the gene's ENTREZ ID
    gene_name <- input$Gene_name
    entrez_lookup <- OrganismDbi::select(org.Hs.eg.db,
                                        keys = gene_name,
                                        columns = "ENTREZID",
                                        keytype = "SYMBOL")
    
    if (!is.null(entrez_lookup) && nrow(entrez_lookup) > 0) {
      entrez_id <- entrez_lookup$ENTREZID[1]
      cat("  Using gene ID", entrez_id, "to create track\n")
      
      # Query by gene ID to get ALL features for this gene
      track <- Gviz::GeneRegionTrack(
        txdb(),
        chromosome = chr,
        gene = entrez_id,  # KEY: Use gene parameter instead of start/end
        name = "Gene Annotation",
        showId = FALSE,
        cex.title = 1.2
      )
      
      # Subset to our region of interest (with margin)
      if (!is.null(track) && length(track@range) > 0) {
        cat("  Gene query found", length(track@range), "features\n")
        
        # Keep only features in our display region
        feature_ranges <- track@range
        overlaps <- GenomicRanges::start(feature_ranges) <= plot_end &
                    GenomicRanges::end(feature_ranges) >= plot_start
        
        if (sum(overlaps) > 0) {
          track@range <- feature_ranges[overlaps]
          cat("  Filtered to", length(track@range), "features in display region\n")
          track
        } else {
          cat("  WARNING: No features in display region\n")
          track
        }
      } else {
        # Fallback to coordinate query
        cat("  Gene query returned nothing, trying coordinates...\n")
        Gviz::GeneRegionTrack(
          txdb(),
          chromosome = chr,
          start = plot_start,
          end = plot_end,
          showId = FALSE,
          name = "Gene Annotation",
          cex.title = 1.2
        )
      }
    } else {
      # Can't find gene ID, use coordinate query
      cat("  Could not find ENTREZ ID, using coordinate query\n")
      Gviz::GeneRegionTrack(
        txdb(),
        chromosome = chr,
        start = plot_start,
        end = plot_end,
        showId = FALSE,
        name = "Gene Annotation",
        cex.title = 1.2
      )
    }
  }, error = function(e) {
    cat("ERROR creating gene track:", e$message, "\n")
    return(NULL)
  })
  
  cat("DEBUG: gr_ex_track@range has", length(gr_ex_track@range), "elements\n")
  cat("DEBUG: Names:", paste(names(gr_ex_track@range), collapse=", "), "\n")
  
  # ============================================================================
  # USE TRANSCRIPT IDs ALREADY COMPUTED BY coord() - No duplicate query!
  # ============================================================================
  gene_transcript_ids <<- NULL
  
  coord_data <- coord()
  
  if (!is.null(coord_data) && "transcript_ids" %in% colnames(coord_data)) {
    transcript_ids <- coord_data$transcript_ids[[1]]
    
    if (!is.null(transcript_ids) && length(transcript_ids) > 0) {
      cat("Using transcript IDs from coord():", length(transcript_ids), "transcripts\n")
      cat("First 5:", paste(head(transcript_ids, 5), collapse=", "), "\n")
      gene_transcript_ids <<- transcript_ids
    } else {
      cat("WARNING: No transcript IDs found in coord()\n")
    }
  } else {
    cat("WARNING: coord() does not contain transcript_ids column\n")
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
  
  # Set ideogram display parameters only if it exists
  if (!is.null(it)) {
    Gviz::displayPars(it) <- list(
      background.title = "white",
      col.title = "black"
    )
  }
  
  # =========================================
  # STEP 5: Generate plot
  # =========================================
  cat("\nSTEP 5: Rendering plot...\n")
  cat("This may take 30-60 seconds for large regions...\n")
  par(mar = c(5, 15, 4, 2) + 0.1)
  
  # Build track list - include ideogram only if available
  track_list <- if (!is.null(it)) {
    list(it, gtrack, ot, gr_ex_track)
  } else {
    list(gtrack, ot, gr_ex_track)
  }
  
  # Adjust sizes based on whether ideogram is present
  track_sizes <- if (!is.null(it)) {
    c(0.5, 1, 3.5, 3.5)  # With ideogram
  } else {
    c(1, 3.5, 3.5)  # Without ideogram
  }
  
  plot_output <- tryCatch({
    Gviz::plotTracks(
      track_list,
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
      sizes = track_sizes,
      col.axis = "gray50",
      cex.axis = 0.8,
      fontsize = 14,
      cex.title = 0.9
    )
  }, error = function(e) {
    cat("\nâŒ ERROR rendering plot:", e$message, "\n")
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

# ============================================================================
# DOWNLOAD: GENE COVERAGE PLOT
# ============================================================================
output$download_plot <- downloadHandler(
  filename = function() {
    gene_name <- if (!is.null(input$Gene_name) && input$Gene_name != "") {
      input$Gene_name
    } else {
      "gene"
    }
    paste0('coverage_plot_', gene_name, '_', Sys.Date(), '.png')
  },
  content = function(file) {
    req(plot_recorded())
    png(file, width = 3600, height = 2400, res = 300)
    on.exit(dev.off(), add = TRUE)
    replayPlot(plot_recorded())
    cat("=== PLOT DOWNLOADED ===\n")
    cat("File:", file, "\n")
  },
  contentType = "application/octet-stream"
)

# ============================================================================
# OUTPUT: TRANSCRIPT LIST (displayed below plot)
# ============================================================================
output$transcript_list_text <- renderUI({
  # CRITICAL: Wait for plot to complete, not just button press
  req(p1())  # This ensures p1() has finished executing
  
  # Small delay to ensure gene_transcript_ids is set
  Sys.sleep(0.1)
  
  transcripts <- tryCatch({
    # Try to get gene_transcript_ids from global environment
    if (exists("gene_transcript_ids", envir = .GlobalEnv)) {
      tlist <- get("gene_transcript_ids", envir = .GlobalEnv)
      
      # Check if it's a function (naming conflict)
      if (is.function(tlist)) {
        cat("WARNING: gene_transcript_ids is a function, not a vector!\n")
        return(NULL)
      }
      
      # Check if it's a vector or list
      if (!is.vector(tlist) && !is.list(tlist)) {
        cat("WARNING: gene_transcript_ids is not a vector/list!\n")
        return(NULL)
      }
      
      # Check for empty
      if (length(tlist) == 0) {
        cat("WARNING: gene_transcript_ids is empty\n")
        return(NULL)
      }
      
      return(tlist)
      
    } else {
      cat("WARNING: gene_transcript_ids does not exist\n")
      return(NULL)
    }
  }, error = function(e) {
    cat("ERROR getting gene_transcript_ids:", e$message, "\n")
    return(NULL)
  })
  
  # Check if transcripts were retrieved successfully
  if (is.null(transcripts) || length(transcripts) == 0) {
    return(shiny::tags$p("No transcript information available", style = "color: gray;"))
  }
  
  # Ensure transcripts is a character vector
  transcripts <- as.character(transcripts)
  
  cat("Displaying", length(transcripts), "transcripts\n")
  
  # Build the UI
  shiny::tagList(
    shiny::tags$div(style = "margin-top: 150px;"),
    shiny::h4(paste0("Transcripts (", length(transcripts), ")")),
    shiny::tags$div(
      style = "border: 1px solid #ddd; 
              padding: 10px; 
              background: #f9f9f9; 
              max-height: 200px; 
              overflow-y: auto;
              width: 80%;
              margin: 0 auto;",
      shiny::tags$ol(
        lapply(seq_along(transcripts), function(i) {
          shiny::tags$li(transcripts[i])
        })
      )
    )
  )
})