# ============================================================================
# compute-reactiveDF.R - GESTIONE DATI BASE
# ============================================================================


# Database reference genome
txdb <- reactive({
  if (input$UCSC_Genome == "hg19") {
    TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    TxDb.Hsapiens.UCSC.hg38.knownGene
  }
})


# ============================================================================
# DATA SOURCE MANAGEMENT
# ============================================================================


data_source <- reactiveVal("none")
raw_upload <- reactiveVal(NULL)


observeEvent(input$file1, {
  req(input$file1)
  raw_upload(input$file1)
  data_source("manual")
})


observeEvent(input$process_coverage, {  
  req(coverage_input())
  data_source("pileup")
})


# ============================================================================
# MAIN DATA REACTIVE
# ============================================================================


mydata <- reactive({
  source <- data_source()
  
  if (source == "none") return(NULL)
  
  if (source == "manual") {
    req(raw_upload())
    file_info <- raw_upload()
    
    cat("\n=== PROCESSING MANUAL UPLOAD ===\n")
    cat("File:", file_info$name, "\n")
    
    tmp <- read.table(file_info$datapath,
                      header = input$header, 
                      stringsAsFactors = FALSE)
    
    cat("Raw dimensions:", dim(tmp), "\n")
    cat("Original columns:", paste(colnames(tmp), collapse=", "), "\n")
    
    # Standardize first 3 columns
    if (ncol(tmp) >= 3) {
      colnames(tmp)[1:3] <- c("chromosome", "start", "end")
    }
    
    tmp$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp$chromosome)))
    tmp$start <- as.integer(tmp$start)
    tmp$end <- as.integer(tmp$end)
    
    # ✅ FIX: Rename ALL columns starting with X to sample_ (in one go!)
    if (ncol(tmp) >= 4) {
      for (i in 4:ncol(tmp)) {
        col_name <- colnames(tmp)[i]
        
        # If starts with X followed by digits → replace X with sample_
        if (grepl("^X\\d", col_name)) {
          colnames(tmp)[i] <- sub("^X", "sample_", col_name)
        } 
        # If doesn't have sample_ prefix → add it
        else if (!grepl("^sample_", col_name)) {
          colnames(tmp)[i] <- paste0("sample_", col_name)
        }
        
        # Convert to integer
        # Convert to numeric and verify
        tmp[[i]] <- suppressWarnings(as.numeric(tmp[[i]]))
        
        # Debug: check if conversion worked
        if (all(tmp[[i]] == 0, na.rm = TRUE)) {
          cat("WARNING: Column", colnames(tmp)[i], "converted to all zeros! Check input format.\n")
        }
      }
    }
    
    cat("Final columns:", paste(colnames(tmp), collapse=", "), "\n")
    return(tmp)
  }
  
  if (source == "pileup") {
    req(coverage_input())
    
    cat("\n=== PROCESSING PILEUP DATA ===\n")
    
    tmp_pileup <- coverage_input()
    cat("Coverage dimensions:", dim(tmp_pileup), "\n")
    if ("SYMBOL" %in% colnames(tmp_pileup)) {
        cat("Removing SYMBOL column (only needed for stat_summ)\n")
        tmp_pileup <- tmp_pileup %>% dplyr::select(-SYMBOL)
    }
    colnames(tmp_pileup)[1:3] <- c("chromosome", "start", "end")
    tmp_pileup$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp_pileup$chromosome)))
    tmp_pileup$chromosome <- as.character(tmp_pileup$chromosome)
    tmp_pileup$start <- as.integer(tmp_pileup$start)
    tmp_pileup$end <- as.integer(tmp_pileup$end)
    
    data_cols <- setdiff(colnames(tmp_pileup), c("chromosome", "start", "end"))
    ncols <- length(data_cols)
    samples <- name_sample()
    k <- length(samples)
    
    cat("Data columns:", ncols, "| Samples:", k, "\n")
    
    if (ncols == k && k > 0) {
      colnames(tmp_pileup)[-1:-3] <- paste0("sample_", samples)
      cat("Applied BED naming\n")
    } else if (ncols == 2 * k && k > 0) {
      interleaved <- as.vector(rbind(
        paste0("sample_", samples), 
        paste0("nucleotide_", samples)
      ))
      colnames(tmp_pileup)[-1:-3] <- interleaved
      cat("Applied BAM naming\n")
    } else {
      colnames(tmp_pileup)[-1:-3] <- paste0("sample_", seq_len(ncols))
      cat("WARNING: Using fallback naming\n")
    }
    
    cat("Final columns:", paste(colnames(tmp_pileup), collapse=", "), "\n")
    
    cov_cols <- grep("^sample_", colnames(tmp_pileup), value = TRUE)
    tmp_pileup[cov_cols] <- lapply(tmp_pileup[cov_cols], as.integer)
    
    return(tmp_pileup)
  }
  
  return(NULL)
})


# ============================================================================
# SAMPLE EXTRACTION (pure function)
# ============================================================================

get_sample_data <- function(dat, sample_name) {
  if (is.null(dat) || is.null(sample_name) || sample_name == "") return(NULL)
  
  cat("\n=== EXTRACTING SAMPLE ===\n")
  cat("Requested:", sample_name, "\n")
  
  sel_base <- tools::file_path_sans_ext(sample_name)
  sel_clean <- sub("^sample_", "", sample_name)
  sel_clean_base <- tools::file_path_sans_ext(sel_clean)
  
  candidates <- unique(c(
    sample_name, sel_base, sel_clean, sel_clean_base,
    paste0("sample_", sample_name),
    paste0("sample_", sel_base),
    paste0("sample_", sel_clean),
    paste0("sample_", sel_clean_base)
  ))
  
  col <- intersect(candidates, names(dat))
  
  if (length(col) == 0) {
    cat("ERROR: No match found\n")
    cat("Available:", paste(names(dat), collapse=", "), "\n")
    return(NULL)
  }
  
  col <- col[1]
  cat("Selected:", col, "\n")
  
  out <- dat[, c("chromosome", "start", "end", col), drop = FALSE]
  names(out)[4] <- "coverage"
  
  cat("Extracted", nrow(out), "rows\n")
  
  return(out)
}

# ============================================================================
# HELPER: Get chromosome from gene name
# ============================================================================


get_chromosome_from_gene <- function(gene_name, genome = "hg19") {
  tryCatch({
    gene_info <- suppressMessages(
      AnnotationDbi::select(Homo.sapiens, 
                           keys = gene_name,
                           columns = c("SYMBOL", "TXCHROM"),
                           keytype = "SYMBOL")
    )
    
    if (!is.null(gene_info) && nrow(gene_info) > 0) {
      chrom <- unique(gene_info$TXCHROM)[1]
      if (!grepl("^chr", chrom)) chrom <- paste0("chr", chrom)
      return(chrom)
    }
    
    return(NULL)
    
  }, error = function(e) {
    cat("Error in get_chromosome_from_gene:", conditionMessage(e), "\n")
    return(NULL)
  })
}


# ============================================================================
# FILTERING REACTIVES (CONTROLLED BY ACTION BUTTON)
# ============================================================================


filtered_low <- eventReactive(input$calc_low_coverage, {
  cat("\n=== BUTTON PRESSED: calc_low_coverage ===\n")
  
  # ✅ 1. MOSTRA WAITER
  show_uncoverapp_waiter(
    message = "Calculating low coverage regions...",
    detail = "This may take a few minutes"
  )
  
  # ✅ 2. DISABILITA BOTTONE
  shinyjs::disable("calc_low_coverage")
  
  # ✅ 3. CAMBIA TAB
  updateTabsetPanel(session, "tabSet", selected = "Low-coverage positions")
  
  # ✅ 4. VALIDAZIONE CON ERROR HANDLING
  if (is.null(input$coverage_co) || input$coverage_co == "") {
    waiter::waiter_hide()
    shinyjs::enable("calc_low_coverage")
    showNotification("Please select a coverage threshold", type = "error", duration = 5)
    return(data.frame())
  }
  
  if (is.null(input$Sample) || input$Sample == "") {
    waiter::waiter_hide()
    shinyjs::enable("calc_low_coverage")
    showNotification("Please enter a sample name", type = "error", duration = 5)
    return(data.frame())
  }
  
  df <- get_sample_data(mydata(), input$Sample)
  
  if (is.null(df) || nrow(df) == 0) {
    waiter::waiter_hide()
    shinyjs::enable("calc_low_coverage")
    showNotification("No data available. Please load a coverage file first.", type = "error", duration = 5)
    return(data.frame())
  }
  
  thr <- input$coverage_co
  
  cat("\n=== FILTERING LOW COVERAGE ===\n")
  cat("Filter by:", input$filter_by, "| Threshold:", thr, "\n")
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 1: ALL CHROMOSOMES
  # ══════════════════════════════════════════════════════════════════════
  if (input$filter_by == "all_chr") {
    cat("Mode: ALL CHROMOSOMES\n")
    
    result <- if (identical(thr, "all")) {
      df
    } else {
      dplyr::filter(df, coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
  }
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 2: GENE NAME
  # ══════════════════════════════════════════════════════════════════════
  else if (input$filter_by == "gene") {
    if (is.null(input$Gene_name) || input$Gene_name == "") {
      waiter::waiter_hide()
      shinyjs::enable("calc_low_coverage")
      showNotification("Please enter a gene name", type = "error", duration = 5)
      return(data.frame())
    }
  
    cat("Mode: GENE | Gene:", input$Gene_name, "\n")
  
    # Seleziona il TxDb corretto in base al genome
    txdb_to_use <- if (input$UCSC_Genome == "hg19") {
      TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
    } else {
      TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
    }
  
    cat("Using genome:", input$UCSC_Genome, "\n")
  
    # Ottieni coordinate ESATTE del gene
    gene_coords <- tryCatch({
      # Step 1: Get ENTREZID from gene symbol
      entrez_info <- AnnotationDbi::select(
        org.Hs.eg.db,
        keys = input$Gene_name,
        columns = "ENTREZID",
        keytype = "SYMBOL"
      )
    
      if (is.null(entrez_info) || nrow(entrez_info) == 0) {
        cat("Gene not found in org.Hs.eg.db\n")
        waiter::waiter_hide()
        shinyjs::enable("calc_low_coverage")
        showNotification(paste("Gene not found:", input$Gene_name), type = "error", duration = 5)
        return(NULL)
      }
    
      entrez_id <- entrez_info$ENTREZID[1]
      cat("ENTREZID:", entrez_id, "\n")
    
      # Step 2: Get coordinates from correct TxDb
      gene_info <- AnnotationDbi::select(
        txdb_to_use,
        keys = entrez_id,
        columns = c("TXCHROM", "TXSTART", "TXEND"),
        keytype = "GENEID"
      )
    
      if (!is.null(gene_info) && nrow(gene_info) > 0) {
        chr <- unique(gene_info$TXCHROM)[1]
        if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
      
        list(
          chr = chr,
          start = min(gene_info$TXSTART, na.rm = TRUE),
          end = max(gene_info$TXEND, na.rm = TRUE)
        )
      } else {
        NULL
      }
    }, error = function(e) {
      cat("ERROR getting gene coordinates:", e$message, "\n")
      NULL
    })
  
    # Applica filtro coordinate gene
    if (is.null(gene_coords)) {
      waiter::waiter_hide()
      shinyjs::enable("calc_low_coverage")
      showNotification(
        paste("Could not find coordinates for gene:", input$Gene_name),
        type = "error",
        duration = 5
      )
      return(data.frame())
    }
  
    cat("Gene region:", gene_coords$chr, ":", 
        gene_coords$start, "-", gene_coords$end, "\n")
  
    # Filtra per gene + threshold
    result <- if (identical(thr, "all")) {
      dplyr::filter(
        df,
        chromosome == gene_coords$chr,
        end >= gene_coords$start,
        start <= gene_coords$end
      )
    } else {
      dplyr::filter(
        df,
        chromosome == gene_coords$chr,
        end >= gene_coords$start,
        start <= gene_coords$end,
        coverage <= as.numeric(thr)
      )
    }

    cat("Result:", nrow(result), "rows\n")
  }

  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 3: SINGLE CHROMOSOME
  # ══════════════════════════════════════════════════════════════════════
  else if (input$filter_by == "chromosome") {
    req(input$Chromosome)
    chr_val <- input$Chromosome
    
    cat("Mode: CHROMOSOME | Chr:", chr_val, "\n")
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, chromosome == chr_val)
    } else {
      dplyr::filter(df, chromosome == chr_val, coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 4: REGION COORDINATES
  # ══════════════════════════════════════════════════════════════════════
  # REGION overlap or exact match filter
  else if (input$filter_by == "region") {
    req(input$query_Database)
  
    region_parts <- tryCatch({
      parts <- strsplit(input$query_Database, "[:-]")[[1]]
      if (length(parts) != 3) return(NULL)
      chr_part <- parts[1]
      if (!grepl("^chr", chr_part)) chr_part <- paste0("chr", chr_part)
      start_pos <- as.integer(parts[2])
      end_pos <- as.integer(parts[3])
      if (is.na(start_pos) || is.na(end_pos) || start_pos > end_pos) return(NULL)
      list(chr = chr_part, start = start_pos, end = end_pos)
    }, error = function(e) NULL)
  
    if (is.null(region_parts)) {
      waiter::waiter_hide()
      shinyjs::enable("calc_low_coverage")
      showNotification("Invalid region format. Use chr:start-end", type = "error", duration = 5)
      return(data.frame())
    }
  
    if (identical(thr, "all")) {
      # PROVA match esatto
      result_exact <- dplyr::filter(
        df,
        chromosome == region_parts$chr,
        start == region_parts$start,
        end == region_parts$end
      )
      result <- if (nrow(result_exact) > 0) {
        result_exact
      } else {
        dplyr::filter(
          df, 
          chromosome == region_parts$chr,
          start <= region_parts$end,
          end >= region_parts$start
        )
      }
    } else {
      result_exact <- dplyr::filter(
        df,
        chromosome == region_parts$chr,
        start == region_parts$start,
        end == region_parts$end,
        coverage <= as.numeric(thr)
      )
      result <- if (nrow(result_exact) > 0) {
        result_exact
      } else {
        dplyr::filter(
          df, 
          chromosome == region_parts$chr,
          start <= region_parts$end,
          end >= region_parts$start,
          coverage <= as.numeric(thr)
        )
      }
    }
  
    cat("Result:", nrow(result), "rows\n")
  }

  
  # ══════════════════════════════════════════════════════════════════════
  # FALLBACK
  # ══════════════════════════════════════════════════════════════════════
  else {
    waiter::waiter_hide()
    shinyjs::enable("calc_low_coverage")
    showNotification("Invalid filter_by selection", type = "error", duration = 5)
    return(data.frame())
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # FINAL: Hide waiter and notify (for all successful branches)
  # ══════════════════════════════════════════════════════════════════════
  waiter::waiter_hide()
  shinyjs::enable("calc_low_coverage")
  
  if (!is.null(result) && nrow(result) > 0) {
    showNotification(
      paste("Found", nrow(result), "low coverage positions"),
      type = "message",
      duration = 5
    )
  } else {
    showNotification(
      "No low coverage positions found with current filters",
      type = "warning",
      duration = 5
    )
  }
  
  return(result)
  
}, ignoreNULL = TRUE, ignoreInit = TRUE)


# ============================================================================
# FILTERED HIGH (for completeness, also button-triggered)
# ============================================================================


filtered_high <- filtered_high <- reactive({
  cat("\n=== FILTERED_HIGH START ===\n")
  
  cat("Checking requirements...\n")
  req(input$coverage_co, input$Sample)
  cat("Requirements OK\n")
  
  cat("Filter by:", input$filter_by, "\n")
  
  # Solo per gene mode
  if (input$filter_by != "gene") {
    cat("Not gene mode, returning empty\n")
    return(data.frame())
  }
  
  # CRITICAL FIX: Use SAME LOGIC as filtered_low() to get gene region data
  cat("Getting sample data from mydata()...\n")
  df <- get_sample_data(mydata(), input$Sample)
  req(df)
  cat("Sample data OK, rows:", nrow(df), "\n")
  
  # Get gene coordinates (SAME as in filtered_low)
  txdb_to_use <- if (input$UCSC_Genome == "hg19") {
    TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene
  } else {
    TxDb.Hsapiens.UCSC.hg38.knownGene::TxDb.Hsapiens.UCSC.hg38.knownGene
  }
  
  gene_coords <- tryCatch({
    entrez_info <- AnnotationDbi::select(
      org.Hs.eg.db,
      keys = input$Gene_name,
      columns = "ENTREZID",
      keytype = "SYMBOL"
    )
    
    if (is.null(entrez_info) || nrow(entrez_info) == 0) {
      return(NULL)
    }
    
    entrez_id <- entrez_info$ENTREZID[1]
    gene_info <- AnnotationDbi::select(
      txdb_to_use,
      keys = entrez_id,
      columns = c("TXCHROM", "TXSTART", "TXEND"),
      keytype = "GENEID"
    )
    
    if (!is.null(gene_info) && nrow(gene_info) > 0) {
      chr <- unique(gene_info$TXCHROM)[1]
      if (!grepl("^chr", chr)) chr <- paste0("chr", chr)
      
      list(
        chr = chr,
        start = min(gene_info$TXSTART, na.rm = TRUE),
        end = max(gene_info$TXEND, na.rm = TRUE)
      )
    } else {
      NULL
    }
  }, error = function(e) {
    cat("ERROR getting gene coordinates:", e$message, "\n")
    return(NULL)
  })
  
  if (is.null(gene_coords)) {
    cat("ERROR: Cannot get gene coordinates\n")
    return(data.frame())
  }
  
  chr <- gene_coords$chr
  start_region <- gene_coords$start
  end_region <- gene_coords$end
  
  cat("Gene region:", chr, ":", start_region, "-", end_region, "\n")
  
  thr <- input$coverage_co
  cat("Threshold for HIGH coverage: >", thr, "\n")
  
  # Filter for ALL coverage ABOVE threshold in ENTIRE gene region
  # This gives us the "high coverage" regions to display in blue on the plot
  result <- dplyr::filter(df, 
                          chromosome == chr,
                          end >= start_region, 
                          start <= end_region,
                          coverage > as.numeric(thr))
  
  cat("HIGH coverage positions:", nrow(result), "\n")
  cat("=== FILTERED_HIGH COMPLETE ===\n")
  
  return(result)
  
})


# ============================================================================
# ALIAS per compatibilità con compute-annotation.R
# ============================================================================


filtered_low_nucl <- filtered_low
filtered_high_nucl <- filtered_high