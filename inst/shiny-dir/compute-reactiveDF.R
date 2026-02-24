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
observeEvent(input$pileup, {
  cat("\n=== PILEUP BUTTON PRESSED ===\n")
  
  # Check if coverage data exists
  cov_data <- tryCatch({
    coverage_input()
  }, error = function(e) {
    NULL
  })
  
  if (!is.null(cov_data) && nrow(cov_data) > 0) {
    data_source("pileup")
    showNotification(
      "Coverage data loaded successfully",
      type = "message",
      duration = 3
    )
  } else {
    showNotification(
      "No coverage data available. Please go to 'Coverage Analysis' and press 'Process Coverage' first.",
      type = "warning",
      duration = 5
    )
  }
})
observeEvent(input$process_coverage, {  
  # Reset state
  data_source("none")
  
  # Try to get coverage data
  cov_data <- tryCatch({
    coverage_input()
  }, error = function(e) {
    showNotification(
      paste("⚠️ Errore elaborazione coverage:", e$message),
      type = "error",
      duration = 10
    )
    return(NULL)
  })
  
  # Set data source only if successful
  if (!is.null(cov_data) && nrow(cov_data) > 0) {
    data_source("pileup")
  } else {
    showNotification(
      "⚠️ Coverage elaboration failed. Check files and retry.",
      type = "error",
      duration = 5
    )
  }
})

# ============================================================================
# MAIN DATA REACTIVE
# ============================================================================

mydata <- reactive({
  source <- data_source()
  
  if (source == "none") return(NULL)
  
  # ========== MANUAL UPLOAD ==========
  if (source == "manual") {
    req(raw_upload())
    file_info <- raw_upload()
    
    cat("\n=== PROCESSING MANUAL UPLOAD ===\n")
    cat("File:", file_info$name, "\n")
    
    # Reading file with error handling
    tmp <- tryCatch({
      read.table(file_info$datapath,
                 header = input$header, 
                 stringsAsFactors = FALSE)
    }, error = function(e) {
      validate(need(FALSE, 
        paste("Impossibile leggere il file. Verifica che sia un file di testo valido.\nErrore:", e$message)))
    })
    
    cat("Raw dimensions:", dim(tmp), "\n")
    cat("Original columns:", paste(colnames(tmp), collapse=", "), "\n")
    
    # Column validation
    validate(
      need(ncol(tmp) >= 3,
           paste("Il file deve avere almeno 3 colonne (chromosome, start, end).",
                 "\nTrovate:", ncol(tmp), "colonne.",
                 "\nSuggerimento: verifica il separatore (tab, virgola, spazio) e se hai selezionato 'Header' correttamente."))
    )
    
    # Rename first 3 columns
    colnames(tmp)[1:3] <- c("chromosome", "start", "end")
    
    # Validate columns exist
    validate(
      need("chromosome" %in% colnames(tmp), "Errore interno: colonna chromosome mancante"),
      need("start" %in% colnames(tmp), "Errore interno: colonna start mancante"),
      need("end" %in% colnames(tmp), "Errore interno: colonna end mancante")
    )
    
    # Standardize chromosome names
    tmp$chromosome <- tryCatch({
      paste0("chr", sub("^chr", "", as.character(tmp$chromosome)))
    }, error = function(e) {
      validate(need(FALSE, "Errore nella colonna chromosome. Deve contenere nomi di cromosomi (es. chr1, chr22, X)"))
    })
    
    # Convert start to integer with validation
    tmp$start <- tryCatch({
      as.integer(tmp$start)
    }, error = function(e) {
      validate(need(FALSE, "La colonna 'start' deve contenere solo numeri interi"))
    })
    
    # Check if conversion produced too many NAs
    validate(
      need(sum(is.na(tmp$start)) < nrow(tmp) * 0.1,
           paste("La colonna 'start' contiene troppi valori non numerici.",
                 "\nRighe con errori:", sum(is.na(tmp$start)), "su", nrow(tmp)))
    )
    
    # Convert end to integer with validation
    tmp$end <- tryCatch({
      as.integer(tmp$end)
    }, error = function(e) {
      validate(need(FALSE, "La colonna 'end' deve contenere solo numeri interi"))
    })
    
    # Check if conversion produced too many NAs
    validate(
      need(sum(is.na(tmp$end)) < nrow(tmp) * 0.1,
           paste("La colonna 'end' contiene troppi valori non numerici.",
                 "\nRighe con errori:", sum(is.na(tmp$end)), "su", nrow(tmp)))
    )
    
    # Validate that start < end
    invalid_intervals <- sum(tmp$start >= tmp$end, na.rm = TRUE)
    validate(
      need(invalid_intervals == 0,
           paste("Trovati", invalid_intervals, "intervalli con start >= end.",
                 "\nGli intervalli devono avere start < end."))
    )
    
    # Rename sample columns
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
  
  # ========== PILEUP DATA ==========
  if (source == "pileup") {
    cov_data <- coverage_input()
    
    validate(
      need(!is.null(cov_data), 
          "⚠️ No coverage data available. Be sure to have done:
          • Uploaded file BAM/BED in 'Coverage Analysis'
          • Uploaded gene names or file target
          • Pressed 'Process Coverage' and waited for completion"),
      need(is.data.frame(cov_data) || is.matrix(cov_data),
          "⚠️ Data format not recognized"),
      need(nrow(cov_data) > 0,
          "⚠️ No coverage regions found"),
      need(ncol(cov_data) >= 3,
          "⚠️ Incomplete data: need at least 3 columns")
    )
    
    cat("\n=== PROCESSING PILEUP DATA ===\n")
    tmp_pileup <- cov_data
    cat("Coverage dimensions:", dim(tmp_pileup), "\n")
    
    # Remove SYMBOL if present
    if ("SYMBOL" %in% colnames(tmp_pileup)) {
      cat("Removing SYMBOL column (only needed for stat_summ)\n")
      tmp_pileup <- tmp_pileup %>% dplyr::select(-SYMBOL)
    }
    
    # Standardize columns
    colnames(tmp_pileup)[1:3] <- c("chromosome", "start", "end")
    tmp_pileup$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp_pileup$chromosome)))
    tmp_pileup$chromosome <- as.character(tmp_pileup$chromosome)
    tmp_pileup$start <- as.integer(tmp_pileup$start)
    tmp_pileup$end <- as.integer(tmp_pileup$end)
    
    # Rename sample columns
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
    tmp_pileup[cov_cols] <- lapply(tmp_pileup[cov_cols], function(x) {
      suppressWarnings(as.numeric(x))
    })
    
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
  
  # Show waiter
  show_uncoverapp_waiter(
    message = "Calculating low coverage regions...",
    detail = "This may take a few minutes"
  )
  
  # Disable button
  shinyjs::disable("calc_low_coverage")
  
  # Change tab
  updateTabsetPanel(session, "tabSet", selected = "Low-coverage positions")
  
  # Validation with error handling
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
      dplyr::filter(df, is.na(coverage) | coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
  }
  # ══════════════════════════════════════════════════════════════════════════
  # BRANCH 2: GENE NAME
  # ══════════════════════════════════════════════════════════════════════════
  else if (input$filter_by == "gene") {
    if (is.null(input$Gene_name) || input$Gene_name == "") {
      waiter::waiter_hide()
      shinyjs::enable("calc_low_coverage")
      showNotification("Please enter a gene name", type = "error", duration = 5)
      return(data.frame())
    }
  
    cat("Mode: GENE | Gene:", input$Gene_name, "\n")
  
    # ✅ USE PRE-CALCULATED COORDINATES FROM gene_coordinates()
    gene_coords <- tryCatch({
      gene_coordinates()
    }, error = function(e) {
      cat("ERROR: Could not get gene_coordinates():", e$message, "\n")
      return(NULL)
    })
  
    if (is.null(gene_coords)) {
      waiter::waiter_hide()
      shinyjs::enable("calc_low_coverage")
      showNotification(
        paste("Please lookup gene coordinates first using 'Lookup UCSC Gene' button"),
        type = "error",
        duration = 5
      )
      return(data.frame())
    }
  
    cat("Using global gene coordinates:", gene_coords$chr, ":", 
        gene_coords$start, "-", gene_coords$end, "\n")
  
    # Filter by gene + threshold
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
      is.na(coverage) | coverage <= as.numeric(thr)
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
      dplyr::filter(df, chromosome == chr_val, is.na(coverage) | coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 4: REGION COORDINATES
  # ══════════════════════════════════════════════════════════════════════
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
        is.na(coverage) | coverage <= as.numeric(thr)
      )
      result <- if (nrow(result_exact) > 0) {
        result_exact
      } else {
        dplyr::filter(
          df, 
          chromosome == region_parts$chr,
          start <= region_parts$end,
          end >= region_parts$start,
          is.na(coverage) | coverage <= as.numeric(thr)
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
  # FINAL: Hide waiter and notify
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
# FILTERED HIGH
# ============================================================================

filtered_high <- reactive({
  cat("\n=== FILTERED_HIGH START ===\n")
  
  cat("Checking requirements...\n")
  req(input$coverage_co, input$Sample)
  cat("Requirements OK\n")
  
  cat("Filter by:", input$filter_by, "\n")
  
  # Only for gene mode
  if (input$filter_by != "gene") {
    cat("Not gene mode, returning empty\n")
    return(data.frame())
  }
  
  cat("Getting sample data from mydata()...\n")
  df <- get_sample_data(mydata(), input$Sample)

  validate(
    need(!is.null(df) && nrow(df) > 0,
        "⚠️ No data available for the selected sample.")
  )
  
  cat("Sample data OK, rows:", nrow(df), "\n")
  
  # ✅ USE PRE-CALCULATED COORDINATES FROM gene_coordinates()
  gene_coords <- tryCatch({
    gene_coordinates()
  }, error = function(e) {
    cat("ERROR: Could not get gene_coordinates():", e$message, "\n")
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
  
  # Filter for coverage ABOVE threshold in gene region
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