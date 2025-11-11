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
    
    if (ncol(tmp) >= 3) {
      colnames(tmp)[1:3] <- c("chromosome", "start", "end")
    }
    
    tmp$chromosome <- paste0("chr", sub("^chr", "", as.character(tmp$chromosome)))
    tmp$start <- as.integer(tmp$start)
    tmp$end <- as.integer(tmp$end)
    
    if (any(grepl("^sample_", names(tmp)))) {
      cat("Found existing sample_* columns\n")
      return(tmp)
    }
    
    if (ncol(tmp) >= 4) {
      sample_label <- tools::file_path_sans_ext(basename(file_info$name))
      cov <- suppressWarnings(as.integer(tmp[[4]]))
      
      new_df <- data.frame(
        chromosome = tmp$chromosome,
        start = tmp$start,
        end = tmp$end,
        cov = cov,
        stringsAsFactors = FALSE
      )
      colnames(new_df)[4] <- paste0("sample_", sample_label)
      
      cat("Created sample column:", colnames(new_df)[4], "\n")
      return(new_df)
    }
    
    return(NULL)
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
  req(input$coverage_co, input$Sample)
  
  df <- get_sample_data(mydata(), input$Sample)
  req(df)
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
    return(result)
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 2: GENE NAME
  # ══════════════════════════════════════════════════════════════════════
  if (input$filter_by == "gene") {
    req(input$Gene_name)
    
    chr_val <- get_chromosome_from_gene(input$Gene_name, input$UCSC_Genome)
    
    if (is.null(chr_val)) {
      showNotification(
        paste("Could not find chromosome for gene:", input$Gene_name),
        type = "error",
        duration = 5
      )
      return(data.frame())
    }
    
    cat("Mode: GENE | Gene:", input$Gene_name, "| Chr:", chr_val, "\n")
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, chromosome == chr_val)
    } else {
      dplyr::filter(df, chromosome == chr_val, coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
    return(result)
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 3: SINGLE CHROMOSOME
  # ══════════════════════════════════════════════════════════════════════
  if (input$filter_by == "chromosome") {
    req(input$Chromosome)
    chr_val <- input$Chromosome
    
    cat("Mode: CHROMOSOME | Chr:", chr_val, "\n")
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, chromosome == chr_val)
    } else {
      dplyr::filter(df, chromosome == chr_val, coverage <= as.numeric(thr))
    }
    
    cat("Result:", nrow(result), "rows\n")
    return(result)
  }
  
  # ══════════════════════════════════════════════════════════════════════
  # BRANCH 4: REGION COORDINATES
  # ══════════════════════════════════════════════════════════════════════
  # REGION overlap or exact match filter
  if (input$filter_by == "region") {
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
  
    req(region_parts)
  
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
    return(result)
  }

  
  # Fallback
  showNotification("Invalid filter_by selection", type = "error")
  return(data.frame())
  
}, ignoreNULL = TRUE, ignoreInit = TRUE)


# ============================================================================
# FILTERED HIGH (for completeness, also button-triggered)
# ============================================================================


filtered_high <- eventReactive(input$calc_low_coverage, {
  req(input$coverage_co, input$Sample)
  
  df <- get_sample_data(mydata(), input$Sample)
  req(df)
  thr <- input$coverage_co
  
  # Same logic as filtered_low but with coverage > threshold
  
  if (input$filter_by == "all_chr") {
    result <- if (identical(thr, "all")) {
      df
    } else {
      dplyr::filter(df, coverage > as.numeric(thr))
    }
    return(result)
  }
  
  if (input$filter_by == "gene") {
    req(input$Gene_name)
    chr_val <- get_chromosome_from_gene(input$Gene_name, input$UCSC_Genome)
    req(chr_val)
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, chromosome == chr_val)
    } else {
      dplyr::filter(df, chromosome == chr_val, coverage > as.numeric(thr))
    }
    return(result)
  }
  
  if (input$filter_by == "chromosome") {
    req(input$Chromosome)
    chr_val <- input$Chromosome
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, chromosome == chr_val)
    } else {
      dplyr::filter(df, chromosome == chr_val, coverage > as.numeric(thr))
    }
    return(result)
  }
  
  if (input$filter_by == "region") {
    req(input$query_Database)
    region_parts <- tryCatch({
      parts <- strsplit(input$query_Database, "[:-]")[[1]]
      chr_part <- parts[1]
      if (!grepl("^chr", chr_part)) chr_part <- paste0("chr", chr_part)
      list(chr = chr_part, start = as.integer(parts[2]), end = as.integer(parts[3]))
    }, error = function(e) NULL)
    
    req(region_parts)
    
    result <- if (identical(thr, "all")) {
      dplyr::filter(df, 
                    chromosome == region_parts$chr,
                    start >= region_parts$start,
                    end <= region_parts$end)
    } else {
      dplyr::filter(df, 
                    chromosome == region_parts$chr,
                    start >= region_parts$start,
                    end <= region_parts$end,
                    coverage > as.numeric(thr))
    }
    return(result)
  }
  
  return(data.frame())
  
}, ignoreNULL = TRUE, ignoreInit = TRUE)


# ============================================================================
# ALIAS per compatibilità con compute-annotation.R
# ============================================================================


filtered_low_nucl <- filtered_low
filtered_high_nucl <- filtered_high