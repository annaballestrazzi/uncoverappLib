server <- function (input, output, session){


  suppressPackageStartupMessages({
    require(dplyr)
    require(Gviz)
    require(Homo.sapiens)
    require(stringr)
    require(bedr)
    require(EnsDb.Hsapiens.v75)
    require(EnsDb.Hsapiens.v86)
    require(DT)
    require(shinycssloaders)
    require(waiter) 
  })

  options(shiny.maxRequestSize=30*1024^2)
    # ============================================================================
  # LOAD OMIM GENE LIST (caricato UNA SOLA VOLTA per tutta la sessione)
  # ============================================================================
  omim_file <- system.file("extdata", "sys_ndd_2025_subset.tsv", 
                           package = "uncoverappLib")
  
  if (file.exists(omim_file)) {
    OMIM_DATA <- read.table(omim_file, header = TRUE, sep = "\t", 
                            stringsAsFactors = FALSE, quote = "", fill = TRUE)
    OMIM_GENES <- unique(OMIM_DATA$SYMBOL)
    cat("✓ Loaded", length(OMIM_GENES), "OMIM genes from database\n")
  } else {
    cat("WARNING: OMIM file not found, OMIM highlighting disabled\n")
    OMIM_DATA <- data.frame()
    OMIM_GENES <- character(0)
  }

  script1 <- system.file(
    "extdata",
    "Rpreprocessing.R",
    package = "uncoverappLib")

  # ============================================================================
  # DOWNLOAD DEPENDENCY SCRIPT
  # ============================================================================
  
  output$dependence = downloadHandler(
    filename="Rpreprocessing.R",
    content=function(file){
      file.copy(script1,file)
    }
  )
  
  # ============================================================================
  # SOURCE REACTIVE SCRIPTS
  # ============================================================================
  
  source('compute-preprocess.R', local= TRUE)
  source('compute-reactiveDF.R', local= TRUE)      # mydata(), mysample(), filtered_low()
  source('compute-annotation.R', local=TRUE)       # intBED(), annotated_variants_data(), condform_table()
  source('compute-tables.R', local= TRUE)
  source('compute-plots.R', local=TRUE)
  source('compute-maxAF.R', local=TRUE)            # data(), uncover_maxaf_data(), uncover_maxaf()
  source('compute-binomial.R', local=TRUE)
  source('waiter-helpers.R', local=TRUE)

  # ============================================================================
  # DYNAMIC SAMPLE SELECTOR - Update dropdown when file is loaded
  # ============================================================================

  observeEvent(mydata(), {
    req(mydata())
    
    cat("\n=== UPDATING SAMPLE DROPDOWN ===\n")
    
    data <- mydata()
    
    # Extract sample column names
    coord_cols <- c("chromosome", "start", "end")
    all_cols <- colnames(data)
    sample_cols <- setdiff(all_cols, coord_cols)
    
    # Filter for sample/count columns
    sample_cols <- sample_cols[grepl("^(sample_|count_)", sample_cols)]
    
    cat("Available samples:", length(sample_cols), "\n")
    if (length(sample_cols) > 0) {
      cat("Sample names:", paste(sample_cols, collapse = ", "), "\n")
    }
    
    if (length(sample_cols) == 0) {
      # No samples found
      cat("WARNING: No sample columns detected\n")
      updateSelectInput(
        session, 
        "Sample",
        label = "Sample (none found)",
        choices = c("No samples available" = ""),
        selected = NULL
      )
      
      showNotification(
        "Warning: No sample columns found in loaded file",
        type = "warning",
        duration = 5
      )
      
    } else {
      # Update dropdown with available samples
      new_label <- paste0("Sample (", length(sample_cols), " available)")
      
      updateSelectInput(
        session, 
        "Sample",
        label = new_label,
        choices = sample_cols,
        selected = sample_cols[1]
      )
      
      cat("✓ Dropdown updated successfully\n")
      cat("✓ Auto-selected:", sample_cols[1], "\n\n")
      
      showNotification(
        paste("Loaded", length(sample_cols), "sample(s)"),
        type = "message",
        duration = 2
      )
    }
    
  }, ignoreInit = TRUE, ignoreNULL = TRUE)

  # ============================================================================
  # OUTPUT: low_coverage_ready (for annotation button visibility)
  # ============================================================================
  output$low_coverage_ready <- reactive({
    result <- !is.null(filtered_low_nucl()) && nrow(filtered_low_nucl()) > 0
    
    # Debug (opzionale - rimuovi dopo test)
    cat("Low coverage ready:", result, "\n")
    if (!is.null(filtered_low_nucl())) {
      cat("  Rows:", nrow(filtered_low_nucl()), "\n")
    }
    
    return(result)
  })
  outputOptions(output, "low_coverage_ready", suspendWhenHidden = FALSE)
  # ============================================================================
  # OUTPUT: annotation_ready (for download button visibility)
  # ============================================================================
  output$annotation_ready <- reactive({
    !is.null(annotated_variants_data()) && nrow(annotated_variants_data()) > 0
  })
  outputOptions(output, "annotation_ready", suspendWhenHidden = FALSE)

  # ============================================================================
  # AUTO-SWITCH TO TAB when buttons are pressed
  # ============================================================================
  
  # Disable download button at start
  shinyjs::disable("summary")
  
  # Handle waiter for process_coverage button
  # Enable download button when coverage_input completes
  observeEvent(coverage_input(), {
    shinyjs::enable("summary")  # ✓ Abilita download solo dopo processing
  })
  
  observeEvent(input$calc_low_coverage, {
    updateTabsetPanel(session, "tabSet", selected = "Low-coverage positions")
  })

  observeEvent(input$calc_annotations, {
    updateTabsetPanel(session, "tabSet", selected = "Annotations on low-coverage positions")
  })
  
  # ============================================================================
  # OUTPUT: COVERAGE INPUT TABLE
  # ============================================================================
  
  output$input1 <- renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    
    # Wait for the process button to be pressed
    req(input$process_coverage > 0)
    
    start_time <- Sys.time()
    
    tryCatch({
      df <- coverage_input()
      
      end_time <- Sys.time()
      cat(paste("Table created in:", 
                round(difftime(end_time, start_time, units = "secs"), 1), 
                "seconds\n"))
      
      # Notifica successo
      showNotification(
        paste0("Input table ready: ", nrow(df), " intervals processed"),
        type = "message",
        duration = 3
      )
      
      return(df)
      
    }, error = function(e) {
      showNotification(
        paste("Error creating input table:", e$message),
        type = "error",
        duration = 5
      )
      return(NULL)
    })
  })

  # ============================================================================
  # DOWNLOAD: STATISTICAL SUMMARY
  # ============================================================================
  
  output$summary <- downloadHandler(
    filename = function() {
      paste('statistical_summary_', Sys.Date(), '.txt', sep='')
    },
    
    content = function(file){
      cat("\n=== DOWNLOAD HANDLER START ===\n")
      
      tryCatch({
        # Step 1: Calcola statistiche
        cat("Calculating statistics...\n")
        original_data <- stat_summ()
        
        if (is.null(original_data) || nrow(original_data) == 0) {
          showNotification("No data to export!", type = "error", duration = 5)
          return(NULL)
        }
        
        cat("  Rows:", nrow(original_data), "\n")
        
        # Step 2: Usa OMIM già caricato
        cat("Using pre-loaded OMIM annotations...\n")
        omim_gene <- OMIM_DATA
        
        cat("  OMIM rows:", nrow(omim_gene), "\n")
        
        # Check for SYMBOL column
        if (!"SYMBOL" %in% colnames(original_data)) {
          stop("ERROR: original_data missing SYMBOL column!")
        }
        if (!"SYMBOL" %in% colnames(omim_gene)) {
          stop("ERROR: omim_gene missing SYMBOL column!")
        }
        
        # Step 3: Merge dati
        cat("Merging data with annotations...\n")
        joined_data <- merge(
          original_data, 
          omim_gene, 
          by = "SYMBOL", 
          all.x = TRUE,
          suffixes = c("", ".omim")
        )
        
        cat("  Rows after merge:", nrow(joined_data), "\n")
        
        # Reorder columns
        stat_cols <- c("SYMBOL", "sample", "Total_bases", "Mean_coverage", 
                      "Median_coverage", "number_of_intervals_under_20x", 
                      "bases_under_20x", "percentage_bases_under_20x")
        omim_cols <- setdiff(colnames(joined_data), stat_cols)
        joined_data <- joined_data[, c(stat_cols, omim_cols)]
        
        # Step 4: Scrivi file
        cat("Writing file...\n")
        write.table(
          joined_data, 
          file, 
          sep = '\t', 
          quote = FALSE, 
          row.names = FALSE, 
          col.names = TRUE,
          na = "NA"
        )
        
        # Notifica successo
        showNotification(
          "Statistical summary created successfully!",
          type = "message",
          duration = 3
        )
        
        cat("=== DOWNLOAD COMPLETE ===\n")
        cat("File saved:", file, "\n\n")
        
      }, error = function(e) {
        showNotification(
          paste("Download error:", e$message),
          type = "error",
          duration = 5
        )
        cat("ERROR:", e$message, "\n")
      })
    }
  )


  # ============================================================================
  # OUTPUT: MYDATA TABLE (full dataset)
  # ============================================================================
  
  output$text <- DT::renderDataTable({
    validate(need(ncol(mydata()) != "0", "Please, upload your file"))
    mydata()
  })

  # ============================================================================
  # OUTPUT: FILTERED COVERAGE TABLE
  # ============================================================================
  
  output$text_cv <- DT::renderDataTable({
    # Non validare input qui, filtered_low() gestisce tutto con req()
    data <- tryCatch(filtered_low(), error = function(e) {
      cat("ERROR in text_cv:", conditionMessage(e), "\n")
      return(NULL)
    })
  
    # Mostra messaggio se non ci sono dati
    validate(
      need(!is.null(data) && nrow(data) > 0,
           "No low coverage data available. Press 'Calculate Low Coverage Regions' button.")
    )
  
    return(data)
  })
  # ============================================================================
  # OUTPUT: ALL GENE COVERAGE PLOT
  # ============================================================================
    output$all_gene <- renderPlot({
      cat("\n=== RENDER PLOT TRIGGERED ===\n")
      p1()
      cat("PLOT RENDERED!\n")
    }, height = 600, width = 800)
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
        # Apri dispositivo grafico PNG
        png(file, width = 1200, height = 800, res = 120)
        
        # Rigenera il plot
        p1()
        
        # Chiudi dispositivo
        dev.off()
      }
    )
    output$transcript_list_text <- renderUI({
      req(input$generate_gene_plot > 0)
      Sys.sleep(0.3)
      
      transcripts <- tryCatch({
        get("transcript_list", envir = .GlobalEnv)
      }, error = function(e) {
        NULL
      })
      
      # Check se esiste e non è vuoto
      if (is.null(transcripts) || length(transcripts) == 0) {
        return(tags$p("No transcript information available", style = "color: gray;"))
      }
      
      tagList(
        tags$div(style = "margin-top: 150px;"),
        h4("Transcripts (", length(transcripts), ")"),
        tags$div(
          style = "border: 1px solid #ddd; 
                  padding: 10px; 
                  background: #f9f9f9; 
                  max-height: 200px; 
                  overflow-y: auto;
                  width: 80%;
                  margin: 0 auto;",
          tags$ol(
            lapply(transcripts, function(t) {
              tags$li(as.character(t))
            })
          )
        )
      )
    })
  # Reduced from 800x1200
  # output$all_gene <- renderPlot({
  #   validate(
  #     need(ncol(mydata()) != "0", "Unrecognized data set: Please
  #         upload your file"))
  #   options(shiny.sanitize.errors = TRUE)
  #   progress <- shiny::Progress$new()
  #   on.exit(progress$close())
  #   progress$set(message = "Please wait a few minutes: Making plot",
  #               detail = 'This may take a while', value = 0)
  #   for (i in 1:40) {
  #     progress$set(message = "Please wait a few minutes: Making plot",
  #                 detail = 'This may take a while', value = i)
  #     Sys.sleep(1.0)
  #   }
  #   Sys.sleep(0.1)
  #   p1()
  # })
  # ============================================================================
  # OUTPUT: ANNOTATED VARIANTS TABLE (DT WIDGET)
  # TABELLA PRINCIPALE CON FORMATTING
  # ============================================================================
  output$uncover_position <- DT::renderDT({
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "Loading annotation table...", value = 0)
    Sys.sleep(0.1)
  
    # RIMUOVI validazione di Gene_name (opzionale nel nuovo UI)
    # Valida SOLO che ci siano dati annotati
    data <- annotated_variants_data()
  
    if(is.null(data) || nrow(data) == 0) {
      showNotification("No annotated variants found. Press 'Calculate Annotations' button.", 
                       type = "warning", duration = 5)
      return(NULL)
    }
  
    cat("Displaying", nrow(data), "annotated variants\n")
  
    # Restituisci il widget DT
    condform_table()
  })



  # ============================================================================
  # DOWNLOAD: ANNOTATED VARIANTS EXCEL (con conditional formatting)
  # ============================================================================
  
  output$downloadData <- downloadHandler(
    filename = function() {
      sample_name <- if (!is.null(input$Sample)) {
        tools::file_path_sans_ext(input$Sample)
      } else {
        "variants"
      }
      paste0('annotated_variants_', sample_name, '_', Sys.Date(), '.xlsx')
    },
    
    content = function(file){
      # ============================================================================
      # WAITER - Show loading screen
      # ============================================================================
      show_uncoverapp_waiter(
        message = "Preparing Excel download...",
        detail = "Formatting annotations and applying styles"
      )
      
      # Ensure waiter closes even if error occurs
      on.exit({
        waiter::waiter_hide()
      })
  
      # USA annotated_variants_data() invece di condform_table()!
      data <- annotated_variants_data()
      
      if (is.null(data) || nrow(data) == 0) {
        showNotification("No data to download!", type = "error")
        return(NULL)
      }
      
      # Rimuovi colonne helper
      data_export <- data
      data_export$highlight_important <- NULL
      data_export$is_omim <- NULL
      data_export$is_pathogenic <- NULL
      
      # Crea workbook
      wb <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb, "Variants")
      
      # Stili
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      yellowStyle <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#FFFFE0")  # Light yellow
      highlighted <- openxlsx::createStyle(fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      
      # Stili per OMIM genes
      omim_light <- openxlsx::createStyle(fgFill = "#E3F2FD")  # azzurro chiaro
      omim_dark <- openxlsx::createStyle(fgFill = "#1976D2", fontColour = "white", 
                                         textDecoration = "bold")  # azzurro scuro + bold
      
      # Scrivi dati
      openxlsx::writeData(wb, "Variants", data_export, headerStyle = hs)
      # Trova indici delle colonne dinamicamente
      col_names <- colnames(data_export)
      col_MutationAssessor <- which(col_names == "MutationAssessor")
      col_M_CAP <- which(col_names == "M_CAP")
      col_CADD_PHED <- which(col_names == "CADD_PHED")
      col_AF_gnomAD <- which(col_names == "AF_gnomAD")
      col_ClinVar <- which(col_names == "ClinVar")
      col_start <- which(col_names == "start")
      col_end <- which(col_names == "end")
      col_GENENAME <- which(col_names == "GENENAME")

      nrows <- nrow(data_export)
      data_export$omim_level <- NULL

      # Conditional formatting: MutationAssessor
      if (length(col_MutationAssessor) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '=="H"', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '=="M"', style = yellowStyle)
                # Use dynamic column reference
        col_letter <- LETTERS[col_MutationAssessor]
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), 
                              rule = paste0('AND($', col_letter, '2<>"H",$', col_letter, '2<>"M")'), 
                              style = posStyle)
              }
              
      # Conditional formatting: M_CAP
      if (length(col_M_CAP) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '=="D"', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '!="D"', style = posStyle)

      }
      
      # Conditional formatting: CADD_PHED
      if (length(col_CADD_PHED) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = ">=20", style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = "<20", style = posStyle)
      }
      
      # Conditional formatting: AF_gnomAD
      if (length(col_AF_gnomAD) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = "<0.01", style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = ">=0.01", style = posStyle)
      }
      
      # Conditional formatting: ClinVar
      if (length(col_ClinVar) > 0) {
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '!="."', style = negStyle)
        openxlsx::conditionalFormatting(wb, "Variants", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '=="."', style = posStyle)
      }
      
      # Highlight important variants (yellow background on start/end)
      if (length(col_start) > 0 && length(col_end) > 0) {
        important_rows <- which(
          grepl("H|M", data_export$MutationAssessor) & 
          data_export$ClinVar != "." & 
          !is.na(data_export$AF_gnomAD) & 
          data_export$AF_gnomAD < 0.01
        )
        
        if (length(important_rows) > 0) {
          openxlsx::addStyle(wb, "Variants", style = highlighted,
                             rows = important_rows + 1, 
                             cols = c(col_start, col_end),
                             gridExpand = TRUE)
        }
      }
      # ============================================
      # Color GENENAME based on OMIM status
      # ============================================
      if (length(col_GENENAME) > 0) {
        # CRITICAL: Usa 'data' (originale) non 'data_export' (già modificato)
        omim_high_rows <- which(data$omim_level == "high")
        omim_medium_rows <- which(data$omim_level == "medium")
        
        # Apply dark blue to high priority (OMIM + pathogenic)
        if (length(omim_high_rows) > 0) {
          openxlsx::addStyle(wb, "Variants", style = omim_dark,
                             rows = omim_high_rows + 1,  # +1 perché riga 1 è header
                             cols = col_GENENAME,
                             gridExpand = TRUE)
        }
        
        # Apply light blue to medium priority (OMIM only)
        if (length(omim_medium_rows) > 0) {
          openxlsx::addStyle(wb, "Variants", style = omim_light,
                             rows = omim_medium_rows + 1,
                             cols = col_GENENAME,
                             gridExpand = TRUE)
        }
      }
      # Auto-size columns
      openxlsx::setColWidths(wb, "Variants", cols = 1:ncol(data_export), widths = "auto")
      
      # Salva
      openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
      
      cat("Excel file saved:", nrow(data_export), "variants\n")
    }
  )

  # ============================================================================
  # OUTPUT: maxAF VALUE
  # ============================================================================
  
  output$maxAF <- renderText({
    maxAF_result <- data()
  
    if (is.null(maxAF_result) || is.na(maxAF_result[[1]])) {
      return("N/A - Check input parameters")
    }
  
    signif(maxAF_result[[1]], 3)
  })
  
  # ============================================================================
  # OUTPUT: maxAF TABLE (DT WIDGET con filtro maxAF)
  # ============================================================================

  output$uncoverPosition <- DT::renderDT({
    # âœ… Remove Gene_name requirement - maxAF works with any filter!
    validate(
      need(!is.null(mydata()) && ncol(mydata()) > 0,
           "Please load your coverage file first")
   )
  
   progress <- shiny::Progress$new()
    on.exit(progress$close())
  progress$set(message = "Filtering by maxAF...", value = 0)
  Sys.sleep(0.1)

  # Check if annotations exist
  base_data <- annotated_variants_data()
  
  if(is.null(base_data) || nrow(base_data) == 0) {
    showNotification(
      "No annotated data available. Please click 'Calculate Annotations' first!", 
      type = "warning",
      duration = 5
    )
    return(NULL)
  }

  cat("Rendering maxAF table with", nrow(base_data), "annotated variants\n")

  # Return the DT widget
  uncover_maxaf()
})
  # ============================================================================
  # DOWNLOAD: maxAF EXCEL
  # ============================================================================
  
  output$download_maxAF <- downloadHandler(
    filename = function() {
      paste0('uncover_maxAF_', Sys.Date(), '.xlsx')
    },
    content = function(file){
      # USA il dataframe pulito, NON il widget!
      data <- uncover_maxaf_data()
      
      if (is.null(data) || nrow(data) == 0) {
        showNotification("No data to download!", type = "error")
        return(NULL)
      }
      
      # Rimuovi colonne helper
      data_export <- data
      data_export$highlight_maxaf <- NULL
      data_export$highlight_important <- NULL
      
      wb1 <- openxlsx::createWorkbook()
      openxlsx::addWorksheet(wb1, "MaxAF")
      
      negStyle <- openxlsx::createStyle(fontColour = "#9C0006", bgFill = "#FFC7CE")
      posStyle <- openxlsx::createStyle(fontColour = "#006100", bgFill = "#C6EFCE")
      highlighted <- openxlsx::createStyle(fgFill = "yellow")
      hs <- openxlsx::createStyle(textDecoration = "BOLD", fontColour = "#FFFFFF",
                        fontSize=12,
                        fontName="Arial Narrow", fgFill = "#4F80BD")
      
      openxlsx::writeData(wb1, "MaxAF", data_export, headerStyle = hs)
      
      # Trova indici colonne dinamicamente
      col_names <- colnames(data_export)
      col_MutationAssessor <- which(col_names == "MutationAssessor")
      col_M_CAP <- which(col_names == "M_CAP")
      col_CADD_PHED <- which(col_names == "CADD_PHED")
      col_AF_gnomAD <- which(col_names == "AF_gnomAD")
      col_ClinVar <- which(col_names == "ClinVar")
      col_start <- which(col_names == "start")
      col_end <- which(col_names == "end")
      
      nrows <- nrow(data_export)
      
      # Conditional formatting (stesso pattern di downloadData)
      if (length(col_MutationAssessor) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '=="H"', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_MutationAssessor,
                              rows = 2:(nrows + 1), rule = '!="H"', style = posStyle)
      }
      
      if (length(col_M_CAP) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '=="D"', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_M_CAP,
                              rows = 2:(nrows + 1), rule = '!="D"', style = posStyle)
      }
      
      if (length(col_CADD_PHED) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = ">=20", style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_CADD_PHED,
                              rows = 2:(nrows + 1), rule = "<20", style = posStyle)
      }
      
      if (length(col_AF_gnomAD) > 0) {
        # Per maxAF usiamo threshold 0.05 invece di 0.5
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = "<=0.05", style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_AF_gnomAD,
                              rows = 2:(nrows + 1), rule = ">0.05", style = posStyle)
      }
      
      if (length(col_ClinVar) > 0) {
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '!="."', style = negStyle)
        openxlsx::conditionalFormatting(wb1, "MaxAF", cols = col_ClinVar,
                              rows = 2:(nrows + 1), rule = '=="."', style = posStyle)
      }
      
      # Highlight important variants
      if (length(col_start) > 0 && length(col_end) > 0) {
        maxAF_value <- signif(data()[[1]], 3)
        
        important_rows <- which(
          grepl("H|M", data_export$MutationAssessor) & 
          data_export$ClinVar != "." & 
          !is.na(data_export$AF_gnomAD) & 
          data_export$AF_gnomAD < maxAF_value
        )
        
        if (length(important_rows) > 0) {
          openxlsx::addStyle(wb1, "MaxAF", style = highlighted,
                             rows = important_rows + 1, 
                             cols = c(col_start, col_end),
                             gridExpand = TRUE)
        }
      }
      
      openxlsx::setColWidths(wb1, "MaxAF", cols = 1:ncol(data_export), widths = "auto")
      openxlsx::saveWorkbook(wb1, file, overwrite = TRUE)
      
      cat("MaxAF Excel saved:", nrow(data_export), "variants\n")
    }
  )
  
  
  # ============================================================================
  # CLOSE APP BUTTON
  # ============================================================================
  
  observe({
    if (input$close > 0) stopApp()
  })
}