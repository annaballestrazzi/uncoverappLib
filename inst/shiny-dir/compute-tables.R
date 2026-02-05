#Tables
tryObserve <- function(x) {
  x <- substitute(x)
  env <- parent.frame()
  observe({
    tryCatch(eval(x, env),
             error = function(e) {
               #showNotification(paste("Error: ", e$message), type = "error")
             })
  })
}


coord= eventReactive(input$ucsc_lookup,{
  cat("\n=== COORD() REACTIVE TRIGGERED ===\n")
  
  disable("ucsc_lookup")
  shinyjs::show("text1")
  Sys.sleep(0.1)
  
  cat("Looking up gene:", input$Gene_name, "\n")
  
  my_gene_name=OrganismDbi::select(org.Hs.eg.db,
                                   key= input$Gene_name,
                                   columns=c("ENTREZID","GENENAME", "ENSEMBL"),
                                   keytype="SYMBOL")
  
  ID=my_gene_name$ENTREZID
  
  if (is.null(ID) || length(ID) == 0) {
    cat("ERROR: No ENTREZID found for gene:", input$Gene_name, "\n")
    shinyjs::enable("ucsc_lookup")
    shinyjs::hide("text1")
    return(NULL)
  }
  
  cat("Found ENTREZID:", ID, "\n")

  gene_coords <- AnnotationDbi::select(
    txdb(),
    keys = ID,
    columns = c("TXCHROM", "TXSTART", "TXEND", "TXNAME"),
    keytype = "GENEID"
  )

  validate(
    need(!is.null(gene_coords) && nrow(gene_coords) > 0,
        paste("Gene not found in", input$UCSC_Genome, "database"))
  )

  # ============================================================================
  # FILTER OUT ALTERNATIVE CONTIGS (chr22_fix, chr1_random, etc.)
  # Same method used in compute-plots.R and compute-preprocess.R
  # ============================================================================
  cat("Raw transcripts:", nrow(gene_coords), "\n")
  cat("Chromosomes found:", paste(unique(gene_coords$TXCHROM), collapse=", "), "\n")
  
  # Filter out any chromosome with underscore (alternative contigs)
  gene_coords_clean <- gene_coords %>%
    dplyr::filter(!grepl("_", TXCHROM))
  
  cat("After filtering alt contigs:", nrow(gene_coords_clean), "transcripts\n")
  
  validate(
    need(nrow(gene_coords_clean) > 0,
        "All transcripts are on alternative contigs - cannot use this gene")
  )
  
  cat("Clean coordinate range:", min(gene_coords_clean$TXSTART), "-", 
      max(gene_coords_clean$TXEND), "\n\n")
  
  # ============================================================================
  

  info <- gene_coords_clean %>%
    dplyr::group_by(GENEID) %>%
    dplyr::summarise(
      seqnames = unique(TXCHROM)[1],
      start = min(TXSTART, na.rm = TRUE),
      end = max(TXEND, na.rm = TRUE),
      transcript_ids = list(unique(TXNAME)),
      .groups = "drop"
    ) %>%
    dplyr::rename(ENTREZID = GENEID) %>%
    dplyr::inner_join(my_gene_name, by = "ENTREZID") %>%
    dplyr::slice(1)  # ✅ Take only first row (removes duplicates from multiple ENSEMBL IDs)

  cat("Gene coordinates:\n")
  print(info[, c("seqnames", "start", "end", "SYMBOL")])
    
    shinyjs::enable("ucsc_lookup")
    shinyjs::hide("text1")
    
    cat("=== COORD() COMPLETE ===\n\n")
    
    return(info)
  })
  # ============================================================================
# GLOBAL GENE COORDINATES - Calculated once, used everywhere
# ============================================================================
gene_coordinates <- reactive({
  # Require that coord() has been called
  coord_data <- coord()
  
  if (is.null(coord_data) || nrow(coord_data) == 0) {
    return(NULL)
  }
  
  list(
    chr = as.character(coord_data$seqnames[1]),
    start = coord_data$start[1],
    end = coord_data$end[1],
    symbol = coord_data$SYMBOL[1]
  )
})


observeEvent(input$ucsc_lookup, {
  output$ccg <-DT::renderDataTable({
    progress <- shiny::Progress$new()
    progress$set(message = "Running", detail = 'This may take a while')
    on.exit(progress$close())
    validate(
      need(input$Gene_name != "",
        "Unrecognized gene name: Please select HGNC gene name \n Click apply"))
    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    coord()


  })
})
# Auto-update Chromosome input when lookup button is pressed
observeEvent(input$ucsc_lookup, {
  cat("\n=== AUTO-UPDATE CHROMOSOME OBSERVER TRIGGERED ===\n")
  
  # Wait a moment for coord() to compute
  Sys.sleep(0.2)
  
  coord_data <- coord()
  
  if (is.null(coord_data)) {
    cat("WARNING: coord() returned NULL, cannot auto-update chromosome\n")
    return()
  }
  
  if (nrow(coord_data) == 0) {
    cat("WARNING: coord() returned empty dataframe\n")
    return()
  }
  
  x <- as.data.frame(coord_data)
  Chrom <- as.character(x$seqnames[1])
  
  cat("Auto-updating Chromosome input to:", Chrom, "\n")
  
  updateTextInput(session, "Chromosome", value = Chrom)
  
  cat("Chromosome input updated successfully\n\n")
})


Chromosome<- reactive({
  xc= as.character(input$Chromosome)
  return(xc)
})


observeEvent(input$remove,{
  shinyjs::hide(coord())
  output$ccg<- NULL
})



exon_gp<-eventReactive(input$ucsc_lookup,{
  #require(rtracklayer)
  #ucsc <- browserSession()
  gname =input$Gene_name
  if (is.null(gname))
    return(NULL)
  else if (gname== "")
    return(NULL)
  if (input$UCSC_Genome == "hg19"){
    edb = EnsDb.Hsapiens.v75}
  else{
    edb= EnsDb.Hsapiens.v86}
  eid <- OrganismDbi::select(org.Hs.eg.db,gname,
                             "ENTREZID", "SYMBOL")[["ENTREZID"]]
  txid <- OrganismDbi::select(txdb(), eid, "TXNAME", "GENEID")[["TXNAME"]]
  cds <- cdsBy(txdb(), by="tx", use.names=TRUE)
  exoncds <- cds[names(cds) %in% txid]
  exon_Id<- as.data.frame(exoncds)
  exon_table=exon_Id  %>%
    dplyr::select(1:6,8,10)
  #exon_table=exon_Id[c(1:6,8,10)]
  colnames(exon_table)=c("number_of_transcript",
                         "type_of_transcript", "chrom", "start","end",
                         "length_of_exon", "cds_id", "exon_rank")
  return(exon_table)
  print(head(exon_table))
})


observeEvent(input$ucsc_lookup, {
  output$exon_pos<- DT::renderDataTable({
    options(shiny.sanitize.errors = TRUE)
    progress <- shiny::Progress$new()
    on.exit(progress$close())
    progress$set(message = "table construction in progress",
                 detail = 'This may take a while', value = 0)
    Sys.sleep(0.1)


    HGNC_org <- keys(org.Hs.eg.db, keytype = "SYMBOL")
    validate(
      need(HGNC_org[HGNC_org %in% input$Gene_name],
           "incorrect HGNC Gene symbol"))
    exon_gp()
  })
})




observeEvent(input$remove,{
  output$exon_pos<- DT::renderDataTable({
    shinyjs::js$reset(exon_gp()) })
})



tryObserve({
  if (is.null(input$exon_number) )
    return()
  message("load inputs")
  #observe({
  # if (is.null(input$exon_number))
  #  return(NULL)
  if (is.null(exon_gp()))
    return(NULL)
  exon_df= as.data.frame(exon_gp())
  one=subset(exon_df,
             exon_df$exon_rank == input$exon_number &
               exon_df$number_of_transcript == input$transcript_id)
  print(one)
  start_exon= as.numeric(one$start)
  end_exon= as.numeric(one$end)
  chr= as.numeric(gsub("\\D", "", one$chrom, perl =TRUE))
  trascript_name= as.character(one$type_of_transcript)
  updateTextInput(session, "id_t", value = trascript_name)
  updateTextInput(session, "Start_genomicPosition", value = start_exon)
  updateTextInput(session, "end_genomicPosition", value = end_exon)
  updateTextInput(session, "query_Database",
                  value = paste0(chr,":",start_exon,"-",end_exon))
})


#I define reactive START and END position


St<- reactive({
  sta= as.numeric(input$Start_genomicPosition)
  return(sta)
})


En<- reactive({
  en= as.numeric(input$end_genomicPosition)
  return(en)
})


id<- reactive({
  id_ucsc= as.character(input$id_t)
  return(id_ucsc)
})


coordinate<- reactive({
  coo=paste0(chr,":",start_exon,"-",end_exon)
  return(coo)
})