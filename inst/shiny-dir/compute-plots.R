
#strack<- reactive({if (input$UCSC_Genome == "hg19"){
#  Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#  else{
#    Gviz::SequenceTrack(BSgenome.Hsapiens.UCSC.hg19::Hsapiens, chromosome = Chromosome())}
#})
dtrack1 <- reactive({
  cat("\n=== DTRACK1 ===\n")
  grcoverage <- filtered_low()
  cat("Rows:", nrow(grcoverage), "\n")
  cat("Columns:", paste(colnames(grcoverage), collapse=", "), "\n")
  cat("First 3 rows:\n")
  print(head(grcoverage, 3))
  
  dt1 <- Gviz::DataTrack(range=grcoverage, type="histogram", 
                         fill.histogram="red", col.histogram="NA",
                         genome=input$UCSC_Genome, name="Seq. Depth")
  cat("DataTrack created OK\n")
  return(dt1)
})

dtrackHigh <- reactive({
  cat("\n=== DTRACKHIGH ===\n")
  
  grcoverage_high <- tryCatch({
    filtered_high()
  }, error = function(e) {
    cat("ERROR calling filtered_high():", e$message, "\n")
    return(NULL)
  })
  
  cat("filtered_high returned\n")
  
  if (is.null(grcoverage_high)) {
    cat("ERROR: filtered_high is NULL\n")
    return(NULL)
  }
  
  cat("Rows:", nrow(grcoverage_high), "\n")
  cat("Columns:", paste(colnames(grcoverage_high), collapse=", "), "\n")
  
  if (nrow(grcoverage_high) == 0) {
    cat("WARNING: filtered_high is empty (no high coverage regions)\n")
    # Ritorna un DataTrack vuoto invece di NULL
    dt2 <- Gviz::DataTrack(
      data.frame(chromosome=character(), start=integer(), end=integer(), coverage=integer()),
      type="histogram", 
      fill.histogram="dodgerblue", 
      col.histogram="NA", 
      genome=input$UCSC_Genome, 
      name="Seq. Depth"
    )
    return(dt2)
  }
  
  cat("First 3 rows:\n")
  print(head(grcoverage_high, 3))
  
  dt2 <- Gviz::DataTrack(range=grcoverage_high, type="histogram", 
                         fill.histogram="dodgerblue", col.histogram="NA", 
                         genome=input$UCSC_Genome, name="Seq. Depth")
  cat("DataTrack HIGH created OK\n")
  return(dt2)
})

itrack <- reactive ({
  i= Gviz::IdeogramTrack(genome =input$UCSC_Genome , chromosome = Chromosome())
})

grtrack  <-reactive({
  ggT<-Gviz::GeneRegionTrack (txdb(), chromosome  = Chromosome() ,  start  =  St() ,  end  =  En(),
                        showId  =  TRUE , transcriptAnnotation="symbol",
                        name  =  " Gene Annotation ")
  return(ggT)
})

#Preparation for all gene coverage plot
p1 <- reactive({
  cat("\n=== P1 START ===\n")
  
  # Check dependencies
  cat("Checking coord()...\n")
  if (is.null(coord())) {
    cat("ERROR: coord() is NULL\n")
    return(NULL)
  }
  
  cat("Checking Chromosome()...\n")
  chr <- Chromosome()
  cat("Chromosome:", chr, "\n")
  
  if (is.null(chr) || chr == "") {
    cat("ERROR: Chromosome is empty\n")
    return(NULL)
  }
  
  cat("Getting start_gene...\n")
  start_gene <- coord()$start[coord()$seqnames == chr]
  if (length(start_gene) > 1)
    start_gene <- start_gene[1]
  print(start_gene)
  
  cat("Getting end_gene...\n")
  end_gene <- coord()$end[coord()$seqnames == chr]
  if (length(end_gene) > 1)
    end_gene <- end_gene[1]
  print(end_gene)
  cat("Creating tracks...\n")
  ot <- Gviz::OverlayTrack(trackList = list(dtrack1(), dtrackHigh()))
  cat("Overlay track OK\n")
  
  gtrack <- Gviz::GenomeAxisTrack()
  cat("Genome axis track OK\n")
  
  ylims <- grDevices::extendrange(range(c(dtrack1()@data), dtrackHigh()@data))
  cat("Y limits OK\n")
  
  cat("Creating itrack...\n")
  it <- itrack()
  cat("itrack OK\n")
  
  cat("Creating gene region track...\n")
  gr_ex_track <- Gviz::GeneRegionTrack(
    txdb(),
    chromosome = chr,
    start = start_gene, 
    end = end_gene,
    showId = TRUE,
    name = "Gene Annotation"
  )
  cat("Gene region track OK\n")
  cat("Calling plotTracks...\n")
  Gviz::plotTracks(
    list(it, gtrack, ot, gr_ex_track),
    from = start_gene,
    to = end_gene,
    reverseStrand = TRUE,
    ylim = ylims,
    type = "histogram",
    baseline = input$coverage_co,
    col.baseline = "black",
    lwd.baseline = 0.3,
    extend.left = 0.5, 
    extend.right = 200
  )
  cat("=== P1 COMPLETE ===\n")
})

#table of uncovered exons

# table1<- reactive({
#   if (is.null(filtered_low()))
#     return(NULL)
#   f.low=filtered_low()
#   f.low[,'new']='NA'
#   f.low$new<- ifelse(sapply(f.low$start, function(p)
#     any(exon_gp()$start <= p & exon_gp()$end >= p)), "YES", "out")
#   l.coverage= as.data.frame(f.low[f.low$new=='YES',])
#   validate(
#     need(nrow(l.coverage) >0, "ALL EXONS ARE COVERED
#          UNDER YOUR CHOOSE THRESHOLD"))

#   x= l.coverage$start
#   getValue3 <- function(x, data) {
#     tmp <- data %>%
#       dplyr::filter(start <= x, x <= end) %>%
#       dplyr::filter(number_of_transcript == input$transcript_id)
#     return(tmp$exon_rank)
#   }

#   a_exon=sapply(x, getValue3, data=exon_gp())
#   exon=unlist(lapply(a_exon, function (x) ifelse(length (x) > 0, x, NA)))
#   exon.df= as.data.frame(exon)
#   if (is.null(exon.df))
#     return(NULL)
#   df.l= cbind(l.coverage, exon.df)
#   t1=table(df.l$exon)
#   df.t1= as.data.frame(t1)
#   colnames(df.t1)= c('exon','uncovered positions')
#   return(df.t1)
# })

#output$df.l<- renderDataTable({
 # table1()
#})



observeEvent(input$btn_run,{
  Sys.sleep(5)
})

