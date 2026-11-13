suppressPackageStartupMessages({
        require(OrganismDbi)
        require(TxDb.Hsapiens.UCSC.hg19.knownGene)
        require(TxDb.Hsapiens.UCSC.hg38.knownGene)
        require(org.Hs.eg.db)
        require(rlist)
        require(Rsamtools)
    })
    cat("Heavy packages loaded!\n")

ff <- reactiveValues()

gene_list <- reactive({
   start_time <- Sys.time()
   file_gene <- input$gene1
   if (is.null(input$gene1)) return(NULL)
   
   if(input$type_file == "target_bed") {
      tmp_gene1 <- read.table(input$gene1$datapath, stringsAsFactors = F)
      tmp_gene <- tmp_gene1[1:4]
      colnames(tmp_gene) <- c('chr', 'start', 'end', 'SYMBOL')
      ff <- tmp_gene
   } else {
      tmp_gene <- scan(file_gene$datapath, character(), quote = "")
      ff <- tmp_gene
   }
    cat(paste("Load input files:", round(difftime(Sys.time(), start_time, units="mins"), 2), "minutes\n\n"))
    return(ff)
})

list_coverage <- reactive({
    if (is.null(input$bam_list)) return(NULL)
    scan(input$bam_list$datapath, character(), quote = "")
})

all_gene <- reactive({
   if (input$Genome == "hg19") {
      TxDb.Hsapiens.UCSC.hg19.knownGene
   } else {
      TxDb.Hsapiens.UCSC.hg38.knownGene
   }
})

# ============================================================================
# NUOVO REACTIVE CENTRALIZZATO: calcola UNA SOLA VOLTA la validazione dei geni
# ============================================================================
validated_genes_data <- reactive({
   cat("\n=== RUNNING ON R SESSION:", Sys.info()["nodename"], "===\n")
   cat("R Process ID:", Sys.getpid(), "\n")
   cat("Working directory:", getwd(), "\n\n")

   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   if (input$type_file == "target_bed") return(NULL)  # Non serve per bed files
   
   cat("\n=== GENE VALIDATION (CENTRALIZED) ===\n")
   
   # STAGE 1: Convert gene names to ENTREZ IDs (UNA SOLA VOLTA!)
   cat("Stage 1: Converting gene names to ENTREZ IDs...\n")
   
   alias_tab <- OrganismDbi::select(
      org.Hs.eg.db,
      key = gene_list(),
      columns = c("ENTREZID", "SYMBOL"),
      keytype = "ALIAS"
   ) %>%
      tibble::as_tibble() %>%
      dplyr::filter(!is.na(ENTREZID))

   # A name can be the official symbol of one gene and an alias of other genes
   # (e.g. ARX is also an alias of UBA2): the gene with that official symbol wins.
   my_gene_name <- alias_tab %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::arrange(dplyr::desc(SYMBOL == ALIAS), .by_group = TRUE) %>%
      dplyr::slice(1) %>%
      dplyr::ungroup()

   # A name that is only an alias, of more than one gene, is ambiguous: warn
   ambiguous <- alias_tab %>%
      dplyr::group_by(ALIAS) %>%
      dplyr::filter(!any(SYMBOL == ALIAS, na.rm = TRUE),
                    dplyr::n_distinct(ENTREZID) > 1) %>%
      dplyr::summarise(others = paste(unique(SYMBOL), collapse = "/"), .groups = "drop") %>%
      dplyr::inner_join(my_gene_name, by = "ALIAS")
   if (nrow(ambiguous) > 0) {
      cat("WARNING: Ambiguous gene names (alias of more than one gene; write the official symbol to choose):\n")
      cat(paste0("  - ", ambiguous$ALIAS, " -> ", ambiguous$SYMBOL, " used (candidates: ",
                 ambiguous$others, ")\n"), sep = "")
   }
   renamed <- my_gene_name %>% dplyr::filter(SYMBOL != ALIAS)
   if (nrow(renamed) > 0) {
      cat("Aliases converted to official symbol:",
          paste0(renamed$ALIAS, " -> ", renamed$SYMBOL, collapse = ", "), "\n")
   }
   my_gene_name <- my_gene_name %>% dplyr::select(ALIAS, ENTREZID)
   found_genes_stage1 <- unique(my_gene_name$ALIAS)
   not_found_stage1 <- base::setdiff(gene_list(), found_genes_stage1)
   
   if (length(not_found_stage1) > 0) {
      cat("WARNING: Genes not found in org.Hs.eg.db:\n")
      for(gene in not_found_stage1) {
         cat(paste("  - Gene '", gene, "'\n", sep=""))
      }
      cat("\n")
   }
   
    validate(
    need(nrow(my_gene_name) > 0, 
        "No valid gene found in org.Hs.eg.db. Verify the entered gene names.")
    )
   
   cat(paste("Stage 1 result:", nrow(my_gene_name), "mappings from", 
             length(found_genes_stage1), "genes\n\n"))
   
   # STAGE 2: Validate against genome (UNA SOLA VOLTA!)
   cat("Stage 2: Validating against reference genome...\n")
   
   all_entrez_ids <- unique(my_gene_name$ENTREZID)
   genome_keys <- AnnotationDbi::keys(all_gene(), "GENEID")  # CACHE THIS!
   valid_entrez_ids <- base::intersect(all_entrez_ids, genome_keys)
   invalid_entrez_ids <- base::setdiff(all_entrez_ids, genome_keys)
   
   invalid_genes_stage2 <- character(0)
   if (length(invalid_entrez_ids) > 0) {
      invalid_genes_stage2 <- my_gene_name %>% 
         dplyr::filter(ENTREZID %in% invalid_entrez_ids) %>%
         dplyr::pull(ALIAS) %>%
         unique()
      
      cat("WARNING: Genes not in reference genome:\n")
      for(gene in invalid_genes_stage2) {
         cat(paste("  - Gene '", gene, "' (", input$Genome, ")\n", sep=""))
      }
      cat("\n")
   }
   
   # Filter to valid genes only
   my_gene_name <- my_gene_name %>%
      dplyr::filter(ENTREZID %in% valid_entrez_ids)
   
   validate(
        need(length(valid_entrez_ids) > 0, 
        "No gene found in the selected reference genome. Try changing the genome (hg19/hg38).")
   )
   
   cat(paste("Stage 2 result:", length(valid_entrez_ids), "validated genes\n"))
   
   total_excluded <- length(not_found_stage1) + length(invalid_genes_stage2)
   if (total_excluded > 0) {
      cat(paste("SUMMARY: ", total_excluded, " genes excluded (", 
                length(not_found_stage1), " stage 1, ", 
                length(invalid_genes_stage2), " stage 2)\n", sep=""))
   }
   cat("=== VALIDATION COMPLETE ===\n\n")
   list(
      my_gene_name = my_gene_name,
      valid_entrez_ids = valid_entrez_ids,
      genome_keys = genome_keys,  # Cache per riuso
      not_found_stage1 = not_found_stage1,
      invalid_genes_stage2 = invalid_genes_stage2
   )

})

# ============================================================================
# SEMPLIFICATO: usa solo i risultati pre-calcolati
# ============================================================================
no_entrID <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   
   if (input$type_file == "target_bed") {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Usa i dati già validati!
   validated <- validated_genes_data()
   
   if (is.null(validated)) {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Se ci sono geni validi, ritorna dataframe vuoto (nessun errore)
   if (length(validated$valid_entrez_ids) > 0) {
      return(data.frame(matrix(ncol = 2, nrow = 0)))
   }
   
   # Altrimenti ritorna errori
   errore <- data.frame(
      Gene = c("No valid genes found"),
      Issue = c("All genes failed validation"),
      stringsAsFactors = FALSE
   )
   return(errore)
})

# ============================================================================
# OTTIMIZZATO: usa dati pre-calcolati + ottimizzazioni loop
# ============================================================================

for_bed <- reactive({
   if (is.null(gene_list())) return(NULL)
   if (is.null(list_coverage())) return(NULL)
   if (nrow(no_entrID()) != 0) return(no_entrID())
   
   if (input$type_file == "target_bed") {
      bed_data <- gene_list()
      cat("Target BED input:", nrow(bed_data), "regions\n")
      
      # Filter alternative contigs from target BED
      cat("Before alt contig filter:", nrow(bed_data), "intervals\n")
      bed_data <- bed_data %>% dplyr::filter(!grepl("_", chr))
      cat("After alt contig filter:", nrow(bed_data), "intervals (removed alt contigs)\n")
      
      validate(
        need(nrow(bed_data) > 0, 
            "all regions target are on alternative contigs. The BED file must contain regions on the main chromosomes (chr1-22, X, Y).")
    )
      
      return(bed_data)
   }

   # USA I DATI GIÀ VALIDATI - zero duplicazione!
   validated <- validated_genes_data()
   if (is.null(validated)) return(NULL)
   
   my_gene_name <- validated$my_gene_name
   ID <- validated$valid_entrez_ids
   # genome_keys già disponibile ma non serve più (vedi sotto)
   
   # OTTIMIZZAZIONE: carica exonsBy UNA SOLA VOLTA (era dentro il loop!)
   cds <- OrganismDbi::exonsBy(all_gene(), by = "tx", use.names = TRUE)
   txid_map <- OrganismDbi::select(all_gene(), ID, "TXNAME", "GENEID")
   
   # OTTIMIZZAZIONE: usa lapply invece di loop con rbind
    pre_list<- lapply(split(txid_map, txid_map$GENEID), function(gene_txs) {
        tx_names <- gene_txs$TXNAME
        coor <- as.data.frame(cds[names(cds) %in% tx_names])
        if (nrow(coor) == 0) return(NULL)
        coor$ENTREZID <- gene_txs$GENEID[1]
        return(coor)
        })

        pre_list <- pre_list[!sapply(pre_list, is.null)]
    
    cat("DEBUG: pre_list has", length(pre_list), "non-NULL elements\n")
    
    # Check se ci sono dati
    validate(
        need(length(pre_list) > 0, 
        "No coordinates found for the entered genes. Verify that the genes exist in the selected genome.")
    )
   
   # Combina tutto in una volta (molto più veloce di rbind iterativo)
   pre <- dplyr::bind_rows(pre_list) %>%
      dplyr::inner_join(my_gene_name, by = "ENTREZID") %>%
      dplyr::mutate(
         start = start - 10,
         end = end + 10
      )
   
   # OTTIMIZZAZIONE: rimuovi alt contigs PRIMA di creare il BED
   for_bed <- pre %>%
      dplyr::transmute(
         chr = as.character(seqnames),
         start = start,
         end = end,
         SYMBOL = ALIAS
      ) %>%
      dplyr::distinct()
    cat("Before alt contig filter:", nrow(for_bed), "intervals\n")
    for_bed <- for_bed %>%
      dplyr::filter(!grepl("_", chr))
    cat("After alt contig filter:", nrow(for_bed), "intervals (removed alt contigs)\n")
   
   # Apply notation transformation
   if (input$notation == "number") {
      for_bed <- for_bed %>%
         dplyr::mutate(chr = sub("^.{0,3}", "", chr, perl = TRUE))
   }
   return(for_bed)
})

# ============================================================================
# COMPLETELY FIXED COVERAGE PROCESSING
# Handles both 0-based and 1-based coordinate systems
# Coverage BED format: chr start end depth (uniform depth per interval)
# ============================================================================
coverage_input <- eventReactive(input$process_coverage, {
    cat("\n=== PROCESS COVERAGE BUTTON PRESSED ===\n")
    
    # Show waiter FIRST
    show_uncoverapp_waiter(
      message = "Processing coverage files...",
      detail = "This may take several minutes for large datasets"
    )
    
    # Ensure waiter is hidden even if error occurs
    on.exit(waiter::waiter_hide(), add = TRUE)
    
    # Validazione iniziale
    
    # Validazione iniziale
    if (is.null(gene_list())) {
        showNotification("Please upload a gene list file first!", 
                        type = "warning", duration = 5)
        return(NULL)
    }
    
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) {
        return(no_entrID())
    }
    
    if (is.null(list_coverage())) {
        showNotification("Please upload a BAM/BED list file first!", 
                        type = "warning", duration = 5)
        return(NULL)
    }

    
    cat("\n=== COVERAGE PROCESSING ===\n")
    
    # Determine input coordinate system (add this as UI input!)
    # For now, default to 0-based (standard BED format)
    input_coord_system <- if (!is.null(input$input_coordinate_system)) {
        input$input_coordinate_system
    } else {
        "0-based"  # Default assumption
    }
    cat(paste("Input coordinate system:", input_coord_system, "\n"))
    
    # Get target regions and convert to 1-based for GenomicRanges
    target_bed <- for_bed()
    # User target BED in 0-based coordinates: convert to 1-based like the BED coverage
    # files (same rule as buildInput). Gene-list regions from UCSC are already 1-based.
    if (input$type_file == "target_bed" && input$type_coverage == "bed" &&
        input_coord_system == "0-based") {
        target_bed$start <- target_bed$start + 1
    }
    for_grange <- GenomicRanges::makeGRangesFromDataFrame(
        target_bed, 
        keep.extra.columns = TRUE
    )
    # Merge overlapping intervals of the same gene (e.g. exons of alternative
    # transcripts) so that every base is counted once per gene
    for_grange <- unlist(GenomicRanges::reduce(split(for_grange, for_grange$SYMBOL)))
    for_grange$SYMBOL <- names(for_grange)
    names(for_grange) <- NULL
    if (input$type_coverage == "bam") {
            cat("Processing BAM files with pileup...\n\n")
    
            idx_stats <- tryCatch({
            # Read the BAM index (.bai) with Rsamtools: same information as "samtools idxstats"
            result <- as.data.frame(Rsamtools::idxstatsBam(list_coverage()[1]))
            result
            }, error = function(e) {
            waiter::waiter_hide()
            showNotification(
                paste("❌ Could not read BAM file:", basename(list_coverage()[1]),
                    "\nVerify that paths in your list file are correct and files are indexed."),
                type = "error", duration = 10
            )
            return(NULL)
            })

            if (is.null(idx_stats)) return(NULL)

            idx_df <- idx_stats
            colnames(idx_df) <- c("seqnames", "length", "mapped", "unmapped")
            idx_df$seqnames <- as.character(idx_df$seqnames)
            chrs_with_data <- idx_df$seqnames[idx_df$mapped > 0]
            target_chrs <- as.character(unique(GenomicRanges::seqnames(for_grange)))
            available_chrs <- intersect(target_chrs, chrs_with_data)

    
            cat(paste("Chromosomes in BAM header:", nrow(idx_df), "\n"))
            cat(paste("Chromosomes with data:", length(chrs_with_data), "\n"))
            cat(paste("Target chromosomes:", length(target_chrs), "\n"))
            cat(paste("Available overlap:", length(available_chrs), "\n\n"))
    
            validate(
                need(length(available_chrs) > 0, 
                    "No chromosome overlap between BAM files and requested genes. Verify that the BAM files have data for the selected genes.")
            )
    
            for_grange_filtered <- for_grange[
                as.character(GenomicRanges::seqnames(for_grange)) %in% available_chrs
            ]
    
            # ORA USA for_grange_filtered:
            # reduce(): regions of different genes may overlap, each base is read once
            # Reads not counted (as samtools depth and mosdepth): PCR/optical duplicates,
            # secondary alignments and reads failing quality checks. Supplementary are counted.
            param <- Rsamtools::ScanBamParam(
                which = GenomicRanges::reduce(for_grange_filtered),
                flag = Rsamtools::scanBamFlag(isDuplicate = FALSE,
                                              isSecondaryAlignment = FALSE,
                                              isNotPassingQualityControls = FALSE)
            )
            p_param <- Rsamtools::PileupParam(
            distinguish_nucleotides = FALSE,
            distinguish_strands = FALSE,
            min_mapq = as.numeric(input$MAPQ),
            min_base_quality = as.numeric(input$base_qual),
            min_nucleotide_depth = 1,
            max_depth = 150000
        )
        
        lst3 <- lapply(list_coverage(), function(i) {
            pileup_df <- Rsamtools::pileup(i, scanBamParam = param, pileupParam = p_param)
            pileup_df <- pileup_df[, c("seqnames", "pos", "count")]
            pileup_df <- pileup_df[!duplicated(pileup_df), ]
            
            # Pileup returns 1-based positions
            pileup_df <- pileup_df %>%
                dplyr::mutate(end = pos) %>%
                dplyr::group_by(seqnames, pos, end) %>%
                dplyr::summarise(count = sum(count), .groups = "drop")

            # Overlapping mates of a pair (same fragment): count the overlapping bases once
            mate_overlap <- if (is.null(input$mate_overlap)) "all" else input$mate_overlap
            if (mate_overlap != "none") {
                ov <- mate_overlap_counts(i, param, as.numeric(input$MAPQ), as.numeric(input$base_qual),
                                          proper_only = (mate_overlap == "proper"))
                pileup_df$seqnames <- as.character(pileup_df$seqnames)
                pileup_df <- dplyr::left_join(pileup_df, ov, by = c("seqnames", "pos"))
                pileup_df$count <- pileup_df$count - ifelse(is.na(pileup_df$n_overlap), 0L, pileup_df$n_overlap)
                pileup_df$n_overlap <- NULL
            }
            pileup_df
        })
        
        # All target bases, each once; a base without reads in a sample gets 0
        all_bases <- GenomicRanges::reduce(GenomicRanges::granges(for_grange))
        w <- GenomicRanges::width(all_bases)
        pp <- data.frame(
          seqnames = rep(as.character(GenomicRanges::seqnames(all_bases)), w),
          pos = rep(GenomicRanges::start(all_bases), w) + sequence(w) - 1L,
          stringsAsFactors = FALSE
        )
        pp$end <- pp$pos
        for (k in seq_along(lst3)) {
          s <- as.data.frame(lst3[[k]])[, c("seqnames", "pos", "count")]
          s$seqnames <- as.character(s$seqnames)
          colnames(s)[3] <- paste0("count", k)
          pp <- merge(pp, s, by = c("seqnames", "pos"), all.x = TRUE)
        }
        pp[is.na(pp)] <- 0
        pp <- pp[, c("seqnames", "pos", "end", paste0("count", seq_along(lst3)))]
        colnames(pp)[1:3] <- c("seqnames", "start", "end")
        # Fix column names: use actual sample names with count_ prefix
        sample_names <- tools::file_path_sans_ext(basename(list_coverage()))
        cat(paste("Merged result:", nrow(pp), "unique intervals\n"))
        cat("NA count after merge:", sum(is.na(pp)), "\n")
        cat("Zero count after merge:", sum(pp == 0, na.rm = TRUE), "\n")
        # Check for duplicates and fix
        duplicates <- sample_names[duplicated(sample_names)]
        if (length(duplicates) > 0) {
            sample_names <- make.unique(sample_names, sep = "_")
            showNotification(
                paste("Duplicate sample names detected:", 
                    paste(unique(duplicates), collapse=", "),
                    "- automatically renamed to avoid conflicts."),
                type = "warning",
                duration = 10
            )
            cat("WARNING: Duplicate sample names renamed to:", paste(sample_names, collapse=", "), "\n")
        }

        count_cols <- colnames(pp)[-(1:3)]

        if (length(sample_names) == length(count_cols)) {
            colnames(pp)[4:ncol(pp)] <- paste0("count_", sample_names)
            cat("Renamed sample columns to:", paste(colnames(pp)[4:ncol(pp)], collapse=", "), "\n")
        }

        # ← AGGIUNGI SYMBOL QUI
        pp_gr <- GenomicRanges::makeGRangesFromDataFrame(pp, keep.extra.columns = TRUE)
        overlaps <- GenomicRanges::findOverlaps(pp_gr, for_grange, type = "any")
        pp_gr_matched <- pp_gr[queryHits(overlaps)]
        mcols(pp_gr_matched)$SYMBOL <- for_grange[subjectHits(overlaps)]$SYMBOL
        pp <- as.data.frame(pp_gr_matched) %>%
            dplyr::select(-width, -strand) %>%
            dplyr::distinct()
        cat("Added SYMBOL to", nrow(pp), "positions\n")
        cat("BAM pileup: per-base coverage (1-based positions)\n\n")

        
    } else if (input$type_coverage == "bed") {
        cat("Processing BED coverage files (format: chr start end depth)...\n")
        cat("Each interval represents a uniform coverage block\n\n")
        
        # Read and process all coverage BED files
        df_list <- lapply(seq_along(list_coverage()), function(idx) {
            file <- list_coverage()[idx]
            cat(paste("Sample", idx, ":", basename(file), "\n"))
            
            # ✅ VALIDAZIONE QUI DENTRO (dove 'file' esiste)
            validate(
                need(file.exists(file),
                    paste("⚠️ File not found:", file, 
                        "\nVerify that the 'BAM/BED list' file contains full paths, not just filenames."))
            )

            # Check file is not binary (BAM files start with "BAM\1")
            first_bytes <- tryCatch({
                con <- file(file, "rb")
                on.exit(close(con))
                readBin(con, "raw", n = 4)
            }, error = function(e) raw(0))

            is_bam <- length(first_bytes) >= 4 && 
                    rawToChar(first_bytes[1:3], multiple = FALSE) == "BAM"

            if (is_bam) {
                waiter::waiter_hide()
                showNotification(
                    paste0("❌ File '", basename(file), "' appears to be a BAM file, not a BED coverage file.",
                        "\nPlease select 'BAM file' as the file format, or provide BED coverage files."),
                    type = "error",
                    duration = 15
                )
                return(NULL)
            }
            
            df <- tryCatch({
                result <- read.table(
                    file,
                    header = FALSE,
                    stringsAsFactors = FALSE,
                    colClasses = c("character", "integer", "integer", "numeric")
                )
                if (ncol(result) != 4) {
                    stop(paste("Expected 4 columns (chr, start, end, depth) but got", ncol(result)))
                }
                colnames(result) <- c("seqnames", "start", "end", "count")
                result
            }, error = function(e) {
                waiter::waiter_hide()
                showNotification(
                    paste0("❌ Could not read file '", basename(file), "' as BED coverage format.",
                        "\nExpected format: chr start end depth (4 columns, tab-separated).",
                        "\nMake sure you selected the correct file type (BAM vs BED).",
                        "\nError: ", e$message),
                    type = "error",
                    duration = 15
                )
                return(NULL)
            })

            if (is.null(df)) return(NULL)
            
            cat(paste("  Raw intervals:", nrow(df), "\n"))
            
            # Convert to 1-based if input is 0-based
            if (input_coord_system == "0-based") {
                df <- df %>% dplyr::mutate(start = start + 1)
                cat("  Converted from 0-based to 1-based\n")
            }
            
            # CRITICAL FIX: Add unique identifier to preserve count through intersect
            df$interval_id <- paste0("interval_", seq_len(nrow(df)))
            
            # Convert to GRanges (now definitely 1-based, with metadata preserved)
            df_gr <- GenomicRanges::makeGRangesFromDataFrame(
                df, 
                keep.extra.columns = TRUE  # Keeps count AND interval_id
            )
            
            # Find overlaps with target regions
            overlaps <- GenomicRanges::findOverlaps(df_gr, for_grange, type = "any")
            
            if (length(overlaps) == 0) {
                warning(paste("No overlaps found for", basename(file)))
                return(data.frame(
                    seqnames = character(), 
                    start = integer(), 
                    end = integer(), 
                    count = numeric(),
                    SYMBOL = character(), 
                    stringsAsFactors = FALSE
                ))
            }
            
            cat(paste("  Overlapping intervals:", length(unique(queryHits(overlaps))), "\n"))
            
            # CRITICAL FIX: Use pintersect to trim AND preserve metadata
            df_gr_overlapping <- df_gr[queryHits(overlaps)]
            for_grange_overlapping <- for_grange[subjectHits(overlaps)]
            
            # pintersect automatically preserves all metadata (count, interval_id)
            df_gr_trimmed <- GenomicRanges::pintersect(
                df_gr_overlapping,
                for_grange_overlapping,
                ignore.strand = TRUE
            )
            
            cat(paste("  Trimmed to", length(df_gr_trimmed), "intervals\n"))
            
            # Convert back to dataframe - count is automatically preserved!
            df_trimmed <- as.data.frame(df_gr_trimmed)
            df_trimmed$SYMBOL <- for_grange_overlapping$SYMBOL
            
            # Clean up: remove metadata columns we don't need
            df_trimmed <- df_trimmed %>%
                dplyr::select(seqnames, start, end, count, SYMBOL) %>%
                dplyr::distinct()  # Remove any duplicate intervals created by overlapping targets
            
            cat(paste("  Final intervals:", nrow(df_trimmed), "\n\n"))
            
            return(df_trimmed)
        })  # ✅ CHIUDI IL LAPPLY QUI!
        
        if (any(sapply(df_list, is.null))) return(NULL)
        validate(
            need(any(sapply(df_list, nrow) > 0),
                "⚠️ No coverage data overlaps with target regions!")
        )

        # Common blocks for all samples: cut at every block boundary of every sample
        # and at every target boundary, so samples are compared block by block.
        # Target bases missing from a coverage file get 0.
        sample_grs <- lapply(df_list, function(d) {
            GenomicRanges::GRanges(as.character(d$seqnames),
                                   IRanges::IRanges(d$start, d$end), count = d$count)
        })
        pieces <- GenomicRanges::disjoin(c(
            GenomicRanges::granges(for_grange),
            unlist(GenomicRanges::GRangesList(lapply(sample_grs, GenomicRanges::granges)))))
        pp <- data.frame(seqnames = as.character(GenomicRanges::seqnames(pieces)),
                         start = GenomicRanges::start(pieces),
                         end = GenomicRanges::end(pieces),
                         stringsAsFactors = FALSE)
        for (k in seq_along(sample_grs)) {
            v <- numeric(length(pieces))
            h <- GenomicRanges::findOverlaps(pieces, sample_grs[[k]])
            v[S4Vectors::queryHits(h)] <- sample_grs[[k]]$count[S4Vectors::subjectHits(h)]
            pp[[paste0("count", k)]] <- v
        }

        # Add SYMBOL (a block shared by two genes appears once per gene)
        h <- GenomicRanges::findOverlaps(pieces, for_grange)
        pp <- pp[S4Vectors::queryHits(h), , drop = FALSE]
        pp$SYMBOL <- for_grange$SYMBOL[S4Vectors::subjectHits(h)]
        rownames(pp) <- NULL
      
        # Replace NA with 0 for missing coverage
        # Fix column names: use actual sample names with count_ prefix
        cat(paste("Merged result:", nrow(pp), "unique intervals\n"))
        cat("NA count after merge:", sum(is.na(pp)), "\n")
        cat("Zero count after merge:", sum(pp == 0, na.rm = TRUE), "\n")

        # Fix column names: use actual sample names with count_ prefix
        sample_names <- tools::file_path_sans_ext(basename(list_coverage()))

        # Check for duplicates and fix
        duplicates <- sample_names[duplicated(sample_names)]
        if (length(duplicates) > 0) {
            sample_names <- make.unique(sample_names, sep = "_")
            showNotification(
                paste("⚠️ Duplicate sample names detected:", 
                    paste(unique(duplicates), collapse=", "),
                    "- automatically renamed to avoid conflicts."),
                type = "warning",
                duration = 10
            )
            cat("WARNING: Duplicate sample names renamed to:", paste(sample_names, collapse=", "), "\n")
        }

        count_cols <- colnames(pp)[!colnames(pp) %in% c("seqnames", "start", "end", "SYMBOL")]

        if (length(sample_names) == length(count_cols)) {
            colnames(pp)[colnames(pp) %in% count_cols] <- paste0("count_", sample_names)
            cat("Renamed sample columns to:", paste(paste0("count_", sample_names), collapse=", "), "\n")
        }
        cat(paste("Merged result:", nrow(pp), "unique intervals\n"))
    }
    cat("\n=== COVERAGE PROCESSING COMPLETE ===\n")
    # Normalize chromosome names
    pp$seqnames <- paste0("chr", sub("^chr", "", as.character(pp$seqnames)))
    
    cat(paste("\nFinal coverage dimensions:", nrow(pp), "intervals x", ncol(pp)-4, "samples\n"))
   
    # Verify SYMBOL column
    validate(
        need("SYMBOL" %in% colnames(pp), 
            " Internal error: incomplete coverage data. Please try again.")
    )

    # Check for missing genes
    # Every gene is in the table now (0x if no reads): list genes with no reads at all
    count_cols_check <- setdiff(colnames(pp), c("seqnames", "start", "end", "SYMBOL"))
    has_reads <- tapply(rowSums(pp[, count_cols_check, drop = FALSE]) > 0, pp$SYMBOL, any)
    missing_genes <- names(has_reads)[!has_reads]

    if (length(missing_genes) > 0) {
        cat(paste("\nWARNING:", length(missing_genes), "genes have NO coverage!\n"))
        cat("First 20 missing genes:", paste(head(missing_genes, 20), collapse=", "), "\n")
    }

    # Store coordinate system metadata
    attr(pp, "coordinate_system") <- "1-based"
    attr(pp, "trimmed_to_targets") <- TRUE
    attr(pp, "original_input_system") <- input_coord_system
    
    cat("\nSUCCESS: Coverage processing completed successfully!\n")
    cat(paste("  Total intervals:", nrow(pp), "\n"))
    cat(paste("  Total genes:", length(unique(pp$SYMBOL)), "\n\n"))
    
    return(pp)
})

name_sample <- reactive({
   if (is.null(list_coverage())) return(character(0))
   tools::file_path_sans_ext(basename(list_coverage()))
})

# ============================================================================
# COMPLETELY FIXED STATISTICS with proper weighted calculations
# ============================================================================

stat_summ <- reactive({
    start_time <- Sys.time()
    if (is.null(gene_list())) return(NULL)
    if (is.null(list_coverage())) return(NULL)
    if (!is.null(no_entrID()) && nrow(no_entrID()) != 0) return(no_entrID())
    
    # Assicurati che coverage_input sia stato calcolato
    cov_data <- coverage_input()

    validate(
        need(!is.null(cov_data), 
            "No data of coverage available. Go to 'Coverage Analysis', upload BAM/BED files and click 'Process Coverage'."),
        need(nrow(cov_data) > 0,
            "Coverage calculated but no data found.")
    )
    cat("\n=== STATISTICS CALCULATION ===\n")
    ppinp <- as.data.frame(cov_data)
    # req(coverage_input())

    # cat("\n=== STATISTICS CALCULATION ===\n")
    # ppinp <- as.data.frame(coverage_input())
    validate(
        need("SYMBOL" %in% colnames(ppinp), 
        "Internal error: incomplete coverage data. Please try again.")
    )

    # RINOMINA
    colnames(ppinp)[colnames(ppinp) == "seqnames"] <- "chromosome"
    cat("\nDopo il rename:\n")
    cat("  Colonne in ppinp:", paste(colnames(ppinp), collapse=", "), "\n")

    # Identify sample columns (already have count_ prefix from coverage_input)
    coord_cols <- c("chromosome", "start", "end", "SYMBOL")
    sample_cols_with_prefix <- setdiff(colnames(ppinp), coord_cols)  # ["count_file1", "count_file2"]

    # Remove count_ prefix for cleaner names
    sample_cols <- sub("^count_", "", sample_cols_with_prefix)  # ["file1", "file2"]

    # Rename columns in ppinp
    colnames(ppinp)[colnames(ppinp) %in% sample_cols_with_prefix] <- sample_cols

    n <- length(sample_cols)

    cat("Sample columns found:", paste(sample_cols, collapse=", "), "\n")
    cat("Number of samples:", n, "\n\n")

    cat(paste("Input data:\n"))
    cat(paste("  Total intervals:", nrow(ppinp), "\n"))
    cat(paste("  Unique genes:", length(unique(ppinp$SYMBOL)), "\n"))
    cat(paste("  Samples:", n, "\n"))
    cat(paste("  Sample names:", paste(sample_cols, collapse=", "), "\n\n"))

    # Interval width for weighted stats
    ppinp$width <- ppinp$end - ppinp$start + 1
    ppinp[sample_cols] <- lapply(ppinp[sample_cols], as.numeric)

    # Prepare per-sample dataframes
    x <- lapply(seq_along(sample_cols), function(j) {
        i <- sample_cols[j]
        if (!(i %in% colnames(ppinp))) {
            cat("ERRORE: Colonna", i, "non trovata!\n")
            cat("Colonne disponibili:", paste(colnames(ppinp), collapse=", "), "\n")
            return(NULL)
        }
        x_df <- ppinp %>%
            dplyr::select("chromosome", "start", "end", "SYMBOL", "width", tidyselect::all_of(i))
        colnames(x_df)[6] <- "value"
        x_df$sample <- i
        x_df
    })

    # Number of contiguous regions below threshold (adjacent low bases/blocks = 1 region)
    count_low_regions <- function(chr, start, end, low) {
        low[is.na(low)] <- FALSE
        if (!any(low)) return(0L)
        o <- order(chr, start)
        chr <- chr[o][low[o]]; start <- start[o][low[o]]; end <- end[o][low[o]]
        n <- length(start)
        sum(c(TRUE, chr[-1] != chr[-n] | start[-1] != end[-n] + 1L))
    }

    # Statistical calculation (weighted mean, correct percentage)
    statistical_operation <- function(df) {
        df <- df %>% dplyr::filter(!is.na(value), !is.na(width), width > 0)
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
                Median_coverage = {
                                    df_sub <- data.frame(value = value, width = width)
                                    df_sub <- df_sub[!is.na(df_sub$value) & !is.na(df_sub$width) & df_sub$width > 0, ]
                                    if (nrow(df_sub) == 0) NA_real_ else {
                                        ord <- order(df_sub$value)
                                        df_sub <- df_sub[ord, ]
                                        cumw <- cumsum(df_sub$width)
                                        df_sub$value[which(cumw >= sum(df_sub$width) / 2)[1]]
                                    }
                                },
                number_of_intervals_under_20x = count_low_regions(chromosome, start, end, value < 20),
                bases_under_20x = sum(width[value < 20], na.rm = TRUE),
                percentage_bases_under_20x = (sum(width[value < 20], na.rm = TRUE) / 
                                              sum(width, na.rm = TRUE)) * 100,
                .groups = "drop"
            ) %>%
            as.data.frame() %>%
            dplyr::mutate_if(is.numeric, round, 3)
    }

    out_r <- do.call(rbind, lapply(x, statistical_operation))

    cat("\n=== RESULTS ===\n")
    cat(paste("Total rows:", nrow(out_r), "\n"))
    cat(paste("Unique genes in output:", length(unique(out_r$SYMBOL)), "\n"))
    cat(paste("\nStatistics computation took:",
              round(difftime(Sys.time(), start_time, units = "mins"), 2),
              "minutes\n\n"))

    # Gene list diagnostic (optional for future debugging)
    if (input$type_file == "target_bed") {
        input_genes <- unique(gene_list()$SYMBOL)
    } else {
        validated <- validated_genes_data()
        input_genes <- unique(validated$my_gene_name$ALIAS)
    }
    output_genes <- unique(out_r$SYMBOL)
    missing_genes <- setdiff(input_genes, output_genes)
    if (length(missing_genes) > 0) {
        cat(paste("\nWARNING:", length(missing_genes), "genes have NO coverage!\n"))
        cat("First 20 missing genes:", paste(head(missing_genes, 20), collapse=", "), "\n")
    } else {
        cat("\nSUCCESS: All validated genes have coverage data!\n")
    }
    return(out_r)
})

# Overlapping mates of a read pair (same DNA fragment read from both ends):
# Rsamtools::pileup() counts the overlapping bases twice. This function returns,
# for each base, how many pairs were counted twice: both mates cover the base and
# both pass the pileup filters (MAPQ, base quality, base not N). Subtracting it
# from the pileup count gives "count once" (as samtools depth -s; as mosdepth with
# proper_only = TRUE). Deletions count as covered, as in the pileup.
# Target regions are read in blocks, to keep memory low with large BAM files.
mate_overlap_counts <- function(bam, param, min_mapq, min_bq, proper_only = FALSE,
                                regions_per_block = 500) {
  empty <- data.frame(seqnames = character(), pos = integer(), n_overlap = integer(),
                      stringsAsFactors = FALSE)
  flag <- Rsamtools::bamFlagAND(
    Rsamtools::bamFlag(param, asInteger = TRUE),
    Rsamtools::scanBamFlag(isPaired = TRUE, isUnmappedQuery = FALSE, hasUnmappedMate = FALSE,
                           isSupplementaryAlignment = FALSE,
                           isProperPair = if (proper_only) TRUE else NA))
  wh <- Rsamtools::bamWhich(param)
  regions <- GenomicRanges::sort(GenomicRanges::GRanges(rep(names(wh), lengths(wh)),
                                                        unlist(wh, use.names = FALSE)))
  blocks <- split(regions, ceiling(seq_along(regions) / regions_per_block))
  # pairs already handled in a previous block (kept only while they can reach the next block)
  done <- data.frame(q = character(0), chr = character(0), end = integer(0), stringsAsFactors = FALSE)
  res <- vector("list", length(blocks))
  for (bl in seq_along(blocks)) {
    b1 <- blocks[[bl]][1]
    done <- done[done$chr == as.character(GenomicRanges::seqnames(b1)) &
                   done$end >= GenomicRanges::start(b1), , drop = FALSE]
    sb <- Rsamtools::scanBam(bam, param = Rsamtools::ScanBamParam(
      which = blocks[[bl]], flag = flag,
      what = c("qname", "flag", "rname", "pos", "cigar", "mapq", "seq", "qual")))
    col <- function(x) do.call(c, unname(lapply(sb, `[[`, x)))
    q <- col("qname")
    if (length(q) == 0) next
    fl <- col("flag"); chr <- as.character(col("rname")); st <- col("pos")
    cig <- col("cigar"); mq <- col("mapq")
    # pair the two mates (a read in two target regions is returned twice: keep it once)
    qid <- match(q, unique(q)); m2 <- bitwAnd(fl, 128L) > 0
    ok <- !is.na(mq) & mq >= min_mapq & !duplicated(qid * 2 + m2) & !(q %in% done$q)
    i <- which(ok); i <- i[qid[i] %in% qid[i][duplicated(qid[i])]]
    if (length(i) == 0) next
    i <- i[order(qid[i], m2[i])]
    a <- i[seq(1, length(i), 2)]; b <- i[seq(2, length(i), 2)]
    en <- st + GenomicAlignments::cigarWidthAlongReferenceSpace(cig) - 1L
    done <- rbind(done, data.frame(q = q[a], chr = chr[a], end = pmax(en[a], en[b]),
                                   stringsAsFactors = FALSE))
    s <- pmax(st[a], st[b]); e <- pmin(en[a], en[b])
    ov <- qid[a] == qid[b] & chr[a] == chr[b] & !is.na(s) & !is.na(e) & s <= e
    if (!any(ov)) next
    a <- a[ov]; b <- b[ov]; s <- s[ov]; e <- e[ov]; w <- e - s + 1L
    qual <- col("qual"); sq <- col("seq")
    # does each mate pass the base filters at each overlapping base?
    pass <- function(m) {
      off <- s - st[m] + 1L
      ql <- GenomicAlignments::sequenceLayer(qual[m], cig[m], from = "query", to = "reference",
                                             D.letter = "~", N.letter = "!")
      ok_b <- as.integer(unlist(Biostrings::subseq(ql, off, width = w))) - 33L >= min_bq
      has_n <- which(Biostrings::vcountPattern("N", sq[m]) > 0)   # pileup ignores N bases
      if (length(has_n) > 0) {
        sl <- GenomicAlignments::sequenceLayer(Biostrings::BStringSet(sq[m][has_n]), cig[m][has_n],
                                               from = "query", to = "reference",
                                               D.letter = "A", N.letter = "A")
        is_n <- as.integer(unlist(Biostrings::subseq(sl, off[has_n], width = w[has_n]))) == 78L
        first <- cumsum(c(0L, w))[has_n]
        idx <- unlist(lapply(seq_along(has_n), function(k) first[k] + seq_len(w[has_n[k]])))
        ok_b[idx[is_n]] <- FALSE
      }
      ok_b
    }
    both <- pass(a) & pass(b)
    if (!any(both)) next
    res[[bl]] <- data.frame(seqnames = rep(chr[a], w)[both],
                            pos = (rep(s, w) + sequence(w) - 1L)[both], stringsAsFactors = FALSE)
  }
  res <- do.call(rbind, res)
  if (is.null(res) || nrow(res) == 0) return(empty)
  out <- lapply(split(res$pos, res$seqnames), function(p) {
    p0 <- min(p); n <- tabulate(p - p0 + 1L); k <- which(n > 0)
    data.frame(pos = p0 + k - 1L, n_overlap = n[k])
  })
  do.call(rbind, Map(function(x, nm) cbind(seqnames = nm, x, stringsAsFactors = FALSE), out, names(out)))
}