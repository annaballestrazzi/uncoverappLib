#' @title unCOVERApp user interface
#' @description This function allows you to open graphical interface.
#' @param agree visualization.
#' @keywords app
#' @export
#' @examples
#' uncoverappLib.run()

suppressPackageStartupMessages({
require(shiny)
require(shinyWidgets)
require(shinyBS)
require(shinyjs)
require(markdown)
require(waiter) 
options("browser" = "xdg-open")
})

first.file <- system.file(
  "extdata",
  "intro.md",
  package = "uncoverappLib"
)
second.file <- system.file(
  "extdata",
  "script.js",
  package = "uncoverappLib"
)

third.file <- system.file(
  "extdata",
  "CONTACTS.md",
  package = "uncoverappLib"
)

intro_processing_short <- system.file(
  "extdata",
  "prep_input_short.md",
  package = "uncoverappLib"
)

intro_processing_full <- system.file(
  "extdata",
  "prep_input_full.md",
  package = "uncoverappLib"
)

mdStyle <- "margin-left: 30px; margin-right: 30px"


intro <- function() {
  tabPanel("Home",
           shiny::includeMarkdown(first.file),
           #includeMarkdown(file.path(".","intro.md")),
           style=mdStyle
  )
}


preprocess <- function() {
  shiny::tabPanel("Processing and Statistical Summary",
                   # Short intro (always visible)
                   shiny::includeMarkdown(intro_processing_short),
                   
                   # Collapsible full guide
                   shinyBS::bsCollapse(
                     id = "prep_guide_collapse",
                     open = NULL,  # Closed by default
                     shinyBS::bsCollapsePanel(
                       title = HTML("<b>📖 Show detailed guide</b>"),
                       shiny::includeMarkdown(intro_processing_full),
                       style = "info"
                     )
                   ),
                   
                   h1(strong("Prepare your input file")),

                       fluidPage(
                         sidebarLayout(
                      sidebarPanel(
                        shiny::selectInput("Genome",
                                           label = "Reference Genome",
                                           choices = c("hg19",
                                                       "hg38"),
                                           selected = "UCSC genome"),
                        hr(),
                        hr(),
                        shiny::selectInput("type_coverage",
                                           label = "File format",
                                           choices = c("BAM file" = "bam",
                                                       "BED coverage file"= "bed"),
                                           selected = "bed"),
                        hr(),
                        # NUOVO CODICE (con conditional panels):

                        # Parametri BAM (visibili solo se type_coverage == "bam")
                        conditionalPanel(
                          condition = "input.type_coverage == 'bam'",
                          hr(),
                          h4(strong("BAM File Parameters")),
                          shinyWidgets::pickerInput("MAPQ",
                            label = "Minimum Mapping Quality (MAPQ)",
                            choices = c(1:1000),
                            options = list(`actions-box` = TRUE), 
                            multiple = FALSE,
                            selected = 20),
  
                          shinyWidgets::pickerInput("base_qual",
                            label = "Minimum Base Quality",
                            choices = c(1:1000),
                            options = list(`actions-box` = TRUE), 
                            multiple = FALSE,
                            selected = 20),
                          helpText(em("These parameters apply only to BAM pileup processing"))
                        ),

                        # Parametri BED (visibili solo se type_coverage == "bed")
                        conditionalPanel(
                          condition = "input.type_coverage == 'bed'",
                          hr(),
                          h4(strong("BED Coverage File Parameters")),
                          radioButtons(
                            "input_coordinate_system",
                            "Coordinate system:",
                            choices = c("0-based (BED standard)" = "0-based",
                                        "1-based (IGV, VCF)" = "1-based"),
                            selected = "0-based"
                          ),
                          helpText(em("BED files are typically 0-based (half-open intervals)"))
                        ),

                        hr(),

                        shinyWidgets::pickerInput("notation",
                                                  label = "Chromosome Notation",
                                                  choices = c("chr", "number"),
                                                  options = list(`actions-box` = TRUE),
                                                  multiple =FALSE),



                        # sliderInput("min_coverage",
                        #             label = "Minimum Coverage Depth (BED only)",
                        #             value = 20,
                        #             min = 1,
                        #             max = 1000,
                        #             step = 1)
                        hr(),

                        fileInput(inputId = "gene1",
                                  label = "Load input filtering file",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".csv")),
                        helpText("Choose the type of your input file"),
                        radioButtons("type_file", "Type of file",
                                     choices = c("List of genes name" = "gene_name",
                                                 "Target Bed "= "target_bed"),
                                     selected = "gene_name"),
                        br(),
                        br(),
                        helpText("The file must contain one path per row"),
                        fileInput(inputId = "bam_list",
                                  label = "Load your genomic data(BAM list or BED coverage list)",
                                  accept = c("text/csv",
                                             ".zip",
                                             ".gz",
                                             ".bed",
                                             "text/comma-separated-values,text/plain",
                                             ".list")),

                        hr(),
                        hr(),
                        actionButton("process_coverage", 
                                     "Process Coverage Files",
                                     icon = icon("play-circle"),
                                     style = "color: #fff; background-color: #3498db; 
                                              border-color: #2980b9; width: 100%; 
                                              margin-bottom: 10px; font-weight: bold;"),
                        
                        helpText(em("Click to generate the input table from your files")),
                        
                        hr(),
                        downloadButton("summary", "Download Statistical Summary", class = "btn-primary",
                                        style='color: #fff !important; background-color: #27ae60;
                                        border-color: #fff;
                                        padding: 15px 14px 15px 14px;
                                        margin: 15px 5px 5px 5px;
                                        font-size: 14px !important;
                                        font-weight: bold !important;')
                      ),
                      shiny::mainPanel(
                        tabsetPanel(
                          tabPanel(title= "input for uncoverapp",
                                   shinycssloaders::withSpinner(
                                     DT::DTOutput("input1"))
                          )
                        ) )
                    )))
}

myHome <- function() {
  tabPanel("Coverage Analysis",
           h1(strong("Coverage Analysis")),
           #Interactive web-application to visualize and annotate low-coverage positions in clinical sequencing
           helpText(
                    em("Note: Select input options"),
                    br(),
                    ("Upload your input bed file with columns: "),
                    br(),
                    ("chromosome, start, end, coverage by sample")
                       ),
           shinyjs::useShinyjs(),
           includeScript(second.file),
           
           sidebarPanel(
             # ═══════════════════════════════════════════════════════════
             # STANDARD FILTERS (always visible)
             # ═══════════════════════════════════════════════════════════
             h4(strong("Standard Filters")),
             selectInput("UCSC_Genome", label = "Reference Genome",
                         choices = c("hg19", "hg38"), selected = "hg19"),
             hr(),
             
             shinyWidgets::pickerInput("coverage_co", label = "Coverage threshold",
                                       choices = c(1:1000, "all"),
                                       options = list(`actions-box` = TRUE), multiple = FALSE),
             helpText(em("Select minimum value as coverage threshold")),
             hr(),
             
             textInput(inputId = "Sample", label = "Sample"),
             helpText(em("Example: example_POLG.bam")),
             hr(),
             
             # ═══════════════════════════════════════════════════════════
             # FILTER BY (radio buttons)
             # ═══════════════════════════════════════════════════════════
             h4(strong("Filter By")),
             radioButtons("filter_by", NULL,
                          choices = c("Gene name" = "gene",
                                      "Chromosome" = "chromosome", 
                                      "All chromosomes" = "all_chr",
                                      "Region coordinates" = "region"),
                          selected = "gene"),
             
             # Conditional inputs based on filter_by selection
             conditionalPanel(
               condition = "input.filter_by == 'gene'",
               textInput(inputId = "Gene_name", label = "Gene name"),
               helpText(em("Write gene name (e.g., POLG)")),
               actionButton("ucsc_lookup", "Lookup UCSC Gene", 
                            icon = icon("search"),
                            style = "color: #fff; background-color: #337ab7; border-color: #2e6da4; margin-bottom: 10px;")
             ),
             
             conditionalPanel(
               condition = "input.filter_by == 'chromosome'",
               selectInput("Chromosome", "Select Chromosome:",
                           choices = c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", 
                                       "chr7", "chr8", "chr9", "chr10", "chr11", "chr12",
                                       "chr13", "chr14", "chr15", "chr16", "chr17", "chr18",
                                       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY", "chrM"),
                           selected = "chr1")
             ),
             
             conditionalPanel(
               condition = "input.filter_by == 'all_chr'",
               helpText(em("Search will be performed across all chromosomes"), 
                        style = "color: #337ab7; font-weight: bold;")
             ),
             
             conditionalPanel(
               condition = "input.filter_by == 'region'",
               textInput(inputId = "query_Database", label = "Region coordinates"),
               helpText(em("Example: chr7:37453745-45627643 or 7:37453745-45627643"))
             ),
             hr(),
             
             # ═══════════════════════════════════════════════════════════
             # ACTION BUTTONS (trigger calculations ONLY when clicked)
             # ═══════════════════════════════════════════════════════════
             h4(strong("Calculate")),
             actionButton("calc_low_coverage", 
                          "Calculate Low Coverage Regions",
                          icon = icon("calculator"),
                          style = "color: #fff; background-color: #5cb85c; border-color: #4cae4c; width: 100%; margin-bottom: 10px;"),
             
             actionButton("calc_annotations", 
                          "Calculate Annotations on Low Coverage",
                          icon = icon("dna"),
                          style = "color: #fff; background-color: #f0ad4e; border-color: #eea236; width: 100%; margin-bottom: 10px;"),
             
             shinyBS::bsButton("remove", label = "Refresh", icon = icon("refresh"),
                               style = "default", size = "small"),
             
             shinyjs::hidden(p(id = "text1", "Running.....")),
             hr(),
             
             # Download button (visible only after annotation calculation)
             conditionalPanel(
               condition = "output.annotation_ready == true",
               downloadButton("downloadData", "Download Annotations", class = "btn-primary",
               style = 'padding: 8px !important; 
                        font-size: 14px !important; 
                        width: 100%;
                        color: #fff !important;
                        font-weight: bold !important;')
             ),
             hr(),
             
             # Close button
             shiny::tags$button(
               id = 'close', type = "button", class = "btn action-button",
               style = 'color: white; background-color: #dd4b39; padding:4px; font-size:120%',
               onclick = "setTimeout(function(){window.close();},500);",
               "Close App"
             )
           ),
           
           mainPanel(
             fileInput(inputId = "file1",
                       label = "Select input file",
                       accept = c("text/csv", ".bedGraph", ".bedGraph.gz", ".zip", 
                                  ".gz", ".bed", "text/comma-separated-values,text/plain", ".csv")),
             checkboxInput("header", "Header", TRUE),
             shinyBS::bsButton("pileup", 
                               label= "load input file", 
                               icon = icon("file"), style = "info",
                               size = "default"),
             hr(),
             
             tabsetPanel(
               tabPanel("bed file", shinycssloaders::withSpinner(DT::DTOutput("text"))),
               tabPanel("Low-coverage positions", DT::DTOutput("text_cv")),
               tabPanel("UCSC gene", shinycssloaders::withSpinner(DT::DTOutput('ccg'))),
               tabPanel("Gene coverage",
                        actionButton("generate_gene_plot", 
                                      "Generate Gene Coverage Plot",
                                      icon = icon("chart-area"),
                                      style = "color: #fff; background-color: #3498db; 
                                              border-color: #2980b9; width: 100%; 
                                              margin-bottom: 15px; font-weight: bold;"),
                        shinycssloaders::withSpinner(plotOutput("all_gene")),
                        DT::DTOutput('df.l')
                ),
               tabPanel("Annotations on low-coverage positions",
                        helpText(em("dbSNP-annotation collects all consequences found in VEP-defined canonical transcripts")),
                        shinycssloaders::withSpinner(DT::DTOutput("uncover_position"))),
               id = "tabSet"
             ),
             hr(),
             
             fluidRow(column(12, DT::DTOutput("x4")))
           )
  )
}


myTab1 <- function() {
  tabPanel("Calculate AF by allele frequency app",

           # Application title
           titlePanel("Maximum credible population allele frequency"),

           ##### Bootstrap method for page costruction
           fluidPage(
             fluidRow(
               ##### Sidebar
               column(8,wellPanel(radioButtons("inh",
                                               "Inheritance:",
                                               choices = list("monoallelic",
                                                              "biallelic"),
                                               selected = "monoallelic"),

                                  numericInput("prev","Prevalence = 1 in ...
                                               (people)",
                                               min = 1,max = 1e8,value = 500),
                                  options = NULL),
                      br(),
                      sliderInput("hetA","Allelic heterogeneity:",
                                  min = 0, max = 1,value = 0.1),
                      sliderInput("hetG",
                                  "Genetic heterogeneity:",
                                  min = 0, max = 1,value = 1),
                      br(),
                      sliderInput("pen", "Penetrance:",
                                  min = 0, max = 1, value = 0.5))),
             br(),
             column(8,
                    h3("Maximum credible population AF:"),
                    h2(textOutput("maxAF"),align="center",style = "color:blue")),
             column(8,
                    h3("Uncover position",
                       helpText(em("Low-coverage positions excluding sites
                                   annotated as variants with AF> maxAF
                                   (default maxAF value: 5%)"),align="center",
                                style="color:blue"),
                       style = "font-size: 100%; width: 100%",
                       shinycssloaders::withSpinner(
                       DT::DTOutput("uncoverPosition")))),
             br(),
             br(),
             downloadButton("download_maxAF", "Download_maxAF",
                            class = "btn-primary",
                            style='padding:4px; font-size:80%',
                            helpText("download low coverage
                                     position dbSNFP-annotation filtered by
                                     maximum allele frequency",
                                     class = "btn-primary",
                                     style='padding:4px; font-size:60%'))
             #)
           ))}


myTab2 <- function() {
  tabPanel("Binomial distribution",
           titlePanel("Binomial distribution "),
           fluidRow(
             column(4,(numericInput("p",
                                    "Allele Fraction",
                                    min = 0,
                                    max = 1,
                                    value = 0.05)),
                    helpText(em("the expected fraction of variant reads
                    (probability of success)",
                    align="center",
                    style="color:gray")),
                    hr(),
                    numericInput("num_all",
                                 "Variant reads",
                                 min=0,
                                 max=100000,
                                 value=10),

                    helpText(em("the minimum number of variant reads required
                                by the user to support variant evidence,
                                (number of successes)"),align="center",
                             style="color:gray"),
                    hr(),

                    textInput(inputId = "start_gp",
                              label = "Genomic position"),

                    #textInput(inputId = "end_gp",
                     #         label = "END genomic position"),
                    helpText(em("Specify a genomic position of interest gene
                                  in which calculating the
                                binomial distribution"))),


             column(4,
                    h2("consideration:"),
                    h3(htmlOutput("ci"))),
             column(4,
                    h2("Binomial Distribution", plotOutput("bd"))),
             column(10, h2("Cumulative distribution function"),
                    h3(plotOutput("pbinom"))))

  )}


myabout <- function() {
  tabPanel("Contacts",
           includeMarkdown(third.file),
           style=mdStyle
  )
}

ui <- shinyUI(
  tagList(
    shinyjs::useShinyjs(), 
    waiter::use_waiter(),
    shiny::tags$head(
      # Import font Inter (usato da Notion)
      shiny::tags$link(
        rel = "stylesheet",
        href = "https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700;800&display=swap"
      ),
      
      shiny::tags$style(HTML("
        /* ============================================ */
        /* RESET BASE - Applica Inter a TUTTO */
        /* ============================================ */
        * {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, 'Segoe UI', 'Roboto', sans-serif !important;
        }
        
        body, .container-fluid, .well, .tab-content, 
        .tab-pane, .markdown-body, p, div, span, li, td, th {
          font-family: 'Inter', -apple-system, BlinkMacSystemFont, sans-serif !important;
          font-size: 15px !important;
          line-height: 1.6 !important;
          color: #37352F !important;
          -webkit-font-smoothing: antialiased !important;
          -moz-osx-font-smoothing: grayscale !important;
        }
        
        /* ============================================ */
        /* TITOLI - Gerarchia chiara e distinta */
        /* ============================================ */
        
        h1 {
          font-size: 36px !important;
          font-weight: 800 !important;
          line-height: 1.2 !important;
          color: #111111 !important;
          margin-top: 40px !important;
          margin-bottom: 24px !important;
          padding-bottom: 16px !important;
          border-bottom: 2px solid #E0E0E0 !important;
          letter-spacing: -0.5px !important;
        }
        
        h2 {
          font-size: 28px !important;
          font-weight: 700 !important;
          line-height: 1.3 !important;
          color: #1a1a1a !important;
          margin-top: 48px !important;
          margin-bottom: 20px !important;
          padding-left: 0px !important;
          letter-spacing: -0.3px !important;
        }
        
        h3 {
          font-size: 20px !important;
          font-weight: 600 !important;
          line-height: 1.4 !important;
          color: #3498DB !important;
          margin-top: 32px !important;
          margin-bottom: 16px !important;
          letter-spacing: -0.2px !important;
        }
        
        h4 {
          font-size: 17px !important;
          font-weight: 600 !important;
          line-height: 1.4 !important;
          color: #555555 !important;
          margin-top: 24px !important;
          margin-bottom: 12px !important;
        }
        
        h5, h6 {
          font-size: 15px !important;
          font-weight: 600 !important;
          color: #666666 !important;
          margin-top: 16px !important;
          margin-bottom: 8px !important;
        }
        
        /* ============================================ */
        /* PARAGRAFI E TESTO */
        /* ============================================ */
        p {
          font-size: 15px !important;
          line-height: 1.7 !important;
          color: #37352F !important;
          margin-bottom: 16px !important;
        }
        
        ul, ol {
          font-size: 15px !important;
          line-height: 1.7 !important;
          color: #37352F !important;
          padding-left: 24px !important;
          margin-bottom: 16px !important;
        }
        
        li {
          margin-bottom: 8px !important;
        }
        
        /* ============================================ */
        /* CODE BLOCKS - Stile Notion */
        /* ============================================ */
        
        code {
          font-family: 'SF Mono', 'Monaco', 'Inconsolata', 'Courier New', monospace !important;
          font-size: 14px !important;
          background-color: rgba(135, 131, 120, 0.15) !important;
          color: #EB5757 !important;
          padding: 3px 6px !important;
          border-radius: 4px !important;
          font-weight: 500 !important;
        }
        
        pre {
          font-family: 'SF Mono', 'Monaco', 'Inconsolata', monospace !important;
          font-size: 13px !important;
          line-height: 1.5 !important;
          background-color: #F7F6F3 !important;
          color: #2C3E50 !important;
          padding: 16px !important;
          border-radius: 8px !important;
          border-left: 4px solid #3498DB !important;
          overflow-x: auto !important;
          margin-bottom: 20px !important;
        }
        
        pre code {
          background-color: transparent !important;
          color: #2C3E50 !important;
          padding: 0 !important;
          font-size: 13px !important;
        }
        
        /* ============================================ */
        /* LINKS - Stile Notion */
        /* ============================================ */
        a {
          color: #3498DB !important;
          text-decoration: none !important;
          border-bottom: 1px solid transparent !important;
          transition: all 0.2s ease !important;
          font-weight: 500 !important;
        }
        
        a:hover {
          border-bottom: 1px solid #3498DB !important;
          color: #2980B9 !important;
        }
        
        /* ============================================ */
        /* TABLES - Pulite e moderne */
        /* ============================================ */
        table {
          width: 100% !important;
          border-collapse: collapse !important;
          margin: 24px 0 !important;
          font-size: 14px !important;
          border-radius: 8px !important;
          overflow: hidden !important;
          box-shadow: 0 1px 3px rgba(0,0,0,0.1) !important;
        }
        
        th {
          background-color: #F7F6F3 !important;
          color: #37352F !important;
          font-weight: 600 !important;
          padding: 12px 16px !important;
          text-align: left !important;
          border-bottom: 2px solid #E0E0E0 !important;
        }
        
        td {
          padding: 12px 16px !important;
          border-bottom: 1px solid #F0F0F0 !important;
          color: #37352F !important;
        }
        
        tr:hover {
          background-color: #FAFAFA !important;
        }
        
        /* ============================================ */
        /* BLOCKQUOTES - Callout stile Notion */
        /* ============================================ */
        blockquote {
          border-left: 4px solid #3498DB !important;
          background-color: #F8F9FA !important;
          padding: 16px 20px !important;
          margin: 20px 0 !important;
          border-radius: 4px !important;
          font-size: 15px !important;
          color: #37352F !important;
        }
        
        blockquote p {
          margin: 0 !important;
        }
        
        /* ============================================ */
        /* STRONG/BOLD - Più visibile */
        /* ============================================ */
        strong, b {
          font-weight: 700 !important;
          color: #111111 !important;
        }
        
        /* ============================================ */
        /* EM/ITALIC */
        /* ============================================ */
        em, i {
          font-style: italic !important;
          color: #555555 !important;
        }
        
        /* ============================================ */
        /* HR - Separatori */
        /* ============================================ */
        hr {
          border: none !important;
          border-top: 1px solid #E0E0E0 !important;
          margin: 32px 0 !important;
        }
        
        /* ============================================ */
        /* LAYOUT - LARGHEZZA PIENA */
        /* ============================================ */
        
        /* Contenuto generale - larghezza piena */
        .container-fluid {
          width: 100% !important;
          padding-left: 15px !important;
          padding-right: 15px !important;
        }
        
        /* Tab panels - larghezza piena */
        .tab-content {
          width: 100% !important;
        }
        
        .tab-pane {
          width: 100% !important;
          padding: 20px !important;
        }
        
        /* SOLO il contenuto MARKDOWN centrato a 900px */
        .tab-pane > p,
        .tab-pane > h1,
        .tab-pane > h2,
        .tab-pane > h3,
        .tab-pane > h4,
        .tab-pane > ul,
        .tab-pane > ol,
        .tab-pane > blockquote,
        .tab-pane > pre {
          max-width: 900px !important;
          margin-left: auto !important;
          margin-right: auto !important;
        }
        
        /* Tabelle e grafici - larghezza piena */
        .tab-pane > .dataTables_wrapper,
        .tab-pane > table,
        .tab-pane > .shiny-plot-output,
        .tab-pane > .well,
        .tab-pane > .form-group {
          max-width: 100% !important;
          width: 100% !important;
        }
        
        /* ============================================ */
        /* NAVBAR - Bianca con bordo blu SOPRA */
        /* ============================================ */
        .navbar { 
          background: white !important;
          min-height: 120px !important;
          box-shadow: 0 2px 8px rgba(0, 0, 0, 0.1) !important;
          border-top: 4px solid #3498DB !important;
          border-bottom: none !important;
          margin-bottom: 0 !important;
          padding-top: 10px !important;
        }
        
        .navbar-brand { 
          padding: 10px 15px !important;
          height: auto !important;
        }
        
        .navbar-brand img { 
          height: 90px !important;
          width: auto !important;
          object-fit: contain !important;
          transition: transform 0.3s ease !important;
          display: block !important;
        }
        
        .navbar-brand img:hover {
          transform: scale(1.05) !important;
        }
        
        /* ============================================ */
        /* MENU ITEMS - Con border-radius fisso */
        /* ============================================ */
        
        .navbar-nav > li > a { 
          color: #2C3E50 !important;
          font-weight: 600 !important;
          font-size: 15px !important;
          margin-top: 35px !important;
          padding: 10px 20px !important;
          transition: all 0.3s ease !important;
          border-radius: 6px !important;
        }
        
        .navbar-nav > li > a:hover { 
          background-color: #ECF0F1 !important;
          color: #3498DB !important;
          transform: translateY(-2px) !important;
        }
        
        /* Menu attivo - mantiene border-radius */
        .navbar-nav > .active > a,
        .navbar-nav > .active > a:hover,
        .navbar-nav > .active > a:focus {
          background-color: #ECF0F1 !important;
          color: #3498DB !important;
          border-radius: 6px !important;
        }
        
        /* ============================================ */
        /* DROPDOWN MENU - Sfondo bianco */
        /* ============================================ */
        
        .dropdown-menu {
          background-color: white !important;
          border: 1px solid #ddd !important;
          box-shadow: 0 4px 12px rgba(0, 0, 0, 0.15) !important;
          border-radius: 6px !important;
          margin-top: 5px !important;
        }
        
        .dropdown-menu > li > a {
          color: #2C3E50 !important;
          padding: 10px 20px !important;
          transition: background-color 0.2s !important;
          font-weight: 500 !important;
        }
        
        .dropdown-menu > li > a:hover {
          background-color: #ECF0F1 !important;
          color: #3498DB !important;
        }
        
        /* ============================================ */
        /* BOTTONI - Moderni e consistenti */
        /* ============================================ */
        .btn {
          font-family: 'Inter', sans-serif !important;
          font-weight: 600 !important;
          font-size: 14px !important;
          border-radius: 6px !important;
          padding: 10px 20px !important;
          transition: all 0.2s ease !important;
        }
        
        .btn-primary {
          background-color: #3498DB !important;
          border-color: #3498DB !important;
        }
        
        .btn-primary:hover {
          background-color: #2980B9 !important;
          border-color: #2980B9 !important;
          transform: translateY(-1px) !important;
          box-shadow: 0 4px 8px rgba(52, 152, 219, 0.3) !important;
        }
        
        /* ============================================ */
        /* INPUT FIELDS */
        /* ============================================ */
        input, select, textarea {
          font-family: 'Inter', sans-serif !important;
          font-size: 14px !important;
          border-radius: 6px !important;
          border: 1px solid #D0D0D0 !important;
          padding: 8px 12px !important;
        }
        
        input:focus, select:focus, textarea:focus {
          border-color: #3498DB !important;
          box-shadow: 0 0 0 3px rgba(52, 152, 219, 0.1) !important;
          outline: none !important;
        }
      "))
    ),
    navbarPage(
      title = shiny::tags$a(
        href = "#",
        class = "navbar-brand",
        img(src = "logo.png")
      ),
      windowTitle = "uncoverApp",
      
      # Dropdown menu per alcune funzionalità principali
      navbarMenu("Menu",
                 intro(),
                 preprocess(),
                 myHome(),
                 myTab1(),
                 myTab2()
                 #myabout()
      ),
      
      # Tab separati
      #myTab1(),
      #myTab2(),
      myabout()
    )
  )
)