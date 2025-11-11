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


intro_processing <- system.file(
  "extdata",
  "prep_input.md",
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
                   shiny::includeMarkdown(intro_processing),
                       #shiny::includeMarkdown(file.path(".", "prep_input.md")),
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
                        helpText(em("Select filtering parameters for processing your input file")),
                        radioButtons(
                            "input_coordinate_system",
                            "Coverage file coordinate system:",
                            choices = c("0-based (BED standard)" = "0-based",
                                          "1-based (IGV, VCF)" = "1-based"),
                        selected = "0-based"
                        ),
                        hr(),

                        shinyWidgets::pickerInput("notation",
                                                  label = "Chromosome Notation",
                                                  choices = c("chr", "number"),
                                                  options = list(`actions-box` = TRUE),
                                                  multiple =FALSE),

                        shinyWidgets::pickerInput("MAPQ",
                                                  label = "Minum Mapping Quality (MAPQ)",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),

                        shinyWidgets::pickerInput("base_qual",
                                                  label = "Minimum Base Quality",
                                                  choices = c(1:1000),
                                                  options = list(`actions-box` = TRUE), multiple =FALSE),


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
                        helpText("Choose the type of your input file"),
                        radioButtons("type_coverage", "File format",
                                     choices = c("BAM file" = "bam",
                                                 "BED coverage file"= "bed"),
                                     selected = "bed"),
                        hr(),
                        hr(),
                        downloadButton("summary", "Statistical_Summary", class = "btn-primary",
                                       style='color: #fff; background-color: #27ae60;
                                       border-color: #fff;
                                       padding: 15px 14px 15px 14px;
                                       margin: 15px 5px 5px 5px; ')
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
           h1(strong("Interactive web-application to visualize and annotate low-coverage positions in clinical sequencing")),
           helpText(em("Note: Select input options",
                       span("Upload your input bed file with columns: chromosome, start, end, coverage by sample and nucleotide count", style = "color:blue"))),
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
                                       "chr19", "chr20", "chr21", "chr22", "chrX", "chrY"),
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
             
             # ═══════════════════════════════════════════════════════════
             # GENOMIC POSITION (for zooming plots)
             # ═══════════════════════════════════════════════════════════
             h4(strong("Genomic Position (for plots)")),
             textInput(inputId = "Start_genomicPosition", label = "START genomic position"),
             textInput(inputId = "end_genomicPosition", label = "END genomic position"),
             helpText(em("Change genomic interval for zooming")),
             hr(),
             
             # Download button (visible only after annotation calculation)
             conditionalPanel(
               condition = "output.annotation_ready == true",
               downloadButton("downloadData", "Download Annotations", class = "btn-primary",
                              style = 'padding:4px; font-size:120%; width: 100%;')
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
    shiny::tags$head(
      shiny::tags$style(HTML("
        .navbar { background-color: #FFFFFF;height: 150px;}
        .navbar-brand img { height: 150px;margin: 0px; }
        .navbar-nav > li > a { color: #333; font-weight: bold; font-size: 15px; margin-top: 60px }
        .navbar-nav > li > a:hover { color: #007BFF; }
        .navbar-brand { padding: 0px 10px; }
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

# Remove the extra closing brace that was causing the error
# The commented-out alternative UI definition can be removed or kept as a comment
#  ui <- shinyUI(navbarPage(
#    title = div(
#      img(src = "logo.png", 
#          height = "100px", 
#          style = "margin-top: 10px; margin-bottom: -10px;"),
#      style = "margin: 4; padding: 4;"
#    ),
#    windowTitle = "uncoverApp",  # This sets the browser tab title
#    intro(),
#    preprocess(),
#    myHome(),
#    myTab1(),
#    myTab2(),
#    myabout()
#  ))

