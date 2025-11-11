#' @title run.uncoverapp
#' @description
#' Launches \code{unCOVERApp}, a \code{Shiny} application for clinical assessment of sequence coverage.
#' The app can be launched in different environments using the \code{where} parameter:
#' \itemize{
#'   \item \code{"browser"} - launches in the user's default web browser,
#'   \item \code{"viewer"} - launches in the RStudio Viewer pane,
#'   \item \code{"window"} - launches in a new RStudio window.
#' }
#' 
#' @author Emanuela Iovino
#'
#' @examples
#' \dontrun{
#'   # Launch in RStudio viewer
#'   run.uncoverapp(where = "viewer")
#'   
#'   # Launch in browser
#'   run.uncoverapp(where = "browser")
#' }
#'
#' ## Only run this example in interactive R sessions
#' if (interactive()) {
#'   run.uncoverapp()
#' }
#'
#' @return No value is returned. The function launches the Shiny app.
#' 
#' @rawNamespace import(shiny, except = c(renderDataTable, dataTableOutput))
#' @rawNamespace import(Gviz, except = tags)
#' @import shinyWidgets
#' @import shinyBS
#' @import markdown
#' @import Homo.sapiens
#' @import openxlsx
#' @import stringr
#' @import condformat
#' @import processx
#' @import bedr
#' @import org.Hs.eg.db
#' @import TxDb.Hsapiens.UCSC.hg38.knownGene
#' @import TxDb.Hsapiens.UCSC.hg19.knownGene
#' @import EnsDb.Hsapiens.v86
#' @import EnsDb.Hsapiens.v75
#' @import OrganismDbi
#' @importFrom shinyjs useShinyjs hidden enable
#' @importFrom shinycssloaders withSpinner
#' @importFrom DT renderDataTable dataTableOutput
#' @importFrom dplyr filter mutate select arrange group_by summarise
#' @importFrom Rsamtools ScanBamParam PileupParam pileup
#' @importFrom rlist list.append
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom rappdirs user_cache_dir
#' 
#' @export
uncoverAPP <- function() {
  appDir <- system.file("shiny-dir", package = "uncoverappLib")
  if (appDir == "") {
    stop("Could not find shiny app directory. Try reinstalling uncoverappLib.", call. = FALSE)
  }
  shiny::runApp(appDir)
}


#' @title Launch uncoverApp in specific environment
#' @description
#' Controls the \code{shiny.launch.browser} option to launch uncoverApp
#' in an external web browser, the RStudio Viewer pane, or a new RStudio window.
#' 
#' @param where Character. One of \code{"browser"}, \code{"viewer"}, or \code{"window"}.
#'   If missing, the default shiny launch behavior will be used.
#'
#' @return No value is returned. The function launches the Shiny app.
#'
#' @examples
#' if (interactive()) {
#'   run.uncoverapp(where = "window")
#' }
#' 
#' @export
run.uncoverapp <- function(where = c("browser", "viewer", "window")) {
  if (missing(where)) {
    message("uncoverApp is launching with your default option.")
    uncoverAPP()
  } else {
    options(shiny.launch.browser = switch(
      match.arg(where, c("browser", "viewer", "window")),
      browser = get(".rs.invokeShinyWindowExternal", "tools:rstudio"),
      viewer = get(".rs.invokeShinyPaneViewer", "tools:rstudio"),
      window = get(".rs.invokeShinyWindowViewer", "tools:rstudio")
    ))
    uncoverAPP()
  }
}
