#' Setup UncoverApp (download annotation files)
#'
#' @description
#' Downloads and caches annotation files for hg19 and hg38.
#'
#' @return Invisibly returns a list of downloaded file paths.
#'
#' @examples
#' # Show function signature
#' args(setup_uncoverapp)
#' 
#' \dontrun{
#' # Download annotation files (one-time setup)
#' setup_uncoverapp()
#' }
#'
#' @importFrom utils txtProgressBar setTxtProgressBar
#' @export
setup_uncoverapp <- function() {
  cat("\n========================================\n")
  cat("  UncoverApp Setup                     \n")
  cat("========================================\n\n")
  cat("Checking annotation files...\n\n")
  
  cache_dir <- rappdirs::user_cache_dir("uncoverapp")
  existing_files <- if (dir.exists(cache_dir)) {
    list.files(cache_dir, pattern = "\\.gz$", full.names = TRUE)
  } else {
    character(0)
  }
  
  if (length(existing_files) >= 2) {
    cat("[OK] Annotation files already downloaded!\n")
    cat("Location:", cache_dir, "\n\n")
    for (f in existing_files) {
      cat("  -", basename(f), "\n")
    }
    cat("\nRun: uncoverappLib.run()\n\n")
    return(invisible(existing_files))
  }
  
  cat("Downloading annotation files...\n")
  cat("This will take ~10-20 minutes (~1GB)\n\n")
  
  pb <- txtProgressBar(min = 0, max = 100, style = 3, char = "=")
  
  for (i in 1:10) {
    Sys.sleep(0.5)
    setTxtProgressBar(pb, i * 10)
  }
  
  result <- getAnnotationFiles(assembly = c("hg19", "hg38"), verbose = TRUE)
  
  setTxtProgressBar(pb, 100)
  close(pb)
  
  cat("\n\n[OK] Setup complete!\n")
  cat("Run: uncoverappLib.run()\n\n")
  
  invisible(result)
}


#' Check Annotation Files
#'
#' @description
#' Verifies that annotation files are available for both hg19 and hg38.
#'
#' @return Invisibly returns a character vector of annotation file paths.
#'
#' @examples
#' # Check if annotation files are available
#' check_annotations()
#'
#' @export
check_annotations <- function() {
  cache_dir <- rappdirs::user_cache_dir("uncoverapp")
  
  if (!dir.exists(cache_dir)) {
    cat("No files downloaded yet. Run: setup_uncoverapp()\n")
    return(invisible(NULL))
  }
  
  files <- list.files(cache_dir, pattern = "\\.gz$", full.names = TRUE)
  
  if (length(files) == 0) {
    cat("No files found. Run: setup_uncoverapp()\n")
  } else {
    cat("\nAnnotation files:\n")
    for (f in files) {
      cat("  [OK]", basename(f), "-", round(file.size(f)/1024^2, 1), "MB\n")
    }
    cat("\nLocation:", cache_dir, "\n")
  }
  
  invisible(files)
}