#' Get Annotation Files
#'
#' Downloads or retrieves annotation files for the specified genome assembly
#'
#' @param assembly (character) Reference genome: "hg19", "hg38"
#' @param verbose (logical) print messages
#' @examples
#' getAnnotationFiles(assembly = "hg38", verbose = TRUE)

#' @return Named list with paths to annotation files for each assembly
#' @export
getAnnotationFiles <- function(assembly = c("hg19", "hg38"), verbose = FALSE) {
  
  # Match argument
  assembly <- match.arg(assembly, several.ok = TRUE)
  
  if (verbose) {
    cat("=== getAnnotationFiles ===\n")
    cat("Requested assemblies:", paste(assembly, collapse = ", "), "\n")
  }
  
  # Define download URLs and cache directory
  base_url <- "https://zenodo.org/records/17524340/files/"
  cache_dir <- rappdirs::user_cache_dir("uncoverapp")
  
  # Create cache directory if it doesn't exist
  if (!dir.exists(cache_dir)) {
    dir.create(cache_dir, recursive = TRUE)
    if (verbose) cat("Created cache directory:", cache_dir, "\n")
  }
  
  # File definitions
  files <- list(
    hg19 = c(
      gz = "sorted_hg19.bed.gz",
      tbi = "sorted_hg19.bed.gz.tbi"
    ),
    hg38 = c(
      gz = "sorted_hg38.bed.gz",
      tbi = "sorted_hg38.bed.gz.tbi"
    )
  )
  
  # Download function
  download_file <- function(assembly_name, file_type) {
    filename <- files[[assembly_name]][file_type]
    local_path <- file.path(cache_dir, filename)
    
    # Check if file already exists
    if (file.exists(local_path)) {
      if (verbose) cat("File already exists:", local_path, "\n")
      return(local_path)
    }
    
    # Download
    url <- paste0(base_url, filename)
    if (verbose) cat("Downloading:", filename, "\n")
    
    tryCatch({
      # *** ADDED: method = "libcurl" for reliability ***
      download.file(url, local_path, mode = "wb", quiet = !verbose, timeout = 10000000, method = "libcurl")
      if (verbose) cat("Downloaded to:", local_path, "\n")
      return(local_path)
    }, error = function(e) {
      warning("Failed to download ", filename, ": ", e$message)
      return(NULL)
    })
  }
  
  # Download files for requested assemblies
  result <- list()
  
  for (asm in assembly) {
    if (verbose) cat("\nProcessing assembly:", asm, "\n")
    
    # Download both .gz and .tbi files
    gz_path <- download_file(asm, "gz")
    tbi_path <- download_file(asm, "tbi")
    
    if (!is.null(gz_path) && !is.null(tbi_path)) {
      result[[asm]] <- c(gz_path, tbi_path)
    } else {
      warning("Failed to download files for ", asm)
      result[[asm]] <- NULL
    }
  }
  
  if (verbose) {
    cat("\n=== Download complete ===\n")
    cat("Files available:\n")
    print(result)
  }
  
  return(result)
}


