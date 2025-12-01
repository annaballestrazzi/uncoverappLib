# ============================================================================
# waiter-helpers.R - Helper functions per Waiter/Waitress con logo uncoverAPP
# ============================================================================

#' Show uncoverAPP branded waiter
#' 
#' @param message Main message to display
#' @param detail Optional detail message
#' @param color Background color (default: uncoverAPP blue)
#' @param logo_height Height of logo in pixels (default: 120)
show_uncoverapp_waiter <- function(
    message = "Processing...", 
    detail = NULL,
    color = "rgba(27, 79, 114, 0.95)",  # uncoverAPP dark blue
    logo_height = 120
) {
  
  # Logo uncoverAPP (SVG con sfondo trasparente)
  # Il logo deve essere in inst/www/logo.svg
  logo_html <- shiny::tags$div(
    shiny::tags$img(
      src = "logo.svg",  # Path relativo a www/
      alt = "uncoverAPP",
      style = sprintf(
        "height: %dpx; margin-bottom: 30px; filter: drop-shadow(0 2px 8px rgba(0,0,0,0.3));",
        logo_height
      )
    )
  )
  
  waiter_html <- shiny::tagList(
    # Logo
    logo_html,
    
    # Spinner animato
    waiter::spin_folding_cube(),
    
    # Messaggio principale
    shiny::h3(message, style = "color: white; margin-top: 30px; font-weight: 600;"),
    
    # Dettaglio opzionale
    if (!is.null(detail)) {
      shiny::p(detail, style = "color: #E8F4F8; font-size: 16px; margin-top: 10px;")
    }
  )
  
  waiter::waiter_show(
    html = waiter_html,
    color = color
  )
}


#' Create uncoverAPP branded waitress (progress bar)
#' 
#' @param theme Waitress theme (default: "overlay")
#' @param min Minimum value (default: 0)
#' @param max Maximum value (default: 100)
#' @param infinite Use infinite progress bar (default: FALSE)
#' @return Waitress object
create_uncoverapp_waitress <- function(
    theme = "overlay",
    min = 0,
    max = 100,
    infinite = FALSE
) {
  
  waitress <- waiter::Waitress$new(
    theme = theme,
    min = min,
    max = max,
    infinite = infinite
  )
  
  # Return the waitress object without modifying it
  # (customization can be done when calling the methods)
  return(waitress)
}


#' Show progress with steps
#' 
#' @param waitress Waitress object
#' @param step_name Name of current step
#' @param increment Progress increment (percentage)
show_progress_step <- function(waitress, step_name, increment) {
  waitress$inc(increment)
  waitress$notify(paste0("â³ ", step_name))
  Sys.sleep(0.1)  # Brief pause for UX
}


#' Wrapper for plot generation with progress tracking
#' 
#' @param plot_function Function that generates the plot
#' @param steps List of step names and their percentage contributions
#' @return Plot output
generate_plot_with_progress <- function(plot_function, steps = NULL) {
  
  if (is.null(steps)) {
    steps <- list(
      "Preparing data" = 20,
      "Creating gene tracks" = 40,
      "Rendering plot" = 40
    )
  }
  
  waitress <- create_uncoverapp_waitress()
  waitress$start()
  
  tryCatch({
    result <- NULL
    step_idx <- 0
    
    for (step_name in names(steps)) {
      step_idx <- step_idx + 1
      show_progress_step(waitress, step_name, steps[[step_name]])
      
      # Se Ã¨ l'ultimo step, genera il plot
      if (step_idx == length(steps)) {
        result <- plot_function()
      }
    }
    
    waitress$close()
    return(result)
    
  }, error = function(e) {
    waitress$close()
    showNotification(
      paste("Error generating plot:", e$message), 
      type = "error", 
      duration = 5
    )
    return(NULL)
  })
}


#' Wrapper for download with progress tracking
#' 
#' @param download_function Function that performs the download
#' @param filename Output filename
#' @param steps List of step names and their percentage contributions
download_with_progress <- function(download_function, filename, steps = NULL) {
  
  if (is.null(steps)) {
    steps <- list(
      "Loading data" = 15,
      "Preparing export" = 20,
      "Creating file" = 30,
      "Applying formatting" = 25,
      "Saving file" = 10
    )
  }
  
  waitress <- create_uncoverapp_waitress()
  waitress$start()
  
  tryCatch({
    result <- download_function(waitress, steps)
    
    waitress$close()
    showNotification(
      paste0("âœ“ File saved: ", basename(filename)), 
      type = "message", 
      duration = 3
    )
    
    return(result)
    
  }, error = function(e) {
    waitress$close()
    showNotification(
      paste("Download error:", e$message), 
      type = "error", 
      duration = 5
    )
    return(NULL)
  })
}