#' Purpose: Shiny App for DESeqPLUS
#' Author: Katie Lee
#' Date: November 28, 2025
#'
#' This function launches the Shiny app for DESEqPLUS.
#'
#' @return launches Shiny app.
#'
#' @examples
#'
#' DESeqPLUS::runDESeqPLUS()
#'
#' @export
#' @importFrom shiny runApp

runDESeqPLUS <- function() {
  appLocation <- system.file("shiny-scripts", package = "DESeqPLUS")
  actionShiny <- shiny::runApp(appLocation, display.mode = "normal")
  return(actionShiny)
}

# [ END]
