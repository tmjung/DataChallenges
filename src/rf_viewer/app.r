################################################################################
# Viewer
# Archaeological Prediction Map - Eisenzeit in Bayern
################################################################################

library(shiny)

PROJECT_ROOT <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)

if (!dir.exists(file.path(PROJECT_ROOT, "output"))) {
  parent_root <- normalizePath(file.path(getwd(), "../.."), winslash = "/", mustWork = FALSE)

  if (dir.exists(file.path(parent_root, "output"))) {
    PROJECT_ROOT <- parent_root
  }
}

source(file.path(PROJECT_ROOT, "src/rf_viewer/global.r"), local = TRUE)
source(file.path(PROJECT_ROOT, "src/rf_viewer/ui.r"), local = TRUE)
source(file.path(PROJECT_ROOT, "src/rf_viewer/server.r"), local = TRUE)

shinyApp(ui = ui, server = server)
