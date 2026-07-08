################################################################################
# Random Forest Viewer
# Archaeological Prediction Map
################################################################################

library(shiny)

source("src/rf_viewer/global.r", local = TRUE)
source("src/rf_viewer/ui.r", local = TRUE)
source("src/rf_viewer/server.r", local = TRUE)

shinyApp(ui = ui, server = server)