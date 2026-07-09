################################################################################
# ui.R
################################################################################

ui <- fluidPage(
  titlePanel("Viewer - Archaeological Prediction Map: Eisenzeit in Bayern"),

  sidebarLayout(
    sidebarPanel(
      width = 3,

      h4("Kartendarstellung"),

      selectInput(
        inputId = "layer",
        label = "Layer auswaehlen",
        choices = layer_choices,
        selected = default_layer
      ),

      uiOutput("presence_toggle_ui"),

      hr(),

      h4("Informationen"),

      p(
        "Alle .tif-Dateien unter output werden automatisch geladen.",
        "Klicke auf die Karte, um die Rasterwerte am Punkt auszulesen."
      ),

      hr(),

      tableOutput("info_table"),

      hr(),

      h4("Koordinaten"),

      verbatimTextOutput("coordinates"),

      hr(),

      h4("Legende"),

      helpText(
        "Die dargestellte Heatmap basiert auf dem ausgewaehlten GeoTIFF.",
        "Die Tabelle zeigt die Werte aller geladenen Output-Raster am Klickpunkt."
      )
    ),

    mainPanel(
      width = 9,

      leafletOutput(
        outputId = "map",
        height = "850px"
      ),

      uiOutput("preview_title"),

      imageOutput(
        outputId = "heatmap_preview",
        height = "auto"
      )
    )
  )
)
