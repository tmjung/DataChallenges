################################################################################
# ui.R
################################################################################

ui <- fluidPage(

  titlePanel("Random Forest Viewer - Archaeological Prediction"),

  sidebarLayout(

    ###########################################################################
    # Sidebar
    ###########################################################################

    sidebarPanel(

      width = 3,

      h4("Kartendarstellung"),

      selectInput(
        inputId = "layer",
        label = "Heatmap auswählen",
        choices = names(layers),
        selected = "Fundwahrscheinlichkeit"
      ),

      hr(),

      h4("Informationen"),

      p(
        "Klicke auf einen beliebigen Punkt der Karte, um die Werte der",
        "Umweltvariablen und die Fundwahrscheinlichkeit anzuzeigen."
      ),

      hr(),

      tableOutput("info_table"),

      hr(),

      h4("Koordinaten"),

      verbatimTextOutput("coordinates"),

      hr(),

      h4("Legende"),

      helpText(
        "Die dargestellte Heatmap basiert auf dem ausgewählten GeoTIFF.",
        "Alle Werte werden direkt aus den Rasterdaten ausgelesen."
      )

    ),

    ###########################################################################
    # Karte
    ###########################################################################

    mainPanel(

      width = 9,

      leafletOutput(
        outputId = "map",
        height = "850px"
      )

    )

  )

)