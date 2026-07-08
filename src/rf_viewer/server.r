################################################################################
# server.R
################################################################################

server <- function(input, output, session) {

  ##############################################################################
  # Aktueller Raster (vorerst immer die Fundwahrscheinlichkeit)
  ##############################################################################

  current_raster <- reactive({

    prediction

  })

  ##############################################################################
  # Leaflet-Karte erzeugen
  ##############################################################################

  output$map <- renderLeaflet({

    req(current_raster())

    r <- current_raster()

    pal <- palette_fun(r)

    leaflet(options = leafletOptions(
      zoomControl = TRUE
    )) |>

      addProviderTiles(providers$Esri.WorldGrayCanvas) |>

      addRasterImage(
        r,
        colors = pal,
        opacity = 0.8,
        project = TRUE
      ) |>

      addLegend(
        pal = pal,
        values = values(r),
        title = "Fundwahrscheinlichkeit"
      ) |>

      fitBounds(
        xmin(r),
        ymin(r),
        xmax(r),
        ymax(r)
      )

  })

  ##############################################################################
  # TODO:
  # Hier später Layerwechsel einbauen
  #
  # observeEvent(input$layer,{
  #    ...
  # })
  ##############################################################################


  ##############################################################################
  # Klick auf Karte
  ##############################################################################

  observeEvent(input$map_click, {

    click <- input$map_click

    req(click)

    lng <- click$lng
    lat <- click$lat

    ###########################################################################
    # Werte auslesen
    ###########################################################################

    values <- extract_all_values(lng, lat)

    print("\n VALUES: \n")
    print(values)
    ###########################################################################
    # Marker setzen
    ###########################################################################

    leafletProxy("map") |>

      clearMarkers() |>

      addCircleMarkers(

        lng = lng,
        lat = lat,

        radius = 6,

        stroke = TRUE,

        color = "red",

        fillOpacity = 1

      )

    ###########################################################################
    # Koordinaten anzeigen
    ###########################################################################

    output$coordinates <- renderText({

      sprintf(

        "Longitude: %.6f\nLatitude : %.6f",

        lng,
        lat

      )

    })

    ###########################################################################
    # Tabelle erzeugen
    ###########################################################################

    output$info_table <- renderTable({

      data.frame(

        Variable = c(

          "Fundwahrscheinlichkeit",

          "Niederschlag",

          "Temperatur",

          "Elevation",

          "DEM",

          "Hangneigung",

          "Exposition"

        ),

        Wert = c(

          fmt(values$probability * 100, 1, " %"),

          fmt(values$precipitation, 1, " mm"),

          fmt(values$temperature, 1, " °C"),

          fmt(values$elevation, 1, " m"),

          fmt(values$dem, 1, " m"),

          fmt(values$slope, 1, " °"),

          fmt(values$aspect, 1, " °")

        ),

        check.names = FALSE

      )

    }
  )

##############################################################################
# Popup auf Karte anzeigen
##############################################################################

popup_text <- sprintf(
  paste(
    "<b>Koordinaten</b><br>",
    "Lat: %.5f<br>",
    "Lng: %.5f<br><br>",

    "<b>Fundwahrscheinlichkeit:</b> %s<br>",
    "<b>Niederschlag:</b> %s<br>",
    "<b>Temperatur:</b> %s<br>",
    "<b>Elevation:</b> %s<br>",
    "<b>DEM:</b> %s<br>",
    "<b>Hangneigung:</b> %s<br>",
    "<b>Exposition:</b> %s"
  ),

  lat,
  lng,

  fmt(values$probability * 100, 1, " %"),
  fmt(values$precipitation, 1, " mm"),
  fmt(values$temperature, 1, " °C"),
  fmt(values$elevation, 1, " m"),
  fmt(values$dem, 1, " m"),
  fmt(values$slope, 1, " °"),
  fmt(values$aspect, 1, " °")
)

leafletProxy("map") |>

  clearPopups() |>

  addPopups(
    lng = lng,
    lat = lat,
    popup = popup_text
  )

  })

}