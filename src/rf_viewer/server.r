################################################################################
# server.R
################################################################################

server <- function(input, output, session) {
  presence_points_visible <- reactiveVal(FALSE)

  ##############################################################################
  # Aktueller Raster
  ##############################################################################

  current_raster <- reactive({
    req(input$layer)
    req(layers[[input$layer]])
    layers[[input$layer]]
  })

  current_label <- reactive({
    req(input$layer)
    unname(layer_labels[[input$layer]])
  })

  current_preview <- reactive({
    req(input$layer)

    if (!input$layer %in% names(preview_images)) {
      return(NULL)
    }

    preview_images[[input$layer]]
  })

  current_preview_path <- reactive({
    preview <- current_preview()

    if (is.null(preview)) {
      return(NULL)
    }

    if (isTRUE(presence_points_visible()) && "with_presence" %in% names(preview)) {
      return(preview[["with_presence"]])
    }

    if ("without_presence" %in% names(preview)) {
      return(preview[["without_presence"]])
    }

    preview[[1]]
  })

  ##############################################################################
  # Leaflet-Karte erzeugen
  ##############################################################################

  output$map <- renderLeaflet({
    req(current_raster())

    r <- current_raster()
    pal <- palette_fun(r)

    leaflet(options = leafletOptions(zoomControl = TRUE)) |>
      addProviderTiles(providers$Esri.WorldGrayCanvas) |>
      addRasterImage(
        r,
        colors = pal,
        opacity = 0.8,
        project = TRUE,
        layerId = "selected_raster"
      ) |>
      addLegend(
        pal = pal,
        values = values(r),
        title = current_label()
      ) |>
      fitBounds(
        xmin(r),
        ymin(r),
        xmax(r),
        ymax(r)
      )
  })

  ##############################################################################
  # Presence-Punkte in PNG-Vorschau ein-/ausblenden
  ##############################################################################

  observeEvent(input$layer, {
    presence_points_visible(FALSE)
  })

  observeEvent(input$toggle_presence_points, {
    presence_points_visible(!presence_points_visible())
  })

  output$presence_toggle_ui <- renderUI({
    preview <- current_preview()

    if (is.null(preview) || !"with_presence" %in% names(preview)) {
      return(NULL)
    }

    actionButton(
      inputId = "toggle_presence_points",
      label = if (isTRUE(presence_points_visible())) {
        "Presence-Punkte ausblenden"
      } else {
        "Presence-Punkte anzeigen"
      },
      width = "100%"
    )
  })

  output$preview_title <- renderUI({
    req(current_preview_path())

    tags$div(
      style = "margin-top: 18px;",
      h4(sprintf("%s - PNG-Vorschau", current_label()))
    )
  })

  output$heatmap_preview <- renderImage({
    path <- current_preview_path()
    req(path)

    list(
      src = path,
      contentType = "image/png",
      alt = sprintf("%s Vorschau", current_label()),
      width = "100%"
    )
  }, deleteFile = FALSE)

  ##############################################################################
  # Klick auf Karte
  ##############################################################################

  observeEvent(input$map_click, {
    click <- input$map_click
    req(click)

    lng <- click$lng
    lat <- click$lat

    raster_values <- extract_all_values(lng, lat)
    predictor_values <- extract_predictor_values(lng, lat)
    selected_value <- raster_values[[input$layer]]
    popup_predictor_names <- setdiff(names(predictor_values), "dem")
    predictor_popup_rows <- vapply(popup_predictor_names, function(name) {
      result <- predictor_values[[name]]
      hint <- if (identical(result$source, "nearest")) {
        " <small>(naechster gueltiger Wert)</small>"
      } else {
        ""
      }

      sprintf(
        "<b>%s:</b> %s%s",
        htmlEscape(unname(predictor_labels[[name]])),
        htmlEscape(fmt(result$value, unit = predictor_units[[name]])),
        hint
      )
    }, character(1))
    predictor_popup_text <- paste(predictor_popup_rows, collapse = "<br>")

    leafletProxy("map") |>
      clearMarkers() |>
      clearPopups() |>
      addCircleMarkers(
        lng = lng,
        lat = lat,
        radius = 6,
        stroke = TRUE,
        color = "red",
        fillOpacity = 1
      ) |>
      addPopups(
        lng = lng,
        lat = lat,
        popup = sprintf(
          paste(
            "<b>Koordinaten</b><br>",
            "Lat: %.5f<br>",
            "Lng: %.5f<br><br>",
            "<b>%s:</b> %s<br><br>",
            "<b>Praediktoren</b><br>",
            "%s"
          ),
          lat,
          lng,
          htmlEscape(current_label()),
          htmlEscape(fmt(selected_value)),
          predictor_popup_text
        )
      )

    output$coordinates <- renderText({
      sprintf(
        "Longitude: %.6f\nLatitude : %.6f",
        lng,
        lat
      )
    })

    output$info_table <- renderTable({
      data.frame(
        Kategorie = "Modell",
        Variable = unname(layer_labels[names(raster_values)]),
        Wert = vapply(raster_values, fmt, character(1)),
        Hinweis = "",
        check.names = FALSE
      )
    })
  })
}
