################################################################################
# global.R
# Laedt alle GeoTIFF-Dateien aus output und stellt Hilfsfunktionen bereit.
################################################################################

library(shiny)
library(leaflet)
library(terra)
library(sf)
library(viridis)
library(htmltools)

################################################################################
# Datenordner
################################################################################

PROJECT_ROOT <- if (exists("PROJECT_ROOT")) {
  PROJECT_ROOT
} else {
  normalizePath(getwd(), winslash = "/", mustWork = TRUE)
}

OUTPUT_DIR <- file.path(PROJECT_ROOT, "output")
DATA_DIR <- file.path(PROJECT_ROOT, "data")

################################################################################
# Raster finden und laden
################################################################################

find_output_tifs <- function(output_dir = OUTPUT_DIR) {
  if (!dir.exists(output_dir)) {
    warning(sprintf("Output-Ordner nicht gefunden: %s", output_dir))
    return(character())
  }

  sort(
    list.files(
      output_dir,
      pattern = "\\.tif$",
      recursive = TRUE,
      full.names = TRUE,
      ignore.case = TRUE
    )
  )
}

make_layer_name <- function(path) {
  rel_path <- normalizePath(path, winslash = "/", mustWork = FALSE)
  output_path <- normalizePath(OUTPUT_DIR, winslash = "/", mustWork = FALSE)
  rel_path <- sub(paste0("^", output_path, "/?"), "", rel_path)
  sub("\\.tif$", "", rel_path, ignore.case = TRUE)
}

make_display_name <- function(layer_name) {
  display_names <- c(
    "ensemble/rf_maxent_agreement_sd" = "RF + MaxEnt Unsicherheit",
    "ensemble/rf_maxent_rf_minus_maxent" = "RF vs. MaxEnt Differenz",
    "ensemble/rf_maxent_weighted_mean" = "RF + MaxEnt Ensemble",
    "ensemble_glm/rf_maxent_glm_agreement_sd" = "RF + MaxEnt + GLM Unsicherheit",
    "ensemble_glm/rf_maxent_glm_maxent_minus_glm" = "MaxEnt vs. GLM Differenz",
    "ensemble_glm/rf_maxent_glm_rf_minus_glm" = "RF vs. GLM Differenz",
    "ensemble_glm/rf_maxent_glm_rf_minus_maxent" = "RF vs. MaxEnt Differenz",
    "ensemble_glm/rf_maxent_glm_weighted_mean" = "RF + MaxEnt + GLM Ensemble",
    "glm/glm_probability" = "GLM Prediction",
    "maxent/maxent_repeated_mean_prediction" = "MaxEnt Prediction",
    "maxent/maxent_repeated_mean_prediction_rf_grid" = "MaxEnt Prediction (RF Grid)",
    "maxent/maxent_repeated_uncertainty" = "MaxEnt Unsicherheit",
    "maxent/maxent_repeated_uncertainty_rf_grid" = "MaxEnt Unsicherheit (RF Grid)",
    "random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit" = "Random Forest Prediction"
  )

  if (layer_name %in% names(display_names)) {
    return(unname(display_names[[layer_name]]))
  }

  fallback <- basename(layer_name)
  fallback <- gsub("_", " ", fallback)
  tools::toTitleCase(fallback)
}

load_raster <- function(path) {
  if (!file.exists(path)) {
    warning(sprintf("Raster nicht gefunden: %s", path))
    return(NULL)
  }

  tryCatch(
    rast(path),
    error = function(e) {
      warning(sprintf("Fehler beim Laden von %s: %s", path, e$message))
      NULL
    }
  )
}

make_single_layer <- function(raster, aggregation, layer_name) {
  if (is.null(raster)) {
    return(NULL)
  }

  if (terra::nlyr(raster) > 1) {
    raster <- switch(
      aggregation,
      sum = sum(raster, na.rm = TRUE),
      mean = mean(raster, na.rm = TRUE),
      raster[[1]]
    )
  }

  names(raster) <- layer_name
  raster
}

tif_files <- find_output_tifs()
layers <- stats::setNames(lapply(tif_files, load_raster), make_layer_name(tif_files))
layers <- layers[!vapply(layers, is.null, logical(1))]
layer_labels <- stats::setNames(vapply(names(layers), make_display_name, character(1)), names(layers))
layer_choices <- stats::setNames(names(layers), layer_labels)

################################################################################
# Praediktor-Raster
################################################################################

predictor_specs <- list(
  precipitation = list(
    label = "Niederschlag",
    path = file.path(DATA_DIR, "climate", "precipitation.tif"),
    aggregation = "sum",
    unit = " mm"
  ),
  temperature = list(
    label = "Temperatur",
    path = file.path(DATA_DIR, "climate", "temperature.tif"),
    aggregation = "mean",
    unit = " deg C"
  ),
  elevation = list(
    label = "Elevation",
    path = file.path(DATA_DIR, "dem", "DEU_elv_msk.tif"),
    aggregation = "first",
    unit = " m"
  ),
  dem = list(
    label = "DEM",
    path = file.path(DATA_DIR, "dem", "dem.tif"),
    aggregation = "first",
    unit = " m"
  ),
  slope = list(
    label = "Hangneigung",
    path = file.path(DATA_DIR, "dem", "slope.tif"),
    aggregation = "first",
    unit = " deg"
  ),
  aspect = list(
    label = "Exposition",
    path = file.path(DATA_DIR, "dem", "aspect.tif"),
    aggregation = "first",
    unit = " deg"
  )
)

load_predictor <- function(name, spec) {
  raster <- load_raster(spec$path)
  make_single_layer(raster, spec$aggregation, name)
}

predictors <- mapply(
  load_predictor,
  names(predictor_specs),
  predictor_specs,
  SIMPLIFY = FALSE
)
predictors <- predictors[!vapply(predictors, is.null, logical(1))]
predictor_labels <- vapply(predictor_specs[names(predictors)], `[[`, character(1), "label")
predictor_units <- vapply(predictor_specs[names(predictors)], `[[`, character(1), "unit")

default_layer <- if ("random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit" %in% names(layers)) {
  "random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit"
} else if (length(layers) > 0) {
  names(layers)[1]
} else {
  NULL
}

################################################################################
# PNG-Vorschau mit optionalen Presence-Punkten
################################################################################

preview_images <- list(
  "glm/glm_probability" = list(
    without_presence = file.path(OUTPUT_DIR, "glm", "glm_probability_heatmap.png"),
    with_presence = file.path(OUTPUT_DIR, "glm", "glm_fundwahrscheinlichkeit_mit_presence_punkten.png")
  ),
  "random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit" = list(
    without_presence = file.path(OUTPUT_DIR, "random_forest_5km_buffer", "random_forest_5km_buffer_fundwahrscheinlichkeit_heatmap.png"),
    with_presence = file.path(OUTPUT_DIR, "random_forest_5km_buffer", "random_forest_5km_buffer_fundwahrscheinlichkeit_mit_presence_punkten.png")
  ),
  "ensemble/rf_maxent_weighted_mean" = list(
    without_presence = file.path(OUTPUT_DIR, "ensemble", "rf_maxent_weighted_mean.png"),
    with_presence = file.path(OUTPUT_DIR, "ensemble", "rf_maxent_weighted_mean_presence_points.png")
  ),
  "ensemble_glm/rf_maxent_glm_weighted_mean" = list(
    without_presence = file.path(OUTPUT_DIR, "ensemble_glm", "rf_maxent_glm_weighted_mean.png"),
    with_presence = file.path(OUTPUT_DIR, "ensemble_glm", "rf_maxent_glm_weighted_mean_presence_points.png")
  )
)

preview_images <- lapply(preview_images, function(paths) {
  existing_paths <- paths[file.exists(unlist(paths))]

  if (length(existing_paths) == 0) {
    return(NULL)
  }

  existing_paths
})
preview_images <- preview_images[!vapply(preview_images, is.null, logical(1))]

################################################################################
# Bounding Box
################################################################################

if (length(layers) > 0) {
  bbox <- ext(layers[[default_layer]])
} else {
  bbox <- ext(
    8.9,
    13.9,
    47.2,
    50.6
  )
}

################################################################################
# Farbpalette
################################################################################

palette_fun <- function(r) {
  colorNumeric(
    palette = viridis(256, option = "magma"),
    domain = values(r),
    na.color = "transparent"
  )
}

################################################################################
# Rasterwert auslesen
################################################################################

extract_value_detail <- function(raster, lng, lat, use_nearest = FALSE, search_radius = 50000) {
  if (is.null(raster)) {
    return(list(value = NA_real_, source = "missing"))
  }

  tryCatch({
    point <- terra::vect(
      data.frame(lng = lng, lat = lat),
      geom = c("lng", "lat"),
      crs = "EPSG:4326"
    )

    raster_crs <- terra::crs(raster)

    if (!is.na(raster_crs) && nzchar(raster_crs)) {
      point <- terra::project(point, raster_crs)
    }

    v <- terra::extract(raster, point)

    value <- extract_first_numeric_value(v)

    if (!is.na(value)) {
      return(list(value = value, source = "direct"))
    }

    if (!isTRUE(use_nearest)) {
      return(list(value = NA_real_, source = "empty"))
    }

    nearest <- terra::extract(raster, point, search_radius = search_radius)
    nearest_value <- extract_first_numeric_value(nearest)

    if (!is.na(nearest_value)) {
      return(list(value = nearest_value, source = "nearest"))
    }

    list(value = NA_real_, source = "empty")
  }, error = function(e) {
    list(value = NA_real_, source = "error")
  })
}

extract_first_numeric_value <- function(extracted) {
  if (is.null(extracted) || nrow(extracted) == 0) {
    return(NA_real_)
  }

  value_columns <- setdiff(names(extracted), "ID")

  if (length(value_columns) == 0) {
    return(NA_real_)
  }

  value <- suppressWarnings(as.numeric(extracted[1, value_columns[1]]))

  if (length(value) == 0 || is.na(value)) {
    return(NA_real_)
  }

  value
}

extract_value <- function(raster, lng, lat, use_nearest = FALSE) {
  extract_value_detail(raster, lng, lat, use_nearest = use_nearest)$value
}

################################################################################
# Zahlen formatieren
################################################################################

fmt <- function(x, digits = 3, unit = "") {
  if (is.na(x)) {
    return("n/a")
  }

  paste0(
    format(round(x, digits), nsmall = digits),
    unit
  )
}

################################################################################
# Alle Rasterwerte eines Punktes
################################################################################

extract_all_values <- function(lng, lat) {
  vapply(layers, extract_value, numeric(1), lng = lng, lat = lat)
}

extract_predictor_values <- function(lng, lat) {
  lapply(predictors, extract_value_detail, lng = lng, lat = lat, use_nearest = TRUE)
}

################################################################################
# Startinformationen
################################################################################

cat("\n")
cat("=========================================\n")
cat(" Random Forest Viewer\n")
cat("=========================================\n")
cat(sprintf("Output-Raster gefunden: %d\n", length(tif_files)))
cat(sprintf("Output-Raster geladen:   %d\n", length(layers)))
cat(sprintf("Praediktor-Raster geladen: %d\n", length(predictors)))

for (n in names(layers)) {
  cat(sprintf("  OK  %s\n", n))
}

if (length(layers) == 0) {
  cat(sprintf("  Keine .tif-Dateien unter %s gefunden.\n", OUTPUT_DIR))
}

cat("=========================================\n\n")
