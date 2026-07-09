# =============================================================================
# Ensemble-Heatmaps aus Random-Forest- und MaxEnt-GeoTIFF-Vorhersagen
# =============================================================================

library(terra)
library(ggplot2)
library(viridis)
library(readr)

# Usage:
# Rscript src/rf_me_mean_heatmap.r [maxent_tif] [rf_tif] [output_prefix] \
#   [rf_weight] [maxent_weight] [output_dir]
#
# Defaults:
# maxent_tif    = output/maxent_rf_grid/maxent_rf_grid_mittlere_standorteignung.tif, if available
# rf_tif        = output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif, if available
# output_prefix = rf_maxent
# rf_weight     = 1
# maxent_weight = 1
# output_dir    = output/ensemble

args <- commandArgs(trailingOnly = TRUE)

default_maxent_tif <- ifelse(
  file.exists("output/maxent_rf_grid/maxent_rf_grid_mittlere_standorteignung.tif"),
  "output/maxent_rf_grid/maxent_rf_grid_mittlere_standorteignung.tif",
  ifelse(
    file.exists("maxent_repeated_mean_prediction_rf_grid.tif"),
    "maxent_repeated_mean_prediction_rf_grid.tif",
    "maxent_repeated_mean_prediction.tif"
  )
)

maxent_tif <- ifelse(length(args) >= 1, args[1], default_maxent_tif)

default_rf_tif <- ifelse(
  file.exists("output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif"),
  "output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif",
  "rf_final_5km_blue_probability.tif"
)

rf_tif <- ifelse(length(args) >= 2, args[2], default_rf_tif)
output_prefix <- ifelse(length(args) >= 3, args[3], "rf_maxent")
rf_weight <- ifelse(length(args) >= 4, as.numeric(args[4]), 1)
maxent_weight <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
output_dir <- ifelse(length(args) >= 6, args[6], "output/ensemble")

if (!file.exists(maxent_tif)) {
  stop("MaxEnt-GeoTIFF nicht gefunden: ", maxent_tif)
}

if (!file.exists(rf_tif)) {
  stop("Random-Forest-GeoTIFF nicht gefunden: ", rf_tif)
}

if (is.na(rf_weight) || is.na(maxent_weight) || rf_weight < 0 || maxent_weight < 0) {
  stop("Die Gewichte muessen nicht-negative Zahlen sein.")
}

if (rf_weight == 0 && maxent_weight == 0) {
  stop("Mindestens ein Modellgewicht muss groesser als 0 sein.")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mean_tif <- file.path(output_dir, paste0(output_prefix, "_weighted_mean.tif"))
mean_png <- file.path(output_dir, paste0(output_prefix, "_weighted_mean.png"))
mean_points_png <- file.path(output_dir, paste0(output_prefix, "_weighted_mean_presence_points.png"))
difference_tif <- file.path(output_dir, paste0(output_prefix, "_rf_minus_maxent.tif"))
difference_png <- file.path(output_dir, paste0(output_prefix, "_rf_minus_maxent.png"))
agreement_tif <- file.path(output_dir, paste0(output_prefix, "_agreement_sd.tif"))
agreement_png <- file.path(output_dir, paste0(output_prefix, "_agreement_sd.png"))
diagnostics_csv <- file.path(output_dir, paste0(output_prefix, "_raster_diagnostics.csv"))

cat("Lade Raster...\n")
maxent_raster <- rast(maxent_tif)[[1]]
rf_raster <- rast(rf_tif)[[1]]

names(maxent_raster) <- "maxent"
names(rf_raster) <- "rf"

raster_summary <- function(r, label, path) {
  e <- ext(r)

  data.frame(
    raster = label,
    file = path,
    nrow = nrow(r),
    ncol = ncol(r),
    nlyr = nlyr(r),
    res_x = res(r)[1],
    res_y = res(r)[2],
    xmin = xmin(e),
    xmax = xmax(e),
    ymin = ymin(e),
    ymax = ymax(e),
    crs = crs(r, describe = TRUE)$code,
    valid_cells = global(!is.na(r), "sum", na.rm = TRUE)[1, 1]
  )
}

cat("\nRasterdiagnose vor der Angleichung:\n")
diagnostics_before <- rbind(
  raster_summary(rf_raster, "rf", rf_tif),
  raster_summary(maxent_raster, "maxent", maxent_tif)
)
print(diagnostics_before)

same_crs <- same.crs(maxent_raster, rf_raster)
same_grid <- compareGeom(
  maxent_raster,
  rf_raster,
  stopOnError = FALSE
)

cat("\nVergleichbarkeit vor der Angleichung:\n")
cat("Gleiches CRS:", same_crs, "\n")
cat("Gleiches Raster/Ausdehnung/Aufloesung:", same_grid, "\n")

# Use the RF raster as target. If train_maxent_5.r was used, this should already
# match exactly and no resampling should happen.
cat("\nGleiche MaxEnt-Raster bei Bedarf an das RF-Raster an...\n")

if (!same_crs) {
  maxent_aligned <- project(
    maxent_raster,
    rf_raster,
    method = "bilinear"
  )
} else if (!same_grid) {
  maxent_aligned <- resample(
    maxent_raster,
    rf_raster,
    method = "bilinear"
  )
} else {
  maxent_aligned <- maxent_raster
}

names(maxent_aligned) <- "maxent"

diagnostics_after <- rbind(
  raster_summary(rf_raster, "rf_aligned", rf_tif),
  raster_summary(maxent_aligned, "maxent_aligned", maxent_tif)
)

prediction_stack <- c(rf_raster, maxent_aligned)
names(prediction_stack) <- c("rf", "maxent")

valid_both <- !is.na(rf_raster) & !is.na(maxent_aligned)
valid_any <- !is.na(rf_raster) | !is.na(maxent_aligned)

common_valid_cells <- global(valid_both, "sum", na.rm = TRUE)[1, 1]
any_valid_cells <- global(valid_any, "sum", na.rm = TRUE)[1, 1]

cat("\nRasterdiagnose nach der Angleichung:\n")
print(diagnostics_after)
cat("Zellen mit Werten beider Modelle:", common_valid_cells, "\n")
cat("Zellen mit mindestens einem Modellwert:", any_valid_cells, "\n")

write_csv(
  rbind(diagnostics_before, diagnostics_after),
  diagnostics_csv
)

cat("\nBerechne gewichteten Mittelwert...\n")
cat("RF-Gewicht:", rf_weight, "\n")
cat("MaxEnt-Gewicht:", maxent_weight, "\n")

weighted_mean_fun <- function(rf, maxent) {
  numerator <- ifelse(is.na(rf), 0, rf * rf_weight) +
    ifelse(is.na(maxent), 0, maxent * maxent_weight)

  denominator <- ifelse(is.na(rf), 0, rf_weight) +
    ifelse(is.na(maxent), 0, maxent_weight)

  ifelse(denominator == 0, NA, numerator / denominator)
}

mean_raster <- lapp(
  prediction_stack,
  fun = weighted_mean_fun
)

names(mean_raster) <- "rf_maxent_weighted_mean"

difference_raster <- rf_raster - maxent_aligned
names(difference_raster) <- "rf_minus_maxent"

agreement_raster <- abs(difference_raster) / sqrt(2)
names(agreement_raster) <- "agreement_sd"

cat("\nWertebereiche der Ausgaben:\n")
print(global(mean_raster, "range", na.rm = TRUE))
print(global(difference_raster, "range", na.rm = TRUE))
print(global(agreement_raster, "range", na.rm = TRUE))

cat("\nSpeichere GeoTIFFs...\n")
writeRaster(mean_raster, mean_tif, overwrite = TRUE)
writeRaster(difference_raster, difference_tif, overwrite = TRUE)
writeRaster(agreement_raster, agreement_tif, overwrite = TRUE)

raster_extent <- ext(mean_raster)

plot_raster_png <- function(r, png_path, title, fill_name, option, limits = NULL) {
  plot_df <- as.data.frame(
    r,
    xy = TRUE,
    na.rm = FALSE
  )

  value_col <- names(r)[1]

  p <- ggplot(plot_df) +
    geom_raster(
      aes(
        x = x,
        y = y,
        fill = .data[[value_col]]
      )
    ) +
    scale_fill_viridis_c(
      option = option,
      limits = limits,
      na.value = "transparent",
      name = fill_name
    ) +
    coord_quickmap(
      xlim = c(xmin(raster_extent), xmax(raster_extent)),
      ylim = c(ymin(raster_extent), ymax(raster_extent)),
      expand = FALSE
    ) +
    labs(
      title = title,
      subtitle = paste(
        "RF:",
        basename(rf_tif),
        "| MaxEnt:",
        basename(maxent_tif)
      ),
      x = "Laengengrad",
      y = "Breitengrad"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 9),
      legend.position = "right"
    )

  ggsave(
    png_path,
    plot = p,
    width = 12,
    height = 9,
    dpi = 300
  )

  p
}

cat("Speichere PNG-Heatmaps...\n")

p_mean <- plot_raster_png(
  mean_raster,
  mean_png,
  "Gewichtete Mittelwert-Heatmap: Random Forest + MaxEnt",
  "Gewichteter Mittelwert",
  "magma",
  limits = c(0, 1)
)

plot_raster_png(
  difference_raster,
  difference_png,
  "Differenz-Heatmap: Random Forest - MaxEnt",
  "RF - MaxEnt",
  "viridis"
)

plot_raster_png(
  agreement_raster,
  agreement_png,
  "Modell-Unstimmigkeit: Standardabweichung zwischen Random Forest und MaxEnt",
  "SD",
  "magma",
  limits = c(0, 1)
)

presence_csv <- "data/ffm_vfpa_eisenzeit.csv"

if (file.exists(presence_csv)) {
  cat("Speichere PNG-Heatmap mit Presence-Punkten...\n")

  presence_raw <- read_csv(
    presence_csv,
    locale = locale(decimal_mark = ","),
    col_types = cols(
      lng_wgs84 = col_character(),
      lat_wgs84 = col_character()
    ),
    show_col_types = FALSE
  )

  presence_df <- data.frame(
    lng = as.numeric(gsub(",", ".", presence_raw$lng_wgs84)),
    lat = as.numeric(gsub(",", ".", presence_raw$lat_wgs84))
  )

  presence_df <- presence_df[complete.cases(presence_df), ]

  p_mean_points <- p_mean +
    geom_point(
      data = presence_df,
      aes(x = lng, y = lat),
      inherit.aes = FALSE,
      color = "dodgerblue3",
      fill = "white",
      shape = 21,
      stroke = 0.25,
      size = 0.8,
      alpha = 0.85
    ) +
    labs(
      title = "Gewichtete Mittelwert-Heatmap mit Presence-Punkten"
    )

  ggsave(
    mean_points_png,
    plot = p_mean_points,
    width = 12,
    height = 9,
    dpi = 300
  )
} else {
  cat("Presence-CSV nicht gefunden, Presence-Punkt-PNG uebersprungen: ", presence_csv, "\n", sep = "")
}

cat("\nGespeicherte Ausgaben:\n")
cat(" - ", diagnostics_csv, "\n", sep = "")
cat(" - ", mean_tif, "\n", sep = "")
cat(" - ", mean_png, "\n", sep = "")
cat(" - ", mean_points_png, "\n", sep = "")
cat(" - ", difference_tif, "\n", sep = "")
cat(" - ", difference_png, "\n", sep = "")
cat(" - ", agreement_tif, "\n", sep = "")
cat(" - ", agreement_png, "\n", sep = "")
