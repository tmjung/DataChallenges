# =============================================================================
# Ensemble-Heatmaps aus Random-Forest-, MaxEnt- und GLM-GeoTIFF-Vorhersagen
# =============================================================================

library(terra)
library(ggplot2)
library(viridis)
library(readr)

# Usage:
# Rscript src/rf_me_glm_heatmap.r [maxent_tif] [rf_tif] [glm_tif] \
#   [output_prefix] [rf_weight] [maxent_weight] [glm_weight] [output_dir]
#
# Defaults:
# maxent_tif    = maxent_repeated_mean_prediction_rf_grid.tif, if available
# rf_tif        = output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif, if available
# glm_tif       = output/glm/glm_probability.tif
# output_prefix = rf_maxent_glm
# rf_weight     = 1
# maxent_weight = 1
# glm_weight    = 1
# output_dir    = output/ensemble_glm

args <- commandArgs(trailingOnly = TRUE)

default_maxent_tif <- ifelse(
  file.exists("maxent_repeated_mean_prediction_rf_grid.tif"),
  "maxent_repeated_mean_prediction_rf_grid.tif",
  "maxent_repeated_mean_prediction.tif"
)

default_rf_tif <- ifelse(
  file.exists("output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif"),
  "output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif",
  "rf_final_5km_blue_probability.tif"
)

maxent_tif <- ifelse(length(args) >= 1, args[1], default_maxent_tif)
rf_tif <- ifelse(length(args) >= 2, args[2], default_rf_tif)
glm_tif <- ifelse(length(args) >= 3, args[3], "output/glm/glm_probability.tif")
output_prefix <- ifelse(length(args) >= 4, args[4], "rf_maxent_glm")
rf_weight <- ifelse(length(args) >= 5, as.numeric(args[5]), 1)
maxent_weight <- ifelse(length(args) >= 6, as.numeric(args[6]), 1)
glm_weight <- ifelse(length(args) >= 7, as.numeric(args[7]), 1)
output_dir <- ifelse(length(args) >= 8, args[8], "output/ensemble_glm")

if (!file.exists(maxent_tif)) {
  stop("MaxEnt-GeoTIFF nicht gefunden: ", maxent_tif)
}

if (!file.exists(rf_tif)) {
  stop("Random-Forest-GeoTIFF nicht gefunden: ", rf_tif)
}

if (!file.exists(glm_tif)) {
  stop("GLM-GeoTIFF nicht gefunden: ", glm_tif)
}

weights <- c(rf = rf_weight, maxent = maxent_weight, glm = glm_weight)

if (any(is.na(weights)) || any(weights < 0)) {
  stop("Die Gewichte muessen nicht-negative Zahlen sein.")
}

if (sum(weights) == 0) {
  stop("Mindestens ein Modellgewicht muss groesser als 0 sein.")
}

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

mean_tif <- file.path(output_dir, paste0(output_prefix, "_weighted_mean.tif"))
mean_png <- file.path(output_dir, paste0(output_prefix, "_weighted_mean.png"))
mean_points_png <- file.path(output_dir, paste0(output_prefix, "_weighted_mean_presence_points.png"))
agreement_tif <- file.path(output_dir, paste0(output_prefix, "_agreement_sd.tif"))
agreement_png <- file.path(output_dir, paste0(output_prefix, "_agreement_sd.png"))
rf_minus_maxent_tif <- file.path(output_dir, paste0(output_prefix, "_rf_minus_maxent.tif"))
rf_minus_maxent_png <- file.path(output_dir, paste0(output_prefix, "_rf_minus_maxent.png"))
rf_minus_glm_tif <- file.path(output_dir, paste0(output_prefix, "_rf_minus_glm.tif"))
rf_minus_glm_png <- file.path(output_dir, paste0(output_prefix, "_rf_minus_glm.png"))
maxent_minus_glm_tif <- file.path(output_dir, paste0(output_prefix, "_maxent_minus_glm.tif"))
maxent_minus_glm_png <- file.path(output_dir, paste0(output_prefix, "_maxent_minus_glm.png"))
diagnostics_csv <- file.path(output_dir, paste0(output_prefix, "_raster_diagnostics.csv"))

cat("Lade Raster...\n")
rf_raster <- rast(rf_tif)[[1]]
maxent_raster <- rast(maxent_tif)[[1]]
glm_raster <- rast(glm_tif)[[1]]

names(rf_raster) <- "rf"
names(maxent_raster) <- "maxent"
names(glm_raster) <- "glm"

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

align_to_rf <- function(r, label) {
  same_crs <- same.crs(r, rf_raster)
  same_grid <- compareGeom(
    r,
    rf_raster,
    stopOnError = FALSE
  )

  cat(label, "- Gleiches CRS:", same_crs, "\n")
  cat(label, "- Gleiches Raster/Ausdehnung/Aufloesung:", same_grid, "\n")

  if (!same_crs) {
    project(
      r,
      rf_raster,
      method = "bilinear"
    )
  } else if (!same_grid) {
    resample(
      r,
      rf_raster,
      method = "bilinear"
    )
  } else {
    r
  }
}

cat("\nRasterdiagnose vor der Angleichung:\n")
diagnostics_before <- rbind(
  raster_summary(rf_raster, "rf", rf_tif),
  raster_summary(maxent_raster, "maxent", maxent_tif),
  raster_summary(glm_raster, "glm", glm_tif)
)
print(diagnostics_before)

cat("\nVergleichbarkeit vor der Angleichung:\n")
maxent_aligned <- align_to_rf(maxent_raster, "MaxEnt")
glm_aligned <- align_to_rf(glm_raster, "GLM")

names(maxent_aligned) <- "maxent"
names(glm_aligned) <- "glm"

diagnostics_after <- rbind(
  raster_summary(rf_raster, "rf_aligned", rf_tif),
  raster_summary(maxent_aligned, "maxent_aligned", maxent_tif),
  raster_summary(glm_aligned, "glm_aligned", glm_tif)
)

prediction_stack <- c(rf_raster, maxent_aligned, glm_aligned)
names(prediction_stack) <- c("rf", "maxent", "glm")

valid_all <- !is.na(rf_raster) & !is.na(maxent_aligned) & !is.na(glm_aligned)
valid_any <- !is.na(rf_raster) | !is.na(maxent_aligned) | !is.na(glm_aligned)

all_valid_cells <- global(valid_all, "sum", na.rm = TRUE)[1, 1]
any_valid_cells <- global(valid_any, "sum", na.rm = TRUE)[1, 1]

cat("\nRasterdiagnose nach der Angleichung:\n")
print(diagnostics_after)
cat("Zellen mit Werten aller Modelle:", all_valid_cells, "\n")
cat("Zellen mit mindestens einem Modellwert:", any_valid_cells, "\n")

write_csv(
  rbind(diagnostics_before, diagnostics_after),
  diagnostics_csv
)

cat("\nBerechne gewichteten Mittelwert...\n")
cat("RF-Gewicht:", rf_weight, "\n")
cat("MaxEnt-Gewicht:", maxent_weight, "\n")
cat("GLM-Gewicht:", glm_weight, "\n")

weighted_mean_fun <- function(rf, maxent, glm) {
  numerator <- ifelse(is.na(rf), 0, rf * rf_weight) +
    ifelse(is.na(maxent), 0, maxent * maxent_weight) +
    ifelse(is.na(glm), 0, glm * glm_weight)

  denominator <- ifelse(is.na(rf), 0, rf_weight) +
    ifelse(is.na(maxent), 0, maxent_weight) +
    ifelse(is.na(glm), 0, glm_weight)

  ifelse(denominator == 0, NA, numerator / denominator)
}

agreement_fun <- function(rf, maxent, glm) {
  values <- cbind(rf, maxent, glm)
  apply(values, 1, sd, na.rm = TRUE)
}

mean_raster <- lapp(
  prediction_stack,
  fun = weighted_mean_fun
)
names(mean_raster) <- "rf_maxent_glm_weighted_mean"

agreement_raster <- lapp(
  prediction_stack,
  fun = agreement_fun
)
names(agreement_raster) <- "agreement_sd"

rf_minus_maxent_raster <- rf_raster - maxent_aligned
names(rf_minus_maxent_raster) <- "rf_minus_maxent"

rf_minus_glm_raster <- rf_raster - glm_aligned
names(rf_minus_glm_raster) <- "rf_minus_glm"

maxent_minus_glm_raster <- maxent_aligned - glm_aligned
names(maxent_minus_glm_raster) <- "maxent_minus_glm"

cat("\nWertebereiche der Ausgaben:\n")
print(global(mean_raster, "range", na.rm = TRUE))
print(global(agreement_raster, "range", na.rm = TRUE))
print(global(rf_minus_maxent_raster, "range", na.rm = TRUE))
print(global(rf_minus_glm_raster, "range", na.rm = TRUE))
print(global(maxent_minus_glm_raster, "range", na.rm = TRUE))

cat("\nSpeichere GeoTIFFs...\n")
writeRaster(mean_raster, mean_tif, overwrite = TRUE)
writeRaster(agreement_raster, agreement_tif, overwrite = TRUE)
writeRaster(rf_minus_maxent_raster, rf_minus_maxent_tif, overwrite = TRUE)
writeRaster(rf_minus_glm_raster, rf_minus_glm_tif, overwrite = TRUE)
writeRaster(maxent_minus_glm_raster, maxent_minus_glm_tif, overwrite = TRUE)

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
        basename(maxent_tif),
        "| GLM:",
        basename(glm_tif)
      ),
      x = "Laengengrad",
      y = "Breitengrad"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 8),
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
  "Gewichtete Mittelwert-Heatmap: Random Forest + MaxEnt + GLM",
  "Gewichteter Mittelwert",
  "magma",
  limits = c(0, 1)
)

plot_raster_png(
  agreement_raster,
  agreement_png,
  "Modell-Unstimmigkeit: Standardabweichung zwischen RF, MaxEnt und GLM",
  "SD",
  "magma",
  limits = c(0, 1)
)

plot_raster_png(
  rf_minus_maxent_raster,
  rf_minus_maxent_png,
  "Differenz-Heatmap: Random Forest - MaxEnt",
  "RF - MaxEnt",
  "viridis"
)

plot_raster_png(
  rf_minus_glm_raster,
  rf_minus_glm_png,
  "Differenz-Heatmap: Random Forest - GLM",
  "RF - GLM",
  "viridis"
)

plot_raster_png(
  maxent_minus_glm_raster,
  maxent_minus_glm_png,
  "Differenz-Heatmap: MaxEnt - GLM",
  "MaxEnt - GLM",
  "viridis"
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
cat(" - ", agreement_tif, "\n", sep = "")
cat(" - ", agreement_png, "\n", sep = "")
cat(" - ", rf_minus_maxent_tif, "\n", sep = "")
cat(" - ", rf_minus_maxent_png, "\n", sep = "")
cat(" - ", rf_minus_glm_tif, "\n", sep = "")
cat(" - ", rf_minus_glm_png, "\n", sep = "")
cat(" - ", maxent_minus_glm_tif, "\n", sep = "")
cat(" - ", maxent_minus_glm_png, "\n", sep = "")
