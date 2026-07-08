# =============================================================================
# Mean heatmap from Random Forest and MaxEnt GeoTIFF predictions
# =============================================================================

library(terra)
library(ggplot2)
library(viridis)

# Usage:
# Rscript src/rf_me_mean_heatmap.r [maxent_tif] [rf_tif] [output_prefix]

args <- commandArgs(trailingOnly = TRUE)

maxent_tif <- ifelse(
  length(args) >= 1,
  args[1],
  "maxent_repeated_mean_prediction.tif"
)

rf_tif <- ifelse(
  length(args) >= 2,
  args[2],
  "rf_final_5km_blue_probability.tif"
)

output_prefix <- ifelse(
  length(args) >= 3,
  args[3],
  "rf_maxent_mean_heatmap"
)

output_tif <- paste0(output_prefix, ".tif")
output_png <- paste0(output_prefix, ".png")

if (!file.exists(maxent_tif)) {
  stop("MaxEnt GeoTIFF not found: ", maxent_tif)
}

if (!file.exists(rf_tif)) {
  stop("Random Forest GeoTIFF not found: ", rf_tif)
}

cat("Loading rasters...\n")
maxent_raster <- rast(maxent_tif)
rf_raster <- rast(rf_tif)

maxent_raster <- maxent_raster[[1]]
rf_raster <- rf_raster[[1]]

names(maxent_raster) <- "maxent"
names(rf_raster) <- "rf"

cat("RF raster:\n")
print(rf_raster)

cat("MaxEnt raster:\n")
print(maxent_raster)

# Use the RF raster as the target grid, because the requested RF file is the
# final 5 km probability raster. MaxEnt is aligned to this grid before averaging.
cat("Aligning MaxEnt raster to RF raster grid...\n")

same_crs <- same.crs(maxent_raster, rf_raster)
same_grid <- compareGeom(
  maxent_raster,
  rf_raster,
  stopOnError = FALSE
)

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

cat("Calculating cell-wise mean...\n")
prediction_stack <- c(rf_raster, maxent_aligned)
names(prediction_stack) <- c("rf", "maxent")

mean_raster <- mean(prediction_stack, na.rm = FALSE)

names(mean_raster) <- "rf_maxent_mean_probability"

cat("Mean raster range:\n")
print(global(mean_raster, "range", na.rm = TRUE))

cat("Saving GeoTIFF...\n")
writeRaster(
  mean_raster,
  output_tif,
  overwrite = TRUE
)

cat("Saving PNG heatmap...\n")
mean_df <- as.data.frame(
  mean_raster,
  xy = TRUE,
  na.rm = TRUE
)

p_mean <- ggplot(mean_df) +
  geom_raster(
    aes(
      x = x,
      y = y,
      fill = rf_maxent_mean_probability
    )
  ) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1),
    na.value = "transparent",
    name = "Mean probability"
  ) +
  coord_equal(expand = FALSE) +
  labs(
    title = "Mean Heatmap: Random Forest + MaxEnt",
    subtitle = paste(
      "RF:",
      basename(rf_tif),
      "| MaxEnt:",
      basename(maxent_tif)
    ),
    x = "Longitude",
    y = "Latitude"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 9),
    legend.position = "right"
  )

ggsave(
  output_png,
  plot = p_mean,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Saved: ", output_tif, "\n", sep = "")
cat("Saved: ", output_png, "\n", sep = "")
