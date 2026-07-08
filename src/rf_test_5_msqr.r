library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(ranger)
library(viridis)

set.seed(42)

# -------------------------------------------------------------------------
# 1. Daten laden
# -------------------------------------------------------------------------

raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  col_types = cols(
    lng_wgs84 = col_character(),
    lat_wgs84 = col_character()
  ),
  show_col_types = FALSE
)

presence_df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84)
  ) %>%
  select(lng, lat) %>%
  drop_na()

# sample up to 10000 presence points for modeling
#presence_df <- presence_df %>%
#  slice_sample(n = min(10000, nrow(presence_df)))

cat("Presence points:", nrow(presence_df), "\n")

# -------------------------------------------------------------------------
# 2. Raster laden und Jahreswerte erzeugen
# -------------------------------------------------------------------------

precip_monthly <- rast("data/climate/precipitation.tif")
temp_monthly   <- rast("data/climate/temperature.tif")
elevation      <- rast("data/dem/DEU_elv_msk.tif")
dem            <- rast("data/dem/dem.tif")
slope          <- rast("data/dem/slope.tif")

# WorldClim precipitation = Monatsniederschlag
precipitation <- sum(precip_monthly, na.rm = TRUE)
names(precipitation) <- "precipitation"

# WorldClim temperature = Monatsmitteltemperatur
temperature <- mean(temp_monthly, na.rm = TRUE)
names(temperature) <- "temperature"

names(elevation) <- "elevation"
names(dem) <- "dem"
names(slope) <- "slope"

cat("Raster loaded and annual predictors created\n")

# -------------------------------------------------------------------------
# 3. Bayern-Grenze laden
# -------------------------------------------------------------------------

states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

bavaria_sf <- states_raw %>%
  filter(admin == "Germany", name_de == "Bayern") %>%
  st_transform(4326)

cat("Bavaria boundary loaded\n")

# -------------------------------------------------------------------------
# MSQR / Müncheberger Soil Quality Rating laden
# -------------------------------------------------------------------------

msqr_raw <- read_csv(
  "data/msqr.csv",
  show_col_types = FALSE
)

msqr_sf <- st_as_sf(
  msqr_raw,
  wkt = "WKT",
  crs = 4326
)

cat("MSQR polygons loaded:", nrow(msqr_sf), "\n")

# Aus Klassen wie "60 - <70 (mittel)" Klassenmittelwert extrahieren
msqr_sf <- msqr_sf %>%
  mutate(
    msqr = case_when(
      sqr == "<35 (äußerst gering)" ~ 17.5,
      sqr == "35 - <50 (sehr gering)" ~ 42.5,
      sqr == "50 - <60 (gering)" ~ 55,
      sqr == "60 - <70 (mittel)" ~ 65,
      sqr == "70 - <85 (hoch)" ~ 77.5,
      sqr == ">=85 (sehr hoch)" ~ 92.5,
      TRUE ~ NA_real_
    )
  )

# Optional: nur Bayern behalten
msqr_sf <- st_intersection(
  st_make_valid(msqr_sf),
  st_make_valid(bavaria_sf)
)


msqr_raster <- rasterize(
  vect(msqr_sf),
  precipitation,
  field = "msqr",
  touches = TRUE
)

names(msqr_raster) <- "msqr"

cat("MSQR raster created\n")

plot(msqr_raster)
global(msqr_raster, "notNA")

# -------------------------------------------------------------------------
# MSQR-Raster als PNG speichern
# -------------------------------------------------------------------------

msqr_df <- as.data.frame(
  msqr_raster,
  xy = TRUE,
  na.rm = TRUE
)

p_msqr <- ggplot() +
  geom_raster(
    data = msqr_df,
    aes(x = x, y = y, fill = msqr)
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    na.value = "transparent",
    name = "MSQR"
  ) +
  geom_sf(
    data = bavaria_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  coord_sf(
    expand = FALSE,
    crs = 4326
  ) +
  labs(
    title = "Muencheberger Soil Quality Rating",
    subtitle = "MSQR raster data for Bavaria"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

print(p_msqr)

ggsave(
  "rf_7_msqr_map.png",
  plot = p_msqr,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Map saved: rf_7_msqr_map.png\n")

# -------------------------------------------------------------------------
# 4. Hilfsfunktion: Rasterwerte an Punkten extrahieren
# -------------------------------------------------------------------------

extract_features <- function(df) {
  pts <- vect(
    as.matrix(df[, c("lng", "lat")]),
    crs = "EPSG:4326"
  )

  df$precipitation <- terra::extract(precipitation, pts)[, 2]
  df$temperature   <- terra::extract(temperature, pts)[, 2]
  df$elevation     <- terra::extract(elevation, pts)[, 2]
  df$dem           <- terra::extract(dem, pts)[, 2]
  df$slope <- terra::extract(slope, pts)[, 2]
  df$msqr  <- terra::extract(msqr_raster, pts)[, 2]

  df
}

presence_df <- extract_features(presence_df)
presence_df$presence <- 1

cat("Presence environmental values extracted\n")


feature_cols <- c(
  "precipitation",
  "temperature",
  "elevation",
  "dem",
  "slope",
  "msqr"
)

# -------------------------------------------------------------------------
# 5. Absence-Punkte erzeugen mit Buffer-Lösung
# -------------------------------------------------------------------------

presence_sf <- st_as_sf(
  presence_df,
  coords = c("lng", "lat"),
  crs = 4326,
  remove = FALSE
)

n_absence <- nrow(presence_df) * 2

# Für Distanzberechnung in metrisches Koordinatensystem transformieren
bavaria_m  <- st_transform(bavaria_sf, 25832)
presence_m <- st_transform(presence_sf, 25832)

# 5-km-Puffer um Presence-Punkte
presence_buffer <- st_buffer(presence_m, dist = 2000)

# Buffer zusammenfassen
presence_buffer_union <- st_union(presence_buffer)

# Gültige Absence-Fläche = Bayern minus 5-km-Buffer
absence_area <- st_difference(bavaria_m, presence_buffer_union)

# Absence-Punkte direkt nur in gültiger Fläche samplen
absence_m <- st_sample(
  absence_area,
  size = n_absence,
  type = "random"
)

absence_sf <- st_as_sf(absence_m)
absence_sf <- st_transform(absence_sf, 4326)

absence_df <- st_coordinates(absence_sf) %>%
  as.data.frame() %>%
  rename(lng = X, lat = Y)

absence_df <- extract_features(absence_df)
absence_df$presence <- 0

cat("Absence points selected:", nrow(absence_df), "\n")

# -------------------------------------------------------------------------
# 6. Trainingsdaten vorbereiten
# -------------------------------------------------------------------------

train_df <- bind_rows(
  presence_df %>% select(lng, lat, all_of(feature_cols), presence),
  absence_df  %>% select(lng, lat, all_of(feature_cols), presence)
) %>%
  drop_na()

train_df$presence <- factor(
  train_df$presence,
  levels = c(0, 1),
  labels = c("absent", "present")
)

cat("Training rows:", nrow(train_df), "\n")
cat("Presence:", sum(train_df$presence == "present"), "\n")
cat("Absence:", sum(train_df$presence == "absent"), "\n")

# -------------------------------------------------------------------------
# 7. Presence-/Absence-Punkte als PNG speichern
# -------------------------------------------------------------------------

points_plot_df <- train_df %>%
  select(lng, lat, presence)

p_points <- ggplot() +
  geom_sf(
    data = bavaria_sf,
    fill = "grey95",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    data = filter(points_plot_df, presence == "absent"),
    aes(x = lng, y = lat, color = presence),
    size = 0.5,
    alpha = 0.55
  ) +
  geom_point(
    data = filter(points_plot_df, presence == "present"),
    aes(x = lng, y = lat, color = presence),
    size = 0.6,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "absent" = "dodgerblue3",
      "present" = "red"
    ),
    labels = c(
      "absent" = "Absence",
      "present" = "Presence"
    ),
    name = "Point type"
  ) +
  coord_sf(
    expand = FALSE,
    crs = 4326
  ) +
  labs(
    title = "Presence and Absence Points",
    subtitle = sprintf(
      "Presence: %d | Absence: %d",
      sum(points_plot_df$presence == "present"),
      sum(points_plot_df$presence == "absent")
    )
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

print(p_points)

ggsave(
  "rf_7_msqr_presence_absence_points.png",
  plot = p_points,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Map saved: rf_7_msqr_presence_absence_points.png\n")

# -------------------------------------------------------------------------
# 8. Random Forest trainieren
# -------------------------------------------------------------------------

rf_model <- ranger(
  presence ~ .,
  data = train_df %>% select(all_of(feature_cols), presence),
  num.trees = 500,
  mtry = 3,
  probability = TRUE,
  importance = "impurity",
  seed = 123
)

cat("\nRandom Forest trained\n")
cat("Prediction error:", rf_model$prediction.error, "\n")
cat("Approx. OOB accuracy:", 1 - rf_model$prediction.error, "\n")

print(rf_model$variable.importance)

# -------------------------------------------------------------------------
# 8. Raster für flächige Vorhersage vorbereiten
# -------------------------------------------------------------------------

cat("\nPreparing raster stack for spatial prediction...\n")

# Referenzraster wählen
# Wichtig: Alle Prädiktoren werden auf dieses Raster gebracht.
# Ich nehme hier precipitation, weil es dein funktionierendes WorldClim-Raster ist.
template <- precipitation

# Funktion: Raster an Template angleichen
align_to_template <- function(r, template, method = "bilinear") {
  r2 <- crop(r, template)
  r2 <- resample(r2, template, method = method)
  return(r2)
}

temperature_aligned <- align_to_template(temperature, template)
elevation_aligned   <- align_to_template(elevation, template)
dem_aligned         <- align_to_template(dem, template)
slope_aligned       <- align_to_template(slope, template)

predictor_stack <- c(
  precipitation,
  temperature_aligned,
  elevation_aligned,
  dem_aligned,
  slope_aligned,
  msqr_raster
)

names(predictor_stack) <- feature_cols

cat("Predictor stack created:\n")
print(predictor_stack)

# -------------------------------------------------------------------------
# 9. Auf Bayern zuschneiden und maskieren
# -------------------------------------------------------------------------

cat("\nCropping and masking predictor stack to Bavaria...\n")

bavaria_vect <- vect(bavaria_sf)

predictor_stack <- crop(predictor_stack, bavaria_vect)
predictor_stack <- mask(predictor_stack, bavaria_vect)

cat("Predictor stack masked to Bavaria\n")

# -------------------------------------------------------------------------
# 10. Flächige Random-Forest-Vorhersage
# -------------------------------------------------------------------------

cat("\nPredicting probability raster...\n")

rf_predict_fun <- function(model, data) {
  data <- as.data.frame(data)

  pred <- predict(
    model,
    data = data
  )$predictions[, "present"]

  return(pred)
}

probability_raster <- terra::predict(
  predictor_stack,
  rf_model,
  fun = rf_predict_fun,
  na.rm = TRUE
)

names(probability_raster) <- "probability"

cat("Prediction finished\n")
cat("Probability range:\n")
print(global(probability_raster, "range", na.rm = TRUE))

# -------------------------------------------------------------------------
# 11. GeoTIFF speichern
# -------------------------------------------------------------------------

writeRaster(
  probability_raster,
  "rf_7_msqr_probability.tif",
  overwrite = TRUE
)

cat("Raster saved: rf_7_msqr_probability.tif\n")

# -------------------------------------------------------------------------
# 12. Heatmap als PNG speichern
# -------------------------------------------------------------------------

prob_df <- as.data.frame(
  probability_raster,
  xy = TRUE,
  na.rm = TRUE
)

n_presence <- sum(train_df$presence == "present")
n_absence  <- sum(train_df$presence == "absent")
n_total    <- nrow(train_df)

predictors_text <- paste(feature_cols, collapse = ", ")
oob_accuracy <- round((1 - rf_model$prediction.error) * 100, 2)

info_text <- sprintf(
  "Presence: %d | Absence: %d | Total: %d\nPredictors: %s\nRF OOB Accuracy: %.2f%%",
  n_presence,
  n_absence,
  n_total,
  predictors_text,
  oob_accuracy
)

p_hm_no_pts <- ggplot() +
  geom_raster(
    data = prob_df,
    aes(x = x, y = y, fill = probability)
  ) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1),
    na.value = "transparent"
  ) +
  geom_sf(
    data = bavaria_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  coord_sf(
    expand = FALSE,
    crs = 4326
  ) +
  labs(
    title = "Random Forest Probability of Archaeological Findings",
    subtitle = info_text,
    fill = "Probability"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

print(p_hm_no_pts)

ggsave(
  "rf_7_msqr_heatmap.png",
  plot = p_hm_no_pts,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Map saved: rf_7_msqr_heatmap.png\n")

# Heatmap mit Presence-Punkten
p_hm_wp <- p_hm_no_pts +
  geom_point(
    data = presence_df,
    aes(x = lng, y = lat),
    color = "yellow",
    size = 0.5,
    alpha = 0.8
  )

print(p_hm_wp)

ggsave(
  "rf_7_msqr_heatmap_wp.png",
  plot = p_hm_wp,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Map saved: rf_7_msqr_heatmap_wp.png\n")

# -------------------------------------------------------------------------
# 13. Modell und Trainingsdaten speichern
# -------------------------------------------------------------------------

saveRDS(rf_model, "rf_7_msqr_model.rds")
write_csv(train_df, "rf_7_msqr_training_data.csv")

cat("\nDone.\n")
cat("Created:\n")
cat(" - rf_7_msqr_map.png\n")
cat(" - rf_7_msqr_probability.tif\n")
cat(" - rf_7_msqr_presence_absence_points.png\n")
cat(" - rf_7_msqr_heatmap.png\n")
cat(" - rf_7_msqr_heatmap_wp.png\n")
cat(" - rf_7_msqr_model.rds\n")
cat(" - rf_7_msqr_training_data.csv\n")
