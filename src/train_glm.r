library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

set.seed(42)

# -------------------------------------------------------------------------
# 1. Einstellungen
# -------------------------------------------------------------------------

output_dir <- "output/glm"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

presence_absence_png <- file.path(output_dir, "glm_presence_absence_points.png")
probability_tif <- file.path(output_dir, "glm_probability.tif")
heatmap_png <- file.path(output_dir, "glm_probability_heatmap.png")
heatmap_points_png <- file.path(output_dir, "glm_fundwahrscheinlichkeit_mit_presence_punkten.png")
coefficients_csv <- file.path(output_dir, "glm_coefficients.csv")
evaluation_csv <- file.path(output_dir, "glm_evaluation.csv")

# -------------------------------------------------------------------------
# 2. Daten laden
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
    lng = as.numeric(gsub(",", ".", lng_wgs84)),
    lat = as.numeric(gsub(",", ".", lat_wgs84))
  ) %>%
  select(lng, lat) %>%
  drop_na()

cat("Presence points:", nrow(presence_df), "\n")

# -------------------------------------------------------------------------
# 3. Raster laden und Jahreswerte erzeugen
# -------------------------------------------------------------------------

precip_monthly <- rast("data/climate/precipitation.tif")
temp_monthly   <- rast("data/climate/temperature.tif")
elevation      <- rast("data/dem/DEU_elv_msk.tif")
dem            <- rast("data/dem/dem.tif")
slope          <- rast("data/dem/slope.tif")

precipitation <- sum(precip_monthly, na.rm = TRUE)
names(precipitation) <- "precipitation"

temperature <- mean(temp_monthly, na.rm = TRUE)
names(temperature) <- "temperature"

names(elevation) <- "elevation"
names(dem) <- "dem"
names(slope) <- "slope"

cat("Raster loaded and annual predictors created\n")

# -------------------------------------------------------------------------
# 4. Bayern-Grenze laden
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
# 5. Hilfsfunktionen
# -------------------------------------------------------------------------

feature_cols <- c(
  "precipitation",
  "temperature",
  "elevation",
  "dem",
  "slope"
)

extract_features <- function(df) {
  pts <- vect(
    as.matrix(df[, c("lng", "lat")]),
    crs = "EPSG:4326"
  )

  df$precipitation <- terra::extract(precipitation, pts)[, 2]
  df$temperature   <- terra::extract(temperature, pts)[, 2]
  df$elevation     <- terra::extract(elevation, pts)[, 2]
  df$dem           <- terra::extract(dem, pts)[, 2]
  df$slope         <- terra::extract(slope, pts)[, 2]

  df
}

calc_auc <- function(actual, scores, positive_class = 1) {
  actual_binary <- actual == positive_class
  n_positive <- sum(actual_binary)
  n_negative <- sum(!actual_binary)

  if (n_positive == 0 || n_negative == 0) {
    return(NA_real_)
  }

  ranks <- rank(scores, ties.method = "average")

  (sum(ranks[actual_binary]) - n_positive * (n_positive + 1) / 2) /
    (n_positive * n_negative)
}

# -------------------------------------------------------------------------
# 6. Presence- und Absence-Punkte vorbereiten
# -------------------------------------------------------------------------

presence_df <- extract_features(presence_df)
presence_df$presence <- 1

presence_sf <- st_as_sf(
  presence_df,
  coords = c("lng", "lat"),
  crs = 4326,
  remove = FALSE
)

n_absence <- nrow(presence_df) * 2

bavaria_m  <- st_transform(bavaria_sf, 25832)
presence_m <- st_transform(presence_sf, 25832)

presence_buffer <- st_buffer(presence_m, dist = 5000)
presence_buffer_union <- st_union(presence_buffer)
absence_area <- st_difference(bavaria_m, presence_buffer_union)

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

train_df <- bind_rows(
  presence_df %>% select(lng, lat, all_of(feature_cols), presence),
  absence_df  %>% select(lng, lat, all_of(feature_cols), presence)
) %>%
  drop_na()

cat("Training rows:", nrow(train_df), "\n")
cat("Presence:", sum(train_df$presence == 1), "\n")
cat("Absence:", sum(train_df$presence == 0), "\n")

# -------------------------------------------------------------------------
# 7. Presence-/Absence-Punkte als PNG speichern
# -------------------------------------------------------------------------

points_plot_df <- train_df %>%
  mutate(
    point_type = factor(
      presence,
      levels = c(0, 1),
      labels = c("Absence", "Presence")
    )
  )

p_points <- ggplot() +
  geom_sf(
    data = bavaria_sf,
    fill = "grey95",
    color = "black",
    linewidth = 0.5
  ) +
  geom_point(
    data = filter(points_plot_df, point_type == "Absence"),
    aes(x = lng, y = lat, color = point_type),
    size = 0.5,
    alpha = 0.55
  ) +
  geom_point(
    data = filter(points_plot_df, point_type == "Presence"),
    aes(x = lng, y = lat, color = point_type),
    size = 0.6,
    alpha = 0.85
  ) +
  scale_color_manual(
    values = c(
      "Absence" = "dodgerblue3",
      "Presence" = "red"
    ),
    name = "Punkttyp"
  ) +
  coord_sf(expand = FALSE, crs = 4326) +
  labs(
    title = "Presence- und Absence-Punkte fuer das GLM",
    subtitle = sprintf(
      "Presence: %d | Absence: %d",
      sum(points_plot_df$point_type == "Presence"),
      sum(points_plot_df$point_type == "Absence")
    ),
    x = "Laengengrad",
    y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

ggsave(
  presence_absence_png,
  plot = p_points,
  width = 12,
  height = 9,
  dpi = 300
)

cat("Map saved:", presence_absence_png, "\n")

# -------------------------------------------------------------------------
# 8. GLM trainieren und evaluieren
# -------------------------------------------------------------------------

set.seed(123)

test_idx <- sample(
  seq_len(nrow(train_df)),
  size = floor(0.25 * nrow(train_df))
)

test_df <- train_df[test_idx, ]
model_df <- train_df[-test_idx, ]

predictor_center <- sapply(model_df[, feature_cols], mean, na.rm = TRUE)
predictor_scale <- sapply(model_df[, feature_cols], sd, na.rm = TRUE)
predictor_scale[predictor_scale == 0] <- 1

scale_predictors <- function(df) {
  scaled <- sweep(df[, feature_cols], 2, predictor_center, "-")
  scaled <- sweep(scaled, 2, predictor_scale, "/")
  as.data.frame(scaled)
}

model_train_df <- bind_cols(
  presence = model_df$presence,
  scale_predictors(model_df)
)

glm_model <- glm(
  presence ~ .,
  data = model_train_df,
  family = binomial(link = "logit")
)

cat("\nGLM trained\n")
print(summary(glm_model))

test_x <- scale_predictors(test_df)
test_pred <- predict(
  glm_model,
  newdata = test_x,
  type = "response"
)

test_auc <- calc_auc(
  actual = test_df$presence,
  scores = test_pred,
  positive_class = 1
)

test_pred_class <- ifelse(test_pred >= 0.5, 1, 0)
test_accuracy <- mean(test_pred_class == test_df$presence)

cat("Test Accuracy:", test_accuracy, "\n")
cat("Test ROC AUC:", test_auc, "\n")

write_csv(
  tibble(
    metric = c("test_accuracy", "test_auc"),
    value = c(test_accuracy, test_auc)
  ),
  evaluation_csv
)

write_csv(
  tibble(
    term = names(coef(glm_model)),
    coefficient = as.numeric(coef(glm_model))
  ),
  coefficients_csv
)

# -------------------------------------------------------------------------
# 9. Raster fuer flaechenhafte Vorhersage vorbereiten
# -------------------------------------------------------------------------

cat("\nPreparing raster stack for spatial prediction...\n")

template <- slope

align_to_template <- function(r, template, method = "bilinear") {
  r2 <- crop(r, template)
  resample(r2, template, method = method)
}

precipitation_aligned <- align_to_template(precipitation, template)
temperature_aligned   <- align_to_template(temperature, template)
elevation_aligned     <- align_to_template(elevation, template)
dem_aligned           <- align_to_template(dem, template)
slope_aligned         <- align_to_template(slope, template)

predictor_stack <- c(
  precipitation_aligned,
  temperature_aligned,
  elevation_aligned,
  dem_aligned,
  slope_aligned
)

names(predictor_stack) <- feature_cols

cat("Predictor stack created:\n")
print(predictor_stack)

cat("\nCropping and masking predictor stack to Bavaria...\n")

bavaria_vect <- vect(bavaria_sf)

predictor_stack <- crop(predictor_stack, bavaria_vect)
predictor_stack <- mask(predictor_stack, bavaria_vect)

cat("Predictor stack masked to Bavaria\n")

# -------------------------------------------------------------------------
# 10. Flaechenhafte GLM-Vorhersage
# -------------------------------------------------------------------------

cat("\nPredicting probability raster...\n")

glm_predict_fun <- function(model, data) {
  data <- as.data.frame(data)
  pred <- rep(NA_real_, nrow(data))
  valid <- complete.cases(data[, feature_cols])

  if (any(valid)) {
    scaled <- sweep(data[valid, feature_cols], 2, predictor_center, "-")
    scaled <- sweep(scaled, 2, predictor_scale, "/")
    scaled <- as.data.frame(scaled)

    pred[valid] <- predict(
      model,
      newdata = scaled,
      type = "response"
    )
  }

  pred
}

probability_raster <- terra::predict(
  predictor_stack,
  glm_model,
  fun = glm_predict_fun,
  na.rm = FALSE
)

names(probability_raster) <- "glm_probability"

cat("Prediction finished\n")
cat("Probability range:\n")
print(global(probability_raster, "range", na.rm = TRUE))

# -------------------------------------------------------------------------
# 11. GeoTIFF und Heatmaps speichern
# -------------------------------------------------------------------------

writeRaster(
  probability_raster,
  probability_tif,
  overwrite = TRUE
)

prob_df <- as.data.frame(
  probability_raster,
  xy = TRUE,
  na.rm = FALSE
)

raster_extent <- ext(probability_raster)

p_heatmap <- ggplot(prob_df) +
  geom_raster(
    aes(
      x = x,
      y = y,
      fill = glm_probability
    )
  ) +
  scale_fill_viridis_c(
    option = "magma",
    limits = c(0, 1),
    na.value = "transparent",
    name = "Wahrscheinlichkeit"
  ) +
  geom_sf(
    data = bavaria_sf,
    fill = NA,
    color = "black",
    linewidth = 0.5
  ) +
  coord_sf(
    xlim = c(xmin(raster_extent), xmax(raster_extent)),
    ylim = c(ymin(raster_extent), ymax(raster_extent)),
    expand = FALSE
  ) +
  labs(
    title = "GLM-Wahrscheinlichkeit archaeologischer Fundstellen",
    subtitle = sprintf(
      "Test AUC = %.3f | Test Accuracy = %.3f",
      test_auc,
      test_accuracy
    ),
    x = "Laengengrad",
    y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, face = "italic"),
    legend.position = "right"
  )

ggsave(
  heatmap_png,
  plot = p_heatmap,
  width = 12,
  height = 9,
  dpi = 300
)

p_heatmap_points <- p_heatmap +
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
    title = "GLM-Wahrscheinlichkeit mit Presence-Punkten"
  )

ggsave(
  heatmap_points_png,
  plot = p_heatmap_points,
  width = 12,
  height = 9,
  dpi = 300
)

cat("\nSaved outputs:\n")
cat(" - ", presence_absence_png, "\n", sep = "")
cat(" - ", probability_tif, "\n", sep = "")
cat(" - ", heatmap_png, "\n", sep = "")
cat(" - ", heatmap_points_png, "\n", sep = "")
cat(" - ", coefficients_csv, "\n", sep = "")
cat(" - ", evaluation_csv, "\n", sep = "")
