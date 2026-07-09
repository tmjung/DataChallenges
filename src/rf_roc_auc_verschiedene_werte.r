suppressPackageStartupMessages({
  library(tidyverse)
  library(sf)
  library(terra)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(ranger)
})

if (file.exists("src/project_paths.r")) {
  source("src/project_paths.r")
} else {
  source("project_paths.r")
}

set.seed(42)

# -------------------------------------------------------------------------
# 1. Daten laden
# -------------------------------------------------------------------------

raw <- suppressMessages(read_csv(
  project_path("data", "ffm_vfpa_eisenzeit.csv"),
  locale = locale(decimal_mark = ","),
  col_types = cols(
    lng_wgs84 = col_character(),
    lat_wgs84 = col_character()
  ),
  show_col_types = FALSE
))

presence_df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84)
  ) %>%
  select(lng, lat) %>%
  drop_na()

# -------------------------------------------------------------------------
# 2. Raster laden und Jahreswerte erzeugen
# -------------------------------------------------------------------------

precip_monthly <- rast(project_path("data", "climate", "precipitation.tif"))
temp_monthly   <- rast(project_path("data", "climate", "temperature.tif"))
elevation      <- rast(project_path("data", "dem", "DEU_elv_msk.tif"))
dem            <- rast(project_path("data", "dem", "dem.tif"))
slope          <- rast(project_path("data", "dem", "slope.tif"))

precipitation <- sum(precip_monthly, na.rm = TRUE)
names(precipitation) <- "precipitation"

temperature <- mean(temp_monthly, na.rm = TRUE)
names(temperature) <- "temperature"

names(elevation) <- "elevation"
names(dem) <- "dem"
names(slope) <- "slope"

# -------------------------------------------------------------------------
# 3. Bayern-Grenze laden
# -------------------------------------------------------------------------

states_raw <- suppressMessages(ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
))

bavaria_sf <- states_raw %>%
  filter(admin == "Germany", name_de == "Bayern") %>%
  st_transform(4326)

# -------------------------------------------------------------------------
# 4. Hilfsfunktionen
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

calc_auc <- function(actual, scores, positive_class = "present") {
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

create_absence_points <- function(n_absence, min_distance_m) {
  presence_buffer <- st_buffer(presence_m, dist = min_distance_m)
  presence_buffer_union <- st_union(presence_buffer)
  absence_area <- st_difference(bavaria_m, presence_buffer_union)

  absence_m <- st_sample(
    absence_area,
    size = n_absence,
    type = "random"
  )

  absence_sf <- st_as_sf(absence_m)
  absence_sf <- st_transform(absence_sf, 4326)

  st_coordinates(absence_sf) %>%
    as.data.frame() %>%
    rename(lng = X, lat = Y)
}

run_experiment <- function(absence_factor, min_distance_m) {
  n_absence <- nrow(presence_df) * absence_factor

  absence_df <- create_absence_points(
    n_absence = n_absence,
    min_distance_m = min_distance_m
  )
  absence_df <- extract_features(absence_df)
  absence_df$presence <- 0

  train_df <- bind_rows(
    presence_df %>% select(lng, lat, all_of(feature_cols), presence),
    absence_df %>% select(lng, lat, all_of(feature_cols), presence)
  ) %>%
    drop_na()

  train_df$presence <- factor(
    train_df$presence,
    levels = c(0, 1),
    labels = c("absent", "present")
  )

  rf_model <- ranger(
    presence ~ .,
    data = train_df %>% select(all_of(feature_cols), presence),
    num.trees = 500,
    mtry = 3,
    probability = TRUE,
    importance = "impurity",
    seed = 123
  )

  oob_probs <- rf_model$predictions
  oob_pred_class <- colnames(oob_probs)[
    max.col(oob_probs, ties.method = "first")
  ]

  oob_accuracy <- mean(oob_pred_class == as.character(train_df$presence))
  oob_auc <- calc_auc(
    actual = train_df$presence,
    scores = oob_probs[, "present"],
    positive_class = "present"
  )

  cat(
    "Presence:", sum(train_df$presence == "present"), "|",
    "Absence:", sum(train_df$presence == "absent"), "|",
    "Min distance:", min_distance_m, "m |",
    "OOB Brier score:", rf_model$prediction.error, "|",
    "OOB Accuracy:", oob_accuracy, "|",
    "OOB ROC AUC:", oob_auc, "\n"
  )
}

# -------------------------------------------------------------------------
# 5. Presence vorbereiten und Experimente ausfuehren
# -------------------------------------------------------------------------

presence_df <- extract_features(presence_df)
presence_df$presence <- 1

presence_sf <- st_as_sf(
  presence_df,
  coords = c("lng", "lat"),
  crs = 4326,
  remove = FALSE
)

bavaria_m <- st_transform(bavaria_sf, 25832)
presence_m <- st_transform(presence_sf, 25832)

absence_factors <- 1:3
min_distances_m <- seq(1000, 10000, by = 1000)

for (absence_factor in absence_factors) {
  for (min_distance_m in min_distances_m) {
    set.seed(42 + absence_factor * 100000 + min_distance_m)
    run_experiment(
      absence_factor = absence_factor,
      min_distance_m = min_distance_m
    )
  }
}
