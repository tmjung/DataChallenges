# =============================================================================
# MaxEnt with CSV-based environmental variables + IDW + repeated evaluation
# No external .tif raster required
# =============================================================================

library(tidyverse)
library(sf)
library(maxnet)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ggplot2)
library(pROC)
library(terra)

# 0. Settings ---------------------------------------------------------------

n_repeats <- 100
base_seed <- 42
presence_sample_n <- 1000
res_deg <- 0.03

# 1. Load CSV ---------------------------------------------------------------

fix_decimal <- function(x) {
  if (is.character(x)) {
    x <- gsub(",", ".", x)
    x <- gsub("[^0-9.\\-]", "", x)
    as.numeric(x)
  } else {
    x
  }
}

raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  col_types = cols(
    lng_wgs84 = col_character(),
    lat_wgs84 = col_character()
  ),
  show_col_types = FALSE
)

feature_cols <- c(
  "Höhe_SRTM1_puffer50m",
  "Neigung_SRTM1_puffer50m",
  "Hangausrichtung_SRTM1_puffer50m",
  "Loess_1zu500k_puffer50m",
  "Wasser_puffer50m",
  "Viewshed_km2",
  "Umfeldanalyse_km2",
  "Reliefenergie",
  "Frosttage_Jahr",
  "Niederschlag_Jahr",
  "Sonnenstunden_Jahr",
  "Temperatur_Jahr"
)

df <- raw %>%
  mutate(across(all_of(feature_cols), fix_decimal)) %>%
  mutate(
    lng = as.numeric(gsub(",", ".", lng_wgs84)),
    lat = as.numeric(gsub(",", ".", lat_wgs84))
  ) %>%
  select(lng, lat, all_of(feature_cols)) %>%
  drop_na()

cat("Usable rows:", nrow(df), "\n")

# 2. Bavaria boundary -------------------------------------------------------

states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

bavaria_sf <- states_raw %>%
  filter(admin == "Germany", name_de == "Bayern") %>%
  st_transform(4326)

bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)

# 3. Prediction grid, fixed for all repetitions -----------------------------

grid <- expand.grid(
  lng = seq(bbox_bavaria["xmin"], bbox_bavaria["xmax"], by = res_deg),
  lat = seq(bbox_bavaria["ymin"], bbox_bavaria["ymax"], by = res_deg)
)

grid_sf <- st_as_sf(
  grid,
  coords = c("lng", "lat"),
  crs = 4326,
  remove = FALSE
)

inside_bavaria <- lengths(st_intersects(grid_sf, bavaria_sf)) > 0
pred_grid_base <- grid[inside_bavaria, ]

cat("Prediction grid points:", nrow(pred_grid_base), "\n")

# 4. IDW function -----------------------------------------------------------

assign_idw_env <- function(target_coords, source_df, k = 8, power = 2) {
  source_coords <- as.matrix(source_df[, c("lng", "lat")])
  source_env <- source_df[, feature_cols] %>%
    mutate(across(everything(), as.numeric))

  result <- lapply(seq_len(nrow(target_coords)), function(i) {
    d <- sqrt(
      (source_coords[, 1] - target_coords[i, 1])^2 +
        (source_coords[, 2] - target_coords[i, 2])^2
    )

    nearest_ids <- order(d)[1:min(k, length(d))]
    nearest_d <- d[nearest_ids]
    nearest_d[nearest_d == 0] <- 1e-10

    weights <- 1 / (nearest_d ^ power)
    weights <- weights / sum(weights)

    env_values <- as.matrix(source_env[nearest_ids, , drop = FALSE])
    weighted <- sweep(env_values, 1, weights, `*`)

    colSums(weighted, na.rm = TRUE)
  })

  result_df <- as.data.frame(do.call(rbind, result))
  colnames(result_df) <- feature_cols
  result_df
}

# 5. Repeated MaxEnt modelling ---------------------------------------------

auc_values <- numeric(n_repeats)
prediction_list <- vector("list", n_repeats)
roc_list <- vector("list", n_repeats)

for (i in seq_len(n_repeats)) {
  cat("\n--- Run", i, "of", n_repeats, "---\n")

  set.seed(base_seed + i)

  # Presence sample
  presence_df <- df %>%
    slice_sample(n = min(presence_sample_n, nrow(df))) %>%
    mutate(presence = 1)

  presence_sf <- st_as_sf(
    presence_df,
    coords = c("lng", "lat"),
    crs = 4326,
    remove = FALSE
  )

  # Background points
  n_background <- nrow(presence_df) * 2

  candidate_bg <- st_sample(
    bavaria_sf,
    size = n_background * 5
  ) %>%
    st_as_sf() %>%
    st_set_crs(4326)

  dist_mat <- st_distance(candidate_bg, presence_sf)

  dist_mat_num <- matrix(
    as.numeric(dist_mat),
    nrow = nrow(candidate_bg),
    ncol = nrow(presence_sf)
  )

  min_dist <- apply(dist_mat_num, 1, min)

  selected <- which(min_dist > 5000)
  selected <- selected[seq_len(min(length(selected), n_background))]

  background_sf <- candidate_bg[selected, ]

  background_df <- st_coordinates(background_sf) %>%
    as.data.frame() %>%
    rename(lng = X, lat = Y)

  background_env <- assign_idw_env(
    target_coords = as.matrix(background_df[, c("lng", "lat")]),
    source_df = presence_df,
    k = 8,
    power = 2
  )

  background_df <- bind_cols(background_df, background_env) %>%
    mutate(presence = 0)

  # Train/test data
  train_all <- bind_rows(
    presence_df %>% select(presence, all_of(feature_cols)),
    background_df %>% select(presence, all_of(feature_cols))
  ) %>%
    drop_na()

  set.seed(base_seed + i + 1000)

  test_idx <- sample(
    seq_len(nrow(train_all)),
    size = floor(0.25 * nrow(train_all))
  )

  test_df <- train_all[test_idx, ]
  train_df <- train_all[-test_idx, ]

  train_p <- train_df$presence
  train_x <- train_df %>%
    select(-presence) %>%
    as.data.frame()

  # Train MaxEnt
  mx_formula <- maxnet.formula(
    p = train_p,
    data = train_x,
    classes = "lqph"
  )

  mx <- maxnet(
    p = train_p,
    data = train_x,
    f = mx_formula,
    regmult = 1.5
  )

  # AUC / ROC
  test_x <- test_df %>%
    select(-presence) %>%
    as.data.frame()

  test_pred <- predict(
    mx,
    newdata = test_x,
    type = "logistic"
  )

  roc_obj <- roc(
    response = test_df$presence,
    predictor = test_pred,
    levels = c(0, 1),
    direction = "<",
    quiet = TRUE
  )

  auc_values[i] <- as.numeric(auc(roc_obj))
  roc_list[[i]] <- roc_obj

  cat("AUC:", round(auc_values[i], 4), "\n")

  # Prediction grid with IDW
  pred_grid <- pred_grid_base

  grid_env <- assign_idw_env(
    target_coords = as.matrix(pred_grid[, c("lng", "lat")]),
    source_df = presence_df,
    k = 8,
    power = 2
  )

  pred_grid <- bind_cols(pred_grid, grid_env)

  pred_x <- pred_grid %>%
    select(all_of(feature_cols)) %>%
    as.data.frame()

  pred_grid$suitability <- predict(
    mx,
    newdata = pred_x,
    type = "logistic"
  )

  prediction_list[[i]] <- pred_grid$suitability
}

# 6. Summarize AUC ----------------------------------------------------------

mean_auc <- mean(auc_values, na.rm = TRUE)
sd_auc <- sd(auc_values, na.rm = TRUE)

cat("\n==============================\n")
cat("Mean AUC:", round(mean_auc, 4), "\n")
cat("SD AUC:", round(sd_auc, 4), "\n")
cat("==============================\n")

auc_summary <- tibble(
  run = seq_len(n_repeats),
  auc = auc_values
)

write_csv(auc_summary, "maxent_repeated_auc_values.csv")

# 7. Average prediction map -------------------------------------------------

prediction_matrix <- do.call(cbind, prediction_list)

pred_grid_final <- pred_grid_base %>%
  mutate(
    suitability_mean = rowMeans(prediction_matrix, na.rm = TRUE),
    suitability_sd = apply(prediction_matrix, 1, sd, na.rm = TRUE)
  )

# 8. Save averaged prediction rasters ---------------------------------------

mean_raster <- rast(
  xmin = bbox_bavaria["xmin"],
  xmax = bbox_bavaria["xmax"],
  ymin = bbox_bavaria["ymin"],
  ymax = bbox_bavaria["ymax"],
  resolution = res_deg,
  crs = "EPSG:4326"
)

mean_cell_idx <- cellFromXY(
  mean_raster,
  as.matrix(pred_grid_final[, c("lng", "lat")])
)

mean_raster[mean_cell_idx] <- pred_grid_final$suitability_mean
names(mean_raster) <- "maxent_suitability_mean"

sd_raster <- rast(mean_raster)
sd_raster[mean_cell_idx] <- pred_grid_final$suitability_sd
names(sd_raster) <- "maxent_suitability_sd"

writeRaster(
  mean_raster,
  "maxent_repeated_mean_prediction.tif",
  overwrite = TRUE
)

writeRaster(
  sd_raster,
  "maxent_repeated_uncertainty.tif",
  overwrite = TRUE
)

# 9. Plot AUC distribution --------------------------------------------------

p_auc <- ggplot(auc_summary, aes(x = auc)) +
  geom_histogram(bins = 20, color = "black", fill = "gray80") +
  geom_vline(xintercept = mean_auc, linetype = "dashed", linewidth = 1) +
  labs(
    title = "AUC-Verteilung über 100 MaxEnt-Läufe",
    subtitle = paste(
      "Mittlerer AUC =",
      round(mean_auc, 3),
      "| SD =",
      round(sd_auc, 3)
    ),
    x = "AUC",
    y = "Anzahl"
  ) +
  theme_minimal()

print(p_auc)

ggsave(
  "maxent_repeated_auc_distribution.png",
  plot = p_auc,
  width = 7,
  height = 6,
  dpi = 300
)

# 10. Plot averaged prediction map -----------------------------------------

p_map <- ggplot() +
  geom_raster(
    data = pred_grid_final,
    aes(x = lng, y = lat, fill = suitability_mean)
  ) +
  geom_sf(
    data = bavaria_sf,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    option = "viridis",
    name = "Standorteignung"
  ) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE
  ) +
  labs(
    title = "Gemittelte MaxEnt-Prediction mit CSV-Umweltvariablen",
    subtitle = paste(
      "100 Läufe | IDW-Interpolation | mittlerer AUC =",
      round(mean_auc, 3),
      "±",
      round(sd_auc, 3)
    ),
    x = "Längengrad",
    y = "Breitengrad"
  ) +
  theme_minimal()

print(p_map)

ggsave(
  "maxent_repeated_mean_prediction.png",
  plot = p_map,
  width = 10,
  height = 8,
  dpi = 300
)

# 11. Plot uncertainty map --------------------------------------------------

p_uncertainty <- ggplot() +
  geom_raster(
    data = pred_grid_final,
    aes(x = lng, y = lat, fill = suitability_sd)
  ) +
  geom_sf(
    data = bavaria_sf,
    fill = NA,
    color = "black",
    linewidth = 0.4
  ) +
  scale_fill_viridis_c(
    option = "magma",
    name = "SD"
  ) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE
  ) +
  labs(
    title = "Unsicherheit der MaxEnt-Prediction",
    subtitle = "Standardabweichung der Standorteignung über 100 Läufe",
    x = "Längengrad",
    y = "Breitengrad"
  ) +
  theme_minimal()

print(p_uncertainty)

ggsave(
  "maxent_repeated_uncertainty.png",
  plot = p_uncertainty,
  width = 10,
  height = 8,
  dpi = 300
)

# 12. Save final prediction table ------------------------------------------

write_csv(pred_grid_final, "maxent_repeated_mean_prediction_grid.csv")

cat("Saved: maxent_repeated_auc_values.csv\n")
cat("Saved: maxent_repeated_auc_distribution.png\n")
cat("Saved: maxent_repeated_mean_prediction.png\n")
cat("Saved: maxent_repeated_mean_prediction.tif\n")
cat("Saved: maxent_repeated_uncertainty.png\n")
cat("Saved: maxent_repeated_uncertainty.tif\n")
cat("Saved: maxent_repeated_mean_prediction_grid.csv\n")
