# =============================================================================
# Train Random Forest and Draw Heatmap for Bavaria Presence Prediction
# Data Challenge SS 2026 – Goethe-Universität Frankfurt
# =============================================================================

# ── 0. Load packages ──────────────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(terra)
library(ranger)
library(caret)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# ── 1. Helper functions ───────────────────────────────────────────────────────
fix_decimal <- function(x) {
  if (is.character(x)) as.numeric(gsub(",", ".", x))
  else x
}

# ── 2. Load CSV data ──────────────────────────────────────────────────────────
raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE
)

cat("Rows loaded:", nrow(raw), "\n")

# ── 3. Feature definitions ────────────────────────────────────────────────────
base_feature_cols <- c(
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
spatial_feature_cols <- c("dist_nearest_presence_km", "presence_count_10km")
model_feature_cols <- c(base_feature_cols, spatial_feature_cols)

# ── 4. Prepare presence observations ─────────────────────────────────────────

df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(across(all_of(base_feature_cols), fix_decimal)) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84),
    presence = 1L
  )

presence_df <- df %>%
  select(lng, lat, all_of(base_feature_cols), presence) %>%
  drop_na()

presence_df <- presence_df %>%
  slice_sample(n = min(1000, nrow(presence_df)))

cat("Presence points sampled:", nrow(presence_df), "\n")

# ── 5. Load Bavaria boundary ──────────────────────────────────────────────────
states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

germany_sf <- states_raw %>%
  filter(admin == "Germany")

bavaria_sf <- germany_sf %>%
  filter(name_de == "Bayern")

cat("Bavaria SF loaded\n")

# Diagnostic: show bounding box and ensure longlat CRS
cat("Bavaria bbox:", paste(names(st_bbox(bavaria_sf)), unname(st_bbox(bavaria_sf))), "\n")
cat("bavaria is longlat:", st_is_longlat(bavaria_sf), "\n")
if (st_crs(bavaria_sf)$epsg != 4326) {
  bavaria_sf <- st_transform(bavaria_sf, 4326)
}

bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)

# ── 6. Generate absence points ────────────────────────────────────────────────
set.seed(42)

presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

n_absence <- nrow(presence_df) * 2
absence_pts <- st_sample(bavaria_sf, size = n_absence * 4) %>%
  st_as_sf() %>%
  st_set_crs(4326)

abs_dist <- st_distance(absence_pts, presence_sf)
abs_dist_mat <- matrix(as.numeric(abs_dist), nrow = nrow(absence_pts), ncol = nrow(presence_sf))
min_dist_m <- apply(abs_dist_mat, 1, min)
selected_rows <- which(min_dist_m > 5000)
selected_rows <- selected_rows[seq_len(min(length(selected_rows), n_absence))]
absence_sf <- absence_pts[selected_rows, ]
absence_df <- st_coordinates(absence_sf) %>%
  as.data.frame() %>%
  rename(lng = X, lat = Y)

cat("Absence points generated:", nrow(absence_df), "\n")

# ── 7. Compute spatial features for presence points ───────────────────────────
pres_dist_matrix <- st_distance(presence_sf, presence_sf)
pres_dist_matrix <- matrix(as.numeric(pres_dist_matrix), nrow = nrow(presence_sf))
diag(pres_dist_matrix) <- Inf

presence_df <- presence_df %>%
  mutate(
    dist_nearest_presence_km = apply(pres_dist_matrix, 1, min) / 1000,
    presence_count_10km = apply(pres_dist_matrix, 1, function(d) sum(d <= 10000))
  )

# ── 8. Derive absence features from presence neighbors ────────────────────────
pres_feat_matrix <- as.matrix(select(presence_df, all_of(base_feature_cols)))
pres_feat_sd <- apply(pres_feat_matrix, 2, sd, na.rm = TRUE)

nearest_full_ids <- t(apply(abs_dist_mat[selected_rows, , drop = FALSE], 1, function(d) order(d)[1:10]))

absence_values <- t(apply(nearest_full_ids, 1, function(ord_row) {
  k_use <- min(3, length(ord_row))
  pick_ids <- sample(ord_row, k_use)
  mu <- colMeans(pres_feat_matrix[pick_ids, , drop = FALSE], na.rm = TRUE)
  noise <- rnorm(length(mu), mean = 0, sd = pmax(1e-6, pres_feat_sd * 0.05))
  mu + noise
}))

absence_values <- as.data.frame(absence_values)
colnames(absence_values) <- base_feature_cols

absence_dist_matrix <- matrix(as.numeric(abs_dist_mat[selected_rows, , drop = FALSE]),
                               nrow = length(selected_rows))
absence_dist_km <- apply(absence_dist_matrix, 1, min) / 1000
absence_count_10km <- apply(absence_dist_matrix, 1, function(d) sum(d <= 10000))

absence_features <- bind_cols(absence_df, absence_values) %>%
  mutate(
    dist_nearest_presence_km = absence_dist_km,
    presence_count_10km = absence_count_10km,
    presence = 0L
  )

# ── 9. Cluster presence points by buffering and union (robust to chaining) ──
# Project to metric CRS, buffer each point, union overlapping buffers, and
# use resulting polygons as clusters.
presence_sf_m <- st_transform(presence_sf, 3035)

# Buffer points by 2 km and union overlapping buffers
buffers <- st_buffer(presence_sf_m, dist = 2000)
unioned <- st_union(buffers)

# Ensure we have individual polygon geometries
polys <- tryCatch(st_cast(unioned, "POLYGON"), error = function(e) unioned)
cluster_polygons_m <- st_as_sf(data.frame(geometry = polys))
cluster_polygons_m$cluster_id <- seq_len(nrow(cluster_polygons_m))

# Assign cluster id to each presence point by spatial intersection
ints <- st_intersects(presence_sf_m, cluster_polygons_m, sparse = FALSE)
cluster_idx <- apply(ints, 1, function(row) {
  ids <- which(row)
  if (length(ids) == 0) NA_integer_ else ids[1]
})
presence_sf$cluster_id <- cluster_idx
presence_sf_m$cluster_id <- cluster_idx
presence_df$cluster_id <- cluster_idx

# Expand polygons slightly for nicer visualization, validate and reproject
cluster_polygons <- st_buffer(cluster_polygons_m, dist = 2500) %>% st_make_valid()
cluster_polygons <- st_transform(cluster_polygons, 4326)

# Use complete-linkage hierarchical clustering (no single-linkage chaining)
# Compute pairwise distances (meters) and cut tree at 5000 m
dist_mat <- st_distance(presence_sf_m, presence_sf_m)
dist_num <- matrix(as.numeric(dist_mat), nrow = nrow(presence_sf_m))
# sanitize distance matrix: set diagonal zero and replace NA/Inf with large finite value
diag(dist_num) <- 0
if (any(!is.finite(dist_num))) {
  max_finite <- max(dist_num[is.finite(dist_num)], na.rm = TRUE)
  if (!is.finite(max_finite) || is.na(max_finite)) max_finite <- 1e6
  dist_num[!is.finite(dist_num)] <- max_finite * 10
}
# Convert to dist object
dist_obj <- as.dist(dist_num)

# Hierarchical clustering with complete linkage to avoid chaining
h <- hclust(dist_obj, method = "complete")
cut_h <- 5000
cluster_id <- cutree(h, h = cut_h)
presence_sf$cluster_id <- cluster_id
presence_sf_m$cluster_id <- cluster_id
presence_df$cluster_id <- cluster_id

# Build cluster polygons using convex hull of member points, then buffer
cluster_polygons_m <- presence_sf_m %>%
  group_by(cluster_id) %>%
  summarize(geometry = st_combine(geometry), .groups = "drop") %>%
  mutate(geometry = st_cast(st_convex_hull(geometry), "POLYGON")) %>%
  st_buffer(dist = 2500) %>%
  st_make_valid()
cluster_polygons <- st_transform(cluster_polygons_m, 4326)

cat("Clusters found:", nrow(cluster_polygons), "\\n")

# ── 10. Plot presence cluster areas and absence points ────────────────────────
p <- ggplot() +
  geom_sf(data = bavaria_sf, fill = "lightgray", color = "black", linewidth = 0.5) +
  geom_sf(data = cluster_polygons,
          aes(fill = factor(cluster_id)),
          color = "black",
          linewidth = 0.3,
          alpha = 0.45,
          show.legend = FALSE) +
  geom_point(data = absence_df,
             aes(x = lng, y = lat),
             color = "gray40",
             fill = "gray80",
             shape = 21,
             size = 1.0,
             stroke = 0.3,
             alpha = 0.8) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE,
    crs = 4326
  ) +
  labs(
    title = "Präsenz-Cluster und Absenzpunkte in Bayern",
    x = "Längengrad",
    y = "Breitengrad",
    fill = "Cluster"
  ) +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    panel.grid = element_line(color = "lightgray", linetype = "dashed")
  )

print(p)
ggsave("train_1_cluster_absence_map.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Saved cluster map: train_1_cluster_absence_map.png\n")

# ── 11. Train random forest model ─────────────────────────────────────────────
train_df <- bind_rows(
  select(presence_df, all_of(model_feature_cols), presence),
  select(absence_features, all_of(model_feature_cols), presence)
) %>%
  mutate(presence = factor(presence, levels = c(0, 1), labels = c("absent", "present"))) %>%
  drop_na()

cat("Training rows:", nrow(train_df), "\n")

set.seed(123)
ctrl <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_model <- train(
  presence ~ .,
  data = train_df,
  method = "ranger",
  trControl = ctrl,
  metric = "ROC",
  tuneGrid = expand.grid(
    mtry = c(3, 5, 7),
    splitrule = "gini",
    min.node.size = c(5, 20)
  ),
  importance = "impurity"
)

print(rf_model)
cat("Best tuning:\n")
print(rf_model$bestTune)
cat("Variable importance:\n")
print(varImp(rf_model))

# ── 12. Build prediction grid and compute features ───────────────────────────
res <- 0.01
pred_grid <- expand.grid(
  lng = seq(bbox_bavaria["xmin"], bbox_bavaria["xmax"], by = res),
  lat = seq(bbox_bavaria["ymin"], bbox_bavaria["ymax"], by = res)
)

grid_sf <- st_as_sf(pred_grid, coords = c("lng", "lat"), crs = 4326)
in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[, 1]
if (!is.logical(in_bavaria)) {
  in_bavaria <- as.logical(in_bavaria)
}
pred_grid <- pred_grid[in_bavaria, ]
grid_sf <- grid_sf[in_bavaria, ]

pres_coords <- as.matrix(select(presence_df, lng, lat))
pres_features <- select(presence_df, all_of(base_feature_cols))

compute_nearest_feature_values <- function(grid_xy, pres_xy, pres_feats, k = 3, chunk_size = 1000) {
  n_grid <- nrow(grid_xy)
  out <- matrix(NA_real_, nrow = n_grid, ncol = ncol(pres_feats))
  colnames(out) <- colnames(pres_feats)

  for (start in seq(1, n_grid, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_grid)
    chunk_xy <- grid_xy[start:end, , drop = FALSE]
    dx <- outer(chunk_xy[, 1], pres_xy[, 1], "-")
    dy <- outer(chunk_xy[, 2], pres_xy[, 2], "-")
    dist2 <- dx^2 + dy^2
    nearest_ids <- t(apply(dist2, 1, function(row) order(row)[1:k]))
    out[start:end, ] <- t(apply(nearest_ids, 1, function(ids) {
      colMeans(pres_feats[ids, , drop = FALSE], na.rm = TRUE)
    }))
  }

  as.data.frame(out)
}

compute_spatial_features <- function(grid_sf, pres_sf, radius_m = 10000, chunk_size = 1000) {
  n_grid <- nrow(grid_sf)
  dist_nearest_km <- numeric(n_grid)
  count_within_radius <- integer(n_grid)

  for (start in seq(1, n_grid, by = chunk_size)) {
    end <- min(start + chunk_size - 1, n_grid)
    dmat <- st_distance(grid_sf[start:end, ], pres_sf)
    dnum <- matrix(as.numeric(dmat), nrow = nrow(dmat), ncol = ncol(dmat))
    dist_nearest_km[start:end] <- apply(dnum, 1, min) / 1000
    count_within_radius[start:end] <- apply(dnum, 1, function(d) sum(d <= radius_m))
  }

  tibble(
    dist_nearest_presence_km = dist_nearest_km,
    presence_count_10km = count_within_radius
  )
}

pred_grid <- bind_cols(
  pred_grid,
  compute_nearest_feature_values(
    as.matrix(select(pred_grid, lng, lat)),
    pres_coords,
    pres_features,
    k = 3,
    chunk_size = 1000
  ),
  compute_spatial_features(grid_sf, presence_sf, radius_m = 10000, chunk_size = 1000)
)

pred_prob <- predict(rf_model,
                     newdata = select(pred_grid, all_of(model_feature_cols)),
                     type = "prob")[, "present"]
pred_grid$prob_present <- pred_prob

# ── 13. Create heatmap raster and save ───────────────────────────────────────
r <- rast(
  xmin = bbox_bavaria["xmin"], xmax = bbox_bavaria["xmax"],
  ymin = bbox_bavaria["ymin"], ymax = bbox_bavaria["ymax"],
  resolution = res,
  crs = "EPSG:4326"
)

cell_idx <- cellFromXY(r, as.matrix(pred_grid[, c("lng", "lat")]))
r[cell_idx] <- pred_grid$prob_present
names(r) <- "fund_wahrscheinlichkeit"

writeRaster(r, "train_1_heatmap.tif", overwrite = TRUE)
cat("Heatmap GeoTIFF saved: train_1_heatmap.tif\n")

heatmap_df <- as.data.frame(r, xy = TRUE) %>%
  rename(lng = x, lat = y, prob = fund_wahrscheinlichkeit) %>%
  filter(!is.na(prob))

bavaria_outline <- germany_sf %>%
  filter(name_de == "Bayern")

plot_subtitle <- sprintf(
  "Model: Random Forest | Präsenz: %d | Absenz: %d | Auflösung: %.3f°",
  sum(train_df$presence == "present"),
  sum(train_df$presence == "absent"),
  res
)

p_heatmap <- ggplot() +
  geom_raster(data = heatmap_df, aes(x = lng, y = lat, fill = prob)) +
  geom_sf(data = bavaria_outline, fill = NA, color = "white", linewidth = 0.4) +
  geom_point(data = presence_df, aes(x = lng, y = lat),
             color = "yellow", size = 0.8, alpha = 0.6, shape = 16) +
  scale_fill_viridis_c(
    option = "magma",
    name = "Auffindungs-\nwahrscheinlichkeit",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent
  ) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE
  ) +
  labs(
    title = "Predictive Modelling – Eisenzeitliche Fundstellen Bayern",
    subtitle = plot_subtitle,
    caption = paste0("Datenbasis: VFPA Eisenzeit | Data Challenge SS2026 | Koo & Jungbeck",
                     " | Rasterzellen: ", nrow(pred_grid)),
    x = "Längengrad",
    y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right")

print(p_heatmap)
ggsave("train_1_heatmap.png", plot = p_heatmap, width = 12, height = 9, dpi = 300)
cat("Heatmap PNG saved: train_1_heatmap.png\n")
