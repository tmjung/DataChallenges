# =============================================================================
# Train MaxEnt model and produce heatmap with presence points
# Data Challenge SS 2026 – Goethe-Universität Frankfurt
# =============================================================================

# NOTE: This script uses a simple set of predictors derived from coordinates
# (longitude, latitude and simple transformations) so it runs without external
# environmental rasters. For a real MaxEnt model, replace these predictors
# with meaningful environmental rasters (soil, elevation, climate, etc.).

# 0. Load packages -----------------------------------------------------------
# Ensure required packages are installed (attempt CRAN install for missing)
pkgs <- c(
  "tidyverse", "sf", "rnaturalearth", "rnaturalearthdata",
  "terra", "maxnet", "viridis", "ggplot2"
)
missing <- pkgs[!(pkgs %in% installed.packages()[, "Package"])]
if (length(missing) > 0) {
  message("Installing missing packages: ", paste(missing, collapse = ", "))
  install.packages(missing, repos = "https://cloud.r-project.org")
}
invisible(lapply(pkgs, function(p) suppressPackageStartupMessages(require(p, character.only = TRUE))))

# 1. Load CSV data ----------------------------------------------------------
raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  col_types = cols(lng_wgs84 = col_character(), lat_wgs84 = col_character()),
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

presence_df <- presence_df %>%
  slice_sample(n = min(1000, nrow(presence_df)))

cat("Presence points sampled:", nrow(presence_df), "\n")

# 2. Load Bavaria boundary --------------------------------------------------
states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

germany_sf <- states_raw %>% filter(admin == "Germany")
bavaria_sf <- germany_sf %>% filter(name_de == "Bayern")

# Define bounding box for plotting / grid
bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)

# 3. Generate background (absence) points ----------------------------------
set.seed(42)

presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

n_absence <- nrow(presence_df) * 2
candidate_pts <- st_sample(bavaria_sf, size = n_absence * 4) %>% st_as_sf() %>% st_set_crs(4326)

# Keep candidate points at least 5 km away from any presence
abs_dist <- st_distance(candidate_pts, presence_sf)
abs_dist_mat <- matrix(as.numeric(abs_dist), nrow = nrow(candidate_pts), ncol = nrow(presence_sf))
min_dist_m <- apply(abs_dist_mat, 1, min)
selected_rows <- which(min_dist_m > 5000)
selected_rows <- selected_rows[seq_len(min(length(selected_rows), n_absence))]
absence_sf <- candidate_pts[selected_rows, ]
absence_df <- st_coordinates(absence_sf) %>% as.data.frame() %>% rename(lng = X, lat = Y)

cat("Background points generated:", nrow(absence_df), "\n")

# 4. (Optional) Cluster presence points only for visualization ----------------
# Clustering does not affect training — it's only for map display.
presence_sf_m <- st_transform(presence_sf, 3035)
neighbors <- st_is_within_distance(presence_sf_m, dist = 2000)
n <- nrow(presence_sf_m)
parent <- seq_len(n)
find_parent <- function(i) {
  while (parent[i] != i) {
    parent[i] <<- parent[parent[i]]
    i <- parent[i]
  }
  i
}
union_parent <- function(a, b) {
  ra <- find_parent(a)
  rb <- find_parent(b)
  if (ra != rb) parent[rb] <<- ra
}
for (i in seq_len(n)) {
  nbrs <- neighbors[[i]]
  if (length(nbrs) > 1) {
    for (j in nbrs[nbrs > i]) union_parent(i, j)
  }
}
cluster_root <- sapply(seq_len(n), find_parent)
cluster_id <- as.integer(factor(cluster_root))
presence_sf$cluster_id <- cluster_id
presence_sf_m$cluster_id <- cluster_id
presence_df$cluster_id <- cluster_id

cluster_polygons <- presence_sf_m %>%
  group_by(cluster_id) %>%
  summarize(geometry = st_union(geometry), .groups = "drop") %>%
  st_buffer(dist = 2500) %>%
  st_make_valid()
cluster_polygons <- st_transform(cluster_polygons, 4326)

# 5. Prepare predictors (simple coordinate-based predictors) -----------------
# We create predictors from lon/lat so the script runs without extra data.
make_predictors <- function(df) {
  df %>%
    mutate(
      lon = lng,
      lat = lat,
      lon2 = lon^2,
      lat2 = lat^2,
      lon_lat = lon * lat
    ) %>%
    select(lon, lat, lon2, lat2, lon_lat)
}

pres_preds <- make_predictors(presence_df)
abs_preds <- make_predictors(absence_df)

# 6. Prepare training data for maxnet --------------------------------------
train_xy <- bind_rows(
  pres_preds %>% mutate(presence = 1),
  abs_preds %>% mutate(presence = 0)
)

train_data <- train_xy %>% select(-presence)
train_p <- train_xy$presence

# 7. Train MaxEnt-like model using maxnet ----------------------------------
cat("Training maxnet model...\n")
mx <- maxnet::maxnet(p = train_p, data = train_data)
cat("Model trained.\n")

# 8. Create prediction grid and predict ------------------------------------
lon_seq <- seq(bbox_bavaria['xmin'], bbox_bavaria['xmax'], length.out = 300)
lat_seq <- seq(bbox_bavaria['ymin'], bbox_bavaria['ymax'], length.out = 300)
grid <- expand.grid(lng = lon_seq, lat = lat_seq)
grid_preds <- make_predictors(grid %>% rename(lng = lng, lat = lat))

grid$pred <- predict(mx, newdata = grid_preds, type = "logistic")

# 9. Convert grid to raster and mask to Bavaria for plotting --------------
# terra does not export rastFromXYZ on all versions; create raster from xyz
grid_xyz <- grid %>% rename(x = lng, y = lat, z = pred)
grid_rast <- terra::rast(grid_xyz, type = "xyz")
crs(grid_rast) <- "EPSG:4326"

# Mask the raster so values only remain inside the Bavaria polygon
bav_vect <- terra::vect(bavaria_sf)
grid_rast_masked <- terra::mask(grid_rast, bav_vect)

# Convert masked raster back to data.frame for ggplot
grid_df <- as.data.frame(grid_rast_masked, xy = TRUE)
names(grid_df) <- c("x", "y", "suitability")

# 10. Plot heatmap with presence points ------------------------------------
p <- ggplot() +
  geom_raster(data = grid_df, aes(x = x, y = y, fill = suitability), interpolate = TRUE) +
  scale_fill_viridis_c(option = "viridis", na.value = "transparent") +
  geom_sf(data = bavaria_sf, fill = NA, color = "black", linewidth = 0.4) +
  geom_sf(data = cluster_polygons, fill = "orange", color = NA, alpha = 0.18, show.legend = FALSE) +
  geom_point(data = presence_df, aes(x = lng, y = lat), color = "red", size = 1.2, alpha = 0.9) +
  coord_sf(xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
           ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
           expand = FALSE, crs = 4326) +
  labs(title = "MaxEnt-Prediction (Koordinaten als Prädiktoren)", x = "Längengrad", y = "Breitengrad", fill = "Eignung") +
  theme_minimal()

print(p)
ggsave("train_maxent_1_map.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Map saved: train_maxent_1_map.png\n")

# 11. (Optional) Save model to disk
saveRDS(mx, file = "train_maxent_1_model.rds")
cat("Model saved: train_maxent_1_model.rds\n")

# End of script
