# =============================================================================
# rf_test_1.r
#
# Minimal script to draw two maps:
#  1) Presence (sampled) and Absence points in Bavaria
#  2) Probability-of-new-findings map derived from distance to nearest presence
#
# Logic borrowed from `draw_points.r` and `draw_abs_and_pre.r`.
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)

# --- 1. Load CSV data -------------------------------------------------------
raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  col_types = cols(lng_wgs84 = col_character(), lat_wgs84 = col_character()),
  show_col_types = FALSE
)

cat("Rows loaded:", nrow(raw), "\n")

# --- 2. Extract presence points and load raster data -----------------------
presence_df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84)
  ) %>%
  select(lng, lat) %>%
  drop_na()

# sample up to 2000 presence points for plotting and calculations
presence_df <- presence_df %>%
  slice_sample(n = min(2000, nrow(presence_df)))

cat("Presence points sampled:", nrow(presence_df), "\n")

# Load raster data
library(terra)
precip_raster <- rast("data/climate/precipitation.tif")
temp_raster <- rast("data/climate/temperature.tif")
elev_raster <- rast("data/dem/DEU_elv_msk.tif")
dem_raster <- rast("data/dem/dem.tif")
slope_raster <- rast("data/dem/slope.tif")
cat("Rasters loaded\n")

# Extract values at presence points (convert to terra SpatVector)
presence_vect <- vect(as.matrix(select(presence_df, lng, lat)), 
                       crs = "EPSG:4326")
presence_df$precipitation <- extract(precip_raster, presence_vect)[, 2]
presence_df$temperature <- extract(temp_raster, presence_vect)[, 2]
cat("Environmental values extracted at presence points\n")

# --- 3. Load Bavaria boundary ----------------------------------------------
states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

germany_sf <- states_raw %>% filter(admin == "Germany")
bavaria_sf <- germany_sf %>% filter(name_de == "Bayern")
cat("Bavaria SF loaded\n")

# --- 4. Generate absence points --------------------------------------------
set.seed(42)
presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

n_absence <- nrow(presence_df) * 2
try_pts <- n_absence * 4
absence_pts <- st_sample(bavaria_sf, size = try_pts) %>% st_as_sf() %>% st_set_crs(4326)

# compute min distance (meters) from each candidate absence point to any presence
abs_dist <- st_distance(absence_pts, presence_sf)
abs_dist_mat <- matrix(as.numeric(abs_dist), nrow = nrow(absence_pts), ncol = nrow(presence_sf))
min_dist_m <- apply(abs_dist_mat, 1, min)

selected_rows <- which(min_dist_m > 5000)
selected_rows <- selected_rows[seq_len(min(length(selected_rows), n_absence))]
if (length(selected_rows) < n_absence) {
  warning(sprintf("Could only find %d absence points >=5km away (requested %d)",
                  length(selected_rows), n_absence))
}
absence_sf <- absence_pts[selected_rows, , drop = FALSE]
absence_df <- st_coordinates(absence_sf) %>% as.data.frame() %>% rename(lng = X, lat = Y)

cat("Absence points generated:", nrow(absence_df), "\n")

# --- 5. Cluster presence points (for nicer visualization) ------------------
# use 2km linking and buffer clusters for polygon display (same logic as draw_points)
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
  ra <- find_parent(a); rb <- find_parent(b)
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
cat("Clusters found:", length(unique(cluster_id)), "\n")

# --- 6. Plot presence and absence map -------------------------------------
p1 <- ggplot() +
  geom_sf(data = bavaria_sf, fill = "lightgray", color = "black", linewidth = 0.5) +
  geom_sf(data = cluster_polygons,
          aes(fill = factor(cluster_id)),
          color = "black", linewidth = 0.3, alpha = 0.45, show.legend = FALSE) +
  geom_point(data = absence_df, aes(x = lng, y = lat),
             color = "gray40", fill = "gray80", shape = 21, size = 1.0, stroke = 0.3, alpha = 0.8) +
  coord_sf(expand = FALSE, crs = 4326,
           xlim = c(st_bbox(bavaria_sf)$xmin, st_bbox(bavaria_sf)$xmax),
           ylim = c(st_bbox(bavaria_sf)$ymin, st_bbox(bavaria_sf)$ymax)) +
  labs(title = "Presence Clusters and Absence Points in Bavaria",
       x = "Longitude", y = "Latitude") +
  scale_fill_viridis_d(option = "turbo") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        panel.grid = element_line(color = "lightgray", linetype = "dashed"))

print(p1)
ggsave("rf_presence_absence_map.png", plot = p1, width = 10, height = 8, dpi = 300)
cat("Map saved: rf_presence_absence_map.png\n")

# --- 7. Probability-of-new-findings map -----------------------------------
# create a regular grid within Bavaria and compute min distance to presence
res_deg <- 0.02 # approx ~2km; increase for faster runs
bbox <- st_bbox(bavaria_sf)
xs <- seq(bbox$xmin, bbox$xmax, by = res_deg)
ys <- seq(bbox$ymin, bbox$ymax, by = res_deg)
grid <- expand.grid(lng = xs, lat = ys)
grid_sf <- st_as_sf(grid, coords = c("lng", "lat"), crs = 4326)

in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[,1]
grid_sf <- grid_sf[in_bavaria, , drop = FALSE]
grid_coords <- st_coordinates(grid_sf) %>% as.data.frame() %>% rename(lng = X, lat = Y)

cat("Grid points for probability map:", nrow(grid_sf), "\n")

# compute min distance (meters) from each grid point to presence points
dist_mat <- st_distance(grid_sf, presence_sf)
min_dist_m_grid <- apply(matrix(as.numeric(dist_mat), nrow = nrow(grid_sf)), 1, min)

# convert distance to probability (simple exponential decay)
scale_m <- 10000
prob <- exp(-min_dist_m_grid / scale_m)
grid_df <- bind_cols(grid_coords, prob = prob)

# plot probability raster
p2 <- ggplot() +
  geom_raster(data = grid_df, aes(x = lng, y = lat, fill = prob)) +
  scale_fill_viridis_c(option = "magma", na.value = "transparent") +
  geom_sf(data = bavaria_sf, fill = NA, color = "black", linewidth = 0.5) +
  geom_point(data = presence_df, aes(x = lng, y = lat), color = "yellow", size = 0.6, alpha = 0.9) +
  coord_sf(expand = FALSE, crs = 4326,
           xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  labs(title = "Probability of New Findings (distance-based)",
       fill = "Prob.") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14))

print(p2)
ggsave("rf_probability_map.png", plot = p2, width = 10, height = 8, dpi = 300)
cat("Map saved: rf_probability_map.png\n")

# --- 8. Environmental maps: Use real raster data ----------------------------
# Extract raster values at a regular grid for visualization
res_deg <- 0.05  # coarser grid for efficiency
bbox <- st_bbox(bavaria_sf)
xs <- seq(bbox$xmin, bbox$xmax, by = res_deg)
ys <- seq(bbox$ymin, bbox$ymax, by = res_deg)
grid_env <- expand.grid(lng = xs, lat = ys)
grid_env_sf <- st_as_sf(grid_env, coords = c("lng", "lat"), crs = 4326)
in_bavaria_env <- st_intersects(grid_env_sf, bavaria_sf, sparse = FALSE)[, 1]
grid_env_sf <- grid_env_sf[in_bavaria_env, , drop = FALSE]
grid_env_coords <- st_coordinates(grid_env_sf) %>% as.data.frame() %>% rename(lng = X, lat = Y)

# Extract all raster values
grid_env_vect <- vect(as.matrix(grid_env_coords), crs = "EPSG:4326")
precip_vals <- extract(precip_raster, grid_env_vect)[, 2]
temp_vals <- extract(temp_raster, grid_env_vect)[, 2]
elev_vals <- extract(elev_raster, grid_env_vect)[, 2]
dem_vals <- extract(dem_raster, grid_env_vect)[, 2]
slope_vals <- extract(slope_raster, grid_env_vect)[, 2]

grid_env_coords$precipitation <- precip_vals
grid_env_coords$temperature <- temp_vals
grid_env_coords$elevation <- elev_vals
grid_env_coords$dem <- dem_vals
grid_env_coords$slope <- slope_vals

cat("Environmental grid points:", nrow(grid_env_coords), "\n")

# Function to create and save both versions (without and with points)
create_env_map <- function(grid_data, var_name, var_col, title, fill_label, file_prefix, 
                          presence_pts = NULL, absent_pts = NULL, bavaria_boundary = NULL) {
  if (all(is.na(grid_data[[var_col]]))) {
    cat(sprintf("Raster %s has no valid values; skipping.\n", var_col))
    return()
  }
  
  # Map without points
  p_no_pts <- ggplot() +
    geom_raster(data = grid_data, aes(x = lng, y = lat, fill = .data[[var_col]])) +
    scale_fill_viridis_c(option = "magma", na.value = "transparent") +
    geom_sf(data = bavaria_boundary, fill = NA, color = "black", linewidth = 0.5) +
    coord_sf(expand = FALSE, crs = 4326,
             xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
    labs(title = title, fill = fill_label) +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14))
  
  ggsave(paste0(file_prefix, ".png"), plot = p_no_pts, width = 10, height = 8, dpi = 300)
  cat(sprintf("Map saved: %s.png\n", file_prefix))
  
  # Map with points
  p_wp <- p_no_pts +
    geom_point(data = presence_pts, aes(x = lng, y = lat),
               color = "yellow", size = 0.5, alpha = 0.8)
  
  ggsave(paste0(file_prefix, "_wp.png"), plot = p_wp, width = 10, height = 8, dpi = 300)
  cat(sprintf("Map saved: %s_wp.png\n", file_prefix))
}

# Create all environmental maps
create_env_map(grid_env_coords, "precipitation", "precipitation", 
               "Annual Precipitation (Niederschlag_Jahr)", "mm",
               "rf_rain_map", presence_df, absence_df, bavaria_sf)

create_env_map(grid_env_coords, "temperature", "temperature",
               "Annual Temperature (Sonnenstunden proxy)", "°C",
               "rf_sunhours_map", presence_df, absence_df, bavaria_sf)

create_env_map(grid_env_coords, "elevation", "elevation",
               "Elevation (DEM)", "m",
               "rf_elevation_map", presence_df, absence_df, bavaria_sf)

create_env_map(grid_env_coords, "dem", "dem",
               "Digital Elevation Model", "m",
               "rf_dem_map", presence_df, absence_df, bavaria_sf)

create_env_map(grid_env_coords, "slope", "slope",
               "Slope", "°",
               "rf_slope_map", presence_df, absence_df, bavaria_sf)

# --- End -------------------------------------------------------------------
