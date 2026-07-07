# =============================================================================
# rf_test_2.r
#
# Train a Random Forest model using environmental predictors and create
# probability heatmaps for new findings in Bavaria.
#
# Predictors: precipitation, temperature, elevation, DEM, slope
# =============================================================================

library(tidyverse)
library(sf)
library(terra)
library(rnaturalearth)
library(rnaturalearthdata)
library(viridis)
library(ranger)
library(caret)

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

# sample up to 10000 presence points for modeling
#presence_df <- presence_df %>%
#  slice_sample(n = min(10000, nrow(presence_df)))

cat("Presence points sampled:", nrow(presence_df), "\n")

# Load raster data
library(terra)
precip_raster <- rast("data/climate/precipitation.tif")
temp_raster <- rast("data/climate/temperature.tif")
elev_raster <- rast("data/dem/DEU_elv_msk.tif")
dem_raster <- rast("data/dem/dem.tif")
slope_raster <- rast("data/dem/slope.tif")
cat("Rasters loaded\n")

# Extract values at presence points
presence_vect <- vect(as.matrix(select(presence_df, lng, lat)), crs = "EPSG:4326")
presence_df$precipitation <- extract(precip_raster, presence_vect)[, 2]
presence_df$temperature <- extract(temp_raster, presence_vect)[, 2]
presence_df$elevation <- extract(elev_raster, presence_vect)[, 2]
presence_df$dem <- extract(dem_raster, presence_vect)[, 2]
presence_df$slope <- extract(slope_raster, presence_vect)[, 2]
presence_df$presence <- 1

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

# --- 4. Generate absence points and extract features ----------------------
set.seed(42)
presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

n_absence <- nrow(presence_df) * 2
try_pts <- n_absence * 4
absence_pts <- st_sample(bavaria_sf, size = try_pts) %>% st_as_sf() %>% st_set_crs(4326)

# compute min distance from each candidate absence point to any presence
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

# Extract environmental values at absence points
absence_vect <- vect(as.matrix(select(absence_df, lng, lat)), crs = "EPSG:4326")
absence_df$precipitation <- extract(precip_raster, absence_vect)[, 2]
absence_df$temperature <- extract(temp_raster, absence_vect)[, 2]
absence_df$elevation <- extract(elev_raster, absence_vect)[, 2]
absence_df$dem <- extract(dem_raster, absence_vect)[, 2]
absence_df$slope <- extract(slope_raster, absence_vect)[, 2]
absence_df$presence <- 0

cat("Absence points generated:", nrow(absence_df), "\n")

# --- 5. Combine training data -----------------------------------------------
train_df <- bind_rows(
  select(presence_df, lng, lat, precipitation, temperature, elevation, dem, slope, presence),
  select(absence_df, lng, lat, precipitation, temperature, elevation, dem, slope, presence)
) %>%
  drop_na()

cat("Training data:", nrow(train_df), "rows\n")
cat("Presence:", sum(train_df$presence == 1), "| Absence:", sum(train_df$presence == 0), "\n")

# Convert presence to factor
train_df$presence <- factor(train_df$presence, levels = c(0, 1), labels = c("absent", "present"))

# --- 6. Train Random Forest model -------------------------------------------
feature_cols <- c("precipitation", "temperature", "elevation", "dem", "slope")

set.seed(123)
rf_model <- ranger(
  presence ~ .,
  data = select(train_df, all_of(feature_cols), presence),
  num.trees = 500,
  mtry = 3,
  probability = TRUE,
  seed = 123
)

cat("\nRandom Forest model trained\n")
cat("OOB error estimate:", 1 - rf_model$prediction.error, "\n")

# --- 7. Create prediction grid within Bavaria ----------------------------
res_deg <- 0.02  # coarser grid for efficiency
bbox <- st_bbox(bavaria_sf)
xs <- seq(bbox$xmin, bbox$xmax, by = res_deg)
ys <- seq(bbox$ymin, bbox$ymax, by = res_deg)
grid <- expand.grid(lng = xs, lat = ys)
grid_sf <- st_as_sf(grid, coords = c("lng", "lat"), crs = 4326)
in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[, 1]
grid_sf <- grid_sf[in_bavaria, , drop = FALSE]
grid_coords <- st_coordinates(grid_sf) %>% as.data.frame() %>% rename(lng = X, lat = Y)

cat("Prediction grid points:", nrow(grid_coords), "\n")

# --- 8. Extract environmental values at grid points -----------------------
grid_vect <- vect(as.matrix(grid_coords), crs = "EPSG:4326")
grid_coords$precipitation <- extract(precip_raster, grid_vect)[, 2]
grid_coords$temperature <- extract(temp_raster, grid_vect)[, 2]
grid_coords$elevation <- extract(elev_raster, grid_vect)[, 2]
grid_coords$dem <- extract(dem_raster, grid_vect)[, 2]
grid_coords$slope <- extract(slope_raster, grid_vect)[, 2]

# Remove rows with missing values
grid_coords <- grid_coords %>% drop_na()
cat("Grid points after removing NA:", nrow(grid_coords), "\n")

# --- 9. Make predictions ---------------------------------------------------
pred_data <- select(grid_coords, all_of(feature_cols))
pred_prob <- predict(rf_model, data = pred_data)$predictions[, "present"]
grid_coords$prob <- pred_prob

cat("Probability range:", min(pred_prob, na.rm = TRUE), "to", max(pred_prob, na.rm = TRUE), "\n")

# --- 10. Save probability raster as GeoTIFF -------------------------------
# Create a SpatRaster from the prediction grid and write it to disk
prob_grid <- grid_coords %>%
  select(lng, lat, prob)

prob_rast <- rast(xmin = min(prob_grid$lng), xmax = max(prob_grid$lng),
                  ymin = min(prob_grid$lat), ymax = max(prob_grid$lat),
                  nrows = length(unique(prob_grid$lat)),
                  ncols = length(unique(prob_grid$lng)))

# Convert to raster cells and fill values
cell_ids <- cellFromXY(prob_rast, as.matrix(prob_grid[, c("lng", "lat")]))
prob_rast[cell_ids] <- prob_grid$prob
names(prob_rast) <- "probability"

writeRaster(prob_rast, "rf_all_probability.tif", overwrite = TRUE)
cat("Raster saved: rf_all_probability.tif\n")

# --- 11. Create heatmaps with model information ----------------------------
# Prepare info text for the maps
n_presence <- sum(train_df$presence == "present")
n_absence <- sum(train_df$presence == "absent")
n_total <- nrow(train_df)
predictors_text <- paste(feature_cols, collapse = ", ")
oob_accuracy <- round((1 - rf_model$prediction.error) * 100, 2)

info_text <- sprintf(
  "Presence: %d | Absence: %d | Total: %d\nPredictors: %s\nRF OOB Accuracy: %.2f%%",
  n_presence, n_absence, n_total, predictors_text, oob_accuracy
)

# Heatmap without presence points
p_hm_no_pts <- ggplot() +
  geom_raster(data = grid_coords, aes(x = lng, y = lat, fill = prob)) +
  scale_fill_viridis_c(option = "magma", limits = c(0, 1)) +
  geom_sf(data = bavaria_sf, fill = NA, color = "black", linewidth = 0.5) +
  coord_sf(expand = FALSE, crs = 4326,
           xlim = c(bbox$xmin, bbox$xmax), ylim = c(bbox$ymin, bbox$ymax)) +
  labs(title = "Random Forest Probability of Archaeological Findings",
       subtitle = info_text,
       fill = "Probability") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, face = "italic"),
        legend.position = "right")

print(p_hm_no_pts)
ggsave("rf_all_heatmap.png", plot = p_hm_no_pts, width = 12, height = 9, dpi = 300)
cat("Map saved: rf_all_heatmap.png\n")

# Heatmap with presence points
p_hm_wp <- p_hm_no_pts +
  geom_point(data = presence_df, aes(x = lng, y = lat),
             color = "yellow", size = 0.5, alpha = 0.8)

print(p_hm_wp)
ggsave("rf_all_heatmap_wp.png", plot = p_hm_wp, width = 12, height = 9, dpi = 300)
cat("Map saved: rf_all_heatmap_wp.png\n")

# --- End -------------------------------------------------------------------
