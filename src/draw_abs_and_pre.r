# =============================================================================
# Draw Map of Bavaria with Presence Clusters and Absence Points
# Data Challenge SS 2026 – Goethe-Universität Frankfurt
# =============================================================================

# ── 0. Load packages ──────────────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)

# ── 1. Load CSV data ──────────────────────────────────────────────────────────
raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  col_types = cols(lng_wgs84 = col_character(), lat_wgs84 = col_character()),
  show_col_types = FALSE
)

cat("Rows loaded:", nrow(raw), "\n")

# ── 2. Extract presence points ────────────────────────────────────────────────
presence_df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84)
  ) %>%
  select(lng, lat) %>%
  drop_na()

presence_df <- presence_df %>%
  slice_sample(n = min(2000, nrow(presence_df)))

cat("Presence points sampled:", nrow(presence_df), "\n")
cat("Coordinate ranges:\n")
cat("  lng: ", min(presence_df$lng, na.rm = TRUE), "to", max(presence_df$lng, na.rm = TRUE), "\n")
cat("  lat: ", min(presence_df$lat, na.rm = TRUE), "to", max(presence_df$lat, na.rm = TRUE), "\n")

# ── 3. Load Bavaria boundary ──────────────────────────────────────────────────
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

# ── 4. Define bounding box for Bavaria ───────────────────────────────────────
bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)

# ── 5. Generate absence points ────────────────────────────────────────────────
set.seed(42)

presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

n_absence <- nrow(presence_df) * 2
absence_pts <- st_sample(bavaria_sf, size = n_absence * 4) %>%
  st_as_sf() %>%
  st_set_crs(4326)

# Distance in meters with sf s2 geodesic calculations
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

# ── 5. Cluster presence points within 2 km ─────────────────────────────────────
# (Preserve draw_points.r behavior)
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
  if (ra != rb) {
    parent[rb] <<- ra
  }
}
for (i in seq_len(n)) {
  nbrs <- neighbors[[i]]
  if (length(nbrs) > 1) {
    for (j in nbrs[nbrs > i]) {
      union_parent(i, j)
    }
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

# ── 6. Create map ─────────────────────────────────────────────────────────────
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

# ── 7. Display and save ───────────────────────────────────────────────────────
print(p)
ggsave("draw_abs_and_pre_map.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Map saved: draw_abs_and_pre_map.png\n")
