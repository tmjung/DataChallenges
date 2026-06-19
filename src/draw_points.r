# =============================================================================
# Draw Map of Bavaria with Presence Points
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

cat("Presence points:", nrow(presence_df), "\n")
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

# ── 4. Define bounding box ────────────────────────────────────────────────────
bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)

# ── 5. Cluster presence points within 5 km ─────────────────────────────────────
presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)

# Transform to a projected CRS for meter distances
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
          alpha = 0.5,
          show.legend = FALSE) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE,
    crs = 4326
  ) +
  labs(
    title = "Präsenz-Cluster mit 2-km-Nähe in Bayern",
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
ggsave("draw_clusters_map.png", plot = p, width = 10, height = 8, dpi = 300)
cat("Map saved: draw_clusters_map.png\n")
