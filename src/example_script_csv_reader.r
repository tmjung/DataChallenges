# =============================================================================
# Predictive Modelling – Archäologische Fundstellen (Eisenzeit Bayern)
# Data Challenge SS 2026 – Goethe-Universität Frankfurt
#
# Datenbasis: ffm_vfpa_eisenzeit.csv  (direkt, kein MySQL nötig)
# Methode:    Random Forest (Binärklassifikation: Fund / kein Fund)
#             + Heatmap-Ausgabe via terra
# =============================================================================

# ── 0. Pakete installieren (einmalig) ─────────────────────────────────────────
# install.packages(c(
#   "tidyverse", "sf", "terra", "ranger", "caret",
#   "viridis", "rnaturalearth", "rnaturalearthdata"
# ), dependencies = TRUE, repos = "https://cloud.r-project.org")

# ── 0b. Pakete laden ──────────────────────────────────────────────────────────
library(tidyverse)
library(sf)
library(terra)
library(ranger)
library(caret)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# ── 1. CSV-Dateien laden ──────────────────────────────────────────────────────
# Passe die Pfade an, falls deine CSVs in einem Unterordner liegen.
# Beispiel: "data/ffm_vfpa_eisenzeit.csv"

# Hauptdatensatz (Komma-getrennt, Dezimalkomma in Anführungszeichen)
raw <- read_csv(
  "data/ffm_vfpa_eisenzeit.csv",
  locale = locale(decimal_mark = ","),
  show_col_types = FALSE
)

# Fundort-Geonames (Semikolon-getrennt)
fundorte <- read_csv2(
  "data/export_Immovables_15.4.2026-tchi1_Fundorte_Geonames_Link.csv",
  show_col_types = FALSE
)

# Münzdaten (Semikolon-getrennt)
muenzen <- read_csv2(
  "data/export_Numismatic object_14.4.2026-numisdata4_Muenzen_Fundort_Kontext.csv",
  show_col_types = FALSE
)

cat("Zeilen geladen - ffm:", nrow(raw),
    "| Fundorte:", nrow(fundorte),
    "| Münzen:", nrow(muenzen), "\n")

# ── 2. Daten bereinigen & Features aufbereiten ────────────────────────────────
# Numerische Spalten liegen als VARCHAR mit deutschem Komma vor ("1,23").
# read_csv mit decimal_mark="," löst das bereits für viele Spalten.
# Für Spalten die dennoch als character eingelesen werden: manuell konvertieren.

fix_decimal <- function(x) {
  if (is.character(x)) as.numeric(gsub(",", ".", x))
  else x
}

df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(across(
    c(Höhe_SRTM1_puffer50m, Neigung_SRTM1_puffer50m,
      Hangausrichtung_SRTM1_puffer50m,
      Loess_1zu500k_puffer50m, Wasser_puffer50m,
      Umfeldanalyse_km2, Viewshed_km2,
      Reliefenergie, Frosttage_Jahr,
      Niederschlag_Jahr, Niederschlag_Sommer,
      `Niederschlag_Frühling`, Sonnenstunden_Jahr,
      Sonnenstunden_Frühling, Sonnenstunden_Sommer,
      Temperatur_Jahr, Temperatur_Sommer,
      `Temperatur_Frühling`),
    fix_decimal
  )) %>%
  mutate(
    lng      = as.numeric(lng_wgs84),
    lat      = as.numeric(lat_wgs84),
    presence = 1L
  )

# Feature-Spalten (Prädiktoren)
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
feature_cols <- base_feature_cols
model_feature_cols <- c(base_feature_cols, spatial_feature_cols)

presence_df <- df %>%
  select(lng, lat, all_of(feature_cols), presence) %>%
  drop_na()

cat("Verwertbare Fund-Datensätze:", nrow(presence_df), "\n")

# Reduziere die Anzahl der Präsenzpunkte für schnelleres Testen
presence_df <- presence_df %>%
  slice_sample(n = 2000)

cat("Nach Stichprobe (2000 Präsenzpunkte):", nrow(presence_df), "\n")

# Räumliche Zusatzfeatures für Präsenzpunkte berechnen
presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)
dist_matrix_pres <- st_distance(presence_sf, presence_sf)
dist_matrix_pres <- matrix(as.numeric(dist_matrix_pres), nrow = nrow(presence_df))
diag(dist_matrix_pres) <- Inf

presence_df <- presence_df %>%
  mutate(
    dist_nearest_presence_km = apply(dist_matrix_pres, 1, min) / 1000,
    presence_count_10km      = apply(dist_matrix_pres, 1, function(d) sum(d <= 10000))
  )

# Jetzt die neuen räumlichen Features in die Trainingsspalten aufnehmen
feature_cols <- model_feature_cols

# ── 3. Pseudo-Absence generieren ──────────────────────────────────────────────
set.seed(42)

bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)
n_absence    <- nrow(presence_df) * 2

# PROBLEM: Bei scale=50 ist Deuitschland nicht enthalten

states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

#unique(states_raw$admin)

germany_sf <- states_raw %>%
  filter(admin == "Germany")

bavaria_sf <- germany_sf |>
  filter(name_de == "Bayern")

# Fallback: Nutzt rechteck als Grenzen
#bavaria_sf <- st_as_sfc(st_bbox(c(
#  xmin = 8.9, xmax = 13.9,
#  ymin = 47.2, ymax = 50.6
#), crs = st_crs(4326))) %>%
#  st_as_sf()

print("Bayern SF-Objekt:")
print(bavaria_sf)

absence_pts <- st_sample(bavaria_sf, size = n_absence * 3) %>%
  st_as_sf() %>%
  st_set_crs(4326)

print("Negative Punkte (sf):")
print(absence_pts)

presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)
dist_matrix  <- st_distance(absence_pts, presence_sf)
min_dist_m   <- apply(dist_matrix, 1, min)

print("Abstände der Absence-Punkte zu nächsten Presence-Punkten (m):")
print(summary(min_dist_m))

selected_rows <- which(min_dist_m > 5000)[1:n_absence]
absence_pts    <- absence_pts[selected_rows, ]
absence_coords <- st_coordinates(absence_pts) %>%
  as.data.frame() %>%
  rename(lng = X, lat = Y)

print("Gefilterte Absence-Koordinaten (mind. 5 km von Präsenzpunkten entfernt):")
print(head(absence_coords))

# Für jeden gefilterten Absence-Punkt: robustere, weniger deterministische
# Ableitung der Prädiktor-Werte. Statt immer genau die 3 nächste Präsenz-Mittelwerte
# zu verwenden, ziehen wir die 10 nächsten Präsenzpunkte zur Auswahl und bilden
# für jede Absenz zufällig ein Mittel aus 3 davon. Zusätzlich fügen wir kleinen
# Gaußschen Rauschen (5% der Präsenz-SD) hinzu, um identische Feature-Vektoren
# zu vermeiden, die das RF stark überfitten können.
pres_feat_matrix <- as.matrix(select(presence_df, all_of(base_feature_cols)))
pres_feat_sd <- apply(pres_feat_matrix, 2, sd, na.rm = TRUE)

nearest_full_ids <- t(apply(dist_matrix[selected_rows, , drop = FALSE], 1, function(d) order(d)[1:10]))

absence_values <- t(apply(nearest_full_ids, 1, function(ord_row) {
  # wähle die 3 nächsten Präsenzpunkte aus dem Pool der 10 nächsten aus
  k_use  <- min(3, length(ord_row))
  pick_ids <- sample(ord_row, k_use)
  mu <- colMeans(pres_feat_matrix[pick_ids, , drop = FALSE], na.rm = TRUE)
  # Rauschen proportional zur SD (kleiner Faktor)
  noise <- rnorm(length(mu), mean = 0, sd = pmax(1e-6, pres_feat_sd * 0.05))
  mu + noise
}))

absence_values <- as.data.frame(absence_values)
colnames(absence_values) <- base_feature_cols

absence_dist_matrix <- matrix(as.numeric(dist_matrix[selected_rows, , drop = FALSE]),
                               nrow = length(selected_rows))
absence_dist_km <- apply(absence_dist_matrix, 1, min) / 1000
absence_count_10km <- apply(absence_dist_matrix, 1, function(d) sum(d <= 10000))

absence_features <- bind_cols(absence_coords, absence_values) %>%
  mutate(
    dist_nearest_presence_km = absence_dist_km,
    presence_count_10km      = absence_count_10km,
    presence = 0L
  )

print("Absence-Features (erstes paar Zeilen):")
print(head(absence_features))

# ── 3b. Explorative Karten vor der Modellvorhersage erstellen ────────────────
extra_plot_cols <- c(
  "Temperatur_Sommer", "Niederschlag_Sommer", "Sonnenstunden_Sommer",
  "Temperatur_Frühling", "Niederschlag_Frühling", "Sonnenstunden_Frühling"
)
# Get seasonal columns from df and merge with sampled presence_df
seasonal_data <- df %>%
  select(lng, lat, all_of(extra_plot_cols))

# Use sampled presence_df (2000 points) for consistency with model data
presence_plot_df <- presence_df %>%
  select(lng, lat, all_of(base_feature_cols)) %>%
  left_join(seasonal_data, by = c("lng", "lat")) %>%
  mutate(type = "presence")

# For absence, just use coords (no feature values needed for mapping)
absence_plot_df <- absence_features %>%
  select(lng, lat) %>%
  mutate(type = "absence")

cat(sprintf("Presence points for mapping: %d\n", nrow(presence_plot_df)))
cat(sprintf("Absence points for mapping: %d\n", nrow(absence_plot_df)))

plot_points_df <- bind_rows(presence_plot_df, absence_plot_df)

save_point_map <- function(data, title, subtitle, color_var = NULL,
                           color_label = NULL, file_name) {
  p <- ggplot() +
    geom_sf(data = bavaria_sf, fill = NA, color = "black", size = 0.3) +
    geom_point(data = filter(data, type == "absence"),
               aes(x = lng, y = lat),
               color = "grey40", size = 0.8, alpha = 0.5) +
    {
      if (is.null(color_var)) {
        geom_point(data = filter(data, type == "presence"),
                   aes(x = lng, y = lat, color = type),
                   size = 1.0, alpha = 0.8)
      } else {
        geom_point(data = filter(data, type == "presence"),
                   aes(x = lng, y = lat, color = .data[[color_var]]),
                   size = 1.0, alpha = 0.8)
      }
    } +
    coord_sf(xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
             ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
             expand = FALSE) +
    labs(title = title,
         subtitle = subtitle,
         color = ifelse(is.null(color_label), "type", color_label),
         x = "Längengrad",
         y = "Breitengrad") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold"))

  if (!is.null(color_var)) {
    p <- p + scale_color_viridis_c(option = "magma", na.value = "grey50")
  } else {
    p <- p + scale_color_manual(values = c("presence" = "yellow", "absence" = "blue"))
  }

  ggsave(file_name, plot = p, width = 10, height = 8, dpi = 300)
  cat("Saved map:", file_name, "\n")
}

save_feature_pair_map <- function(df, var1, var2, label1, label2, file_name) {
  long_df <- df %>%
    select(lng, lat, type, all_of(c(var1, var2))) %>%
    pivot_longer(cols = all_of(c(var1, var2)),
                 names_to = "feature",
                 values_to = "value") %>%
    mutate(feature = recode(feature,
                            !!var1 := label1,
                            !!var2 := label2))

  p <- ggplot() +
    geom_sf(data = bavaria_sf, fill = NA, color = "black", size = 0.3) +
    geom_point(data = filter(long_df, type == "absence"),
               aes(x = lng, y = lat),
               color = "grey40", size = 0.7, alpha = 0.4) +
    geom_point(data = filter(long_df, type == "presence"),
               aes(x = lng, y = lat, color = value),
               size = 0.9, alpha = 0.8) +
    facet_wrap(~ feature, scales = "fixed") +
    coord_sf(xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
             ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
             expand = FALSE) +
    scale_color_viridis_c(option = "magma", na.value = "grey50") +
    labs(title = paste0("Feature maps: ", label1, " and ", label2),
         subtitle = "Präsenzpunkte koloriert; Absenzpunkte grau",
         color = "Wert",
         x = "Längengrad",
         y = "Breitengrad") +
    theme_minimal(base_size = 12) +
    theme(panel.grid = element_blank(),
          legend.position = "right",
          plot.title = element_text(face = "bold"))

  ggsave(file_name, plot = p, width = 12, height = 8, dpi = 300)
  cat("Saved paired feature map:", file_name, "\n")
}

cat("Aspect:\n")
summary(presence_plot_df$Hangausrichtung_SRTM1_puffer50m)

cat("NAs:\n")
sum(is.na(presence_plot_df$Hangausrichtung_SRTM1_puffer50m))

cat("Rows:\n")
nrow(filter(presence_plot_df,
            !is.na(Hangausrichtung_SRTM1_puffer50m)))

ggplot(
  presence_plot_df,
  aes(lng, lat)
) +
  geom_point(size = 0.5)

save_point_map(plot_points_df,
              title = "Presence and Absence Points",
              subtitle = "Präsenzpunkte gelb, Absenzpunkte grau",
              color_var = NULL,
              file_name = "map_presence_absence.png")

save_point_map(presence_plot_df,
              title = "Height and Presence Points",
              subtitle = "Farbskala zeigt Höhe (m)",
              color_var = "Höhe_SRTM1_puffer50m",
              color_label = "Höhe (m)",
              file_name = "map_height_points.png")

save_point_map(presence_plot_df,
              title = "Temperature and Presence Points",
              subtitle = "Farbskala zeigt jährliche Temperatur (°C)",
              color_var = "Temperatur_Jahr",
              color_label = "Temperatur (°C)",
              file_name = "map_temperature_points.png")

save_point_map(presence_plot_df,
              title = "Summer Temperature and Presence Points",
              subtitle = "Farbskala zeigt Sommertemperatur (°C)",
              color_var = "Temperatur_Sommer",
              color_label = "Temperatur Sommer (°C)",
              file_name = "map_temperature_summer_points.png")

save_point_map(presence_plot_df,
              title = "Spring Temperature and Presence Points",
              subtitle = "Farbskala zeigt Frühlingstemperatur (°C)",
              color_var = "Temperatur_Frühling",
              color_label = "Temperatur Frühling (°C)",
              file_name = "map_temperature_spring_points.png")

save_point_map(presence_plot_df,
              title = "Sun Hours and Presence Points",
              subtitle = "Farbskala zeigt jährliche Sonnenstunden",
              color_var = "Sonnenstunden_Jahr",
              color_label = "Sonnenstunden",
              file_name = "map_sunhours_points.png")

save_point_map(presence_plot_df,
              title = "Summer Sun Hours and Presence Points",
              subtitle = "Farbskala zeigt Sonnenstunden im Sommer",
              color_var = "Sonnenstunden_Sommer",
              color_label = "Sonnenstunden Sommer",
              file_name = "map_sunhours_summer_points.png")

save_point_map(presence_plot_df,
              title = "Spring Sun Hours and Presence Points",
              subtitle = "Farbskala zeigt Sonnenstunden im Frühling",
              color_var = "Sonnenstunden_Frühling",
              color_label = "Sonnenstunden Frühling",
              file_name = "map_sunhours_spring_points.png")

save_point_map(presence_plot_df,
              title = "Slope and Presence Points",
              subtitle = "Farbskala zeigt Hangneigung",
              color_var = "Neigung_SRTM1_puffer50m",
              color_label = "Neigung",
              file_name = "map_slope_points.png")

save_point_map(presence_plot_df,
              title = "Aspect and Presence Points",
              subtitle = "Farbskala zeigt Hangausrichtung",
              color_var = "Hangausrichtung_SRTM1_puffer50m",
              color_label = "Hangausrichtung",
              file_name = "map_aspect_points.png")

save_point_map(presence_plot_df,
              title = "Loess and Presence Points",
              subtitle = "Farbskala zeigt Loess-Wert",
              color_var = "Loess_1zu500k_puffer50m",
              color_label = "Loess",
              file_name = "map_loess_points.png")

save_point_map(presence_plot_df,
              title = "Water and Presence Points",
              subtitle = "Farbskala zeigt Wasser-Puffer-Wert",
              color_var = "Wasser_puffer50m",
              color_label = "Wasser",
              file_name = "map_water_points.png")

save_point_map(presence_plot_df,
              title = "Viewshed and Presence Points",
              subtitle = "Farbskala zeigt Viewshed (km²)",
              color_var = "Viewshed_km2",
              color_label = "Viewshed",
              file_name = "map_viewshed_points.png")

save_point_map(presence_plot_df,
              title = "Surroundings and Presence Points",
              subtitle = "Farbskala zeigt Umfeldanalyse (km²)",
              color_var = "Umfeldanalyse_km2",
              color_label = "Umfeldanalyse",
              file_name = "map_surroundings_points.png")

save_point_map(presence_plot_df,
              title = "Relief Energy and Presence Points",
              subtitle = "Farbskala zeigt Reliefenergie",
              color_var = "Reliefenergie",
              color_label = "Reliefenergie",
              file_name = "map_reliefenergy_points.png")

save_point_map(presence_plot_df,
              title = "Frost Days and Presence Points",
              subtitle = "Farbskala zeigt Frosttage im Jahr",
              color_var = "Frosttage_Jahr",
              color_label = "Frosttage",
              file_name = "map_frostdays_points.png")

save_point_map(presence_plot_df,
              title = "Annual Precipitation and Presence Points",
              subtitle = "Farbskala zeigt jährlichen Niederschlag",
              color_var = "Niederschlag_Jahr",
              color_label = "Niederschlag",
              file_name = "map_precipitation_points.png")

save_point_map(presence_plot_df,
              title = "Summer Precipitation and Presence Points",
              subtitle = "Farbskala zeigt Niederschlag im Sommer",
              color_var = "Niederschlag_Sommer",
              color_label = "Niederschlag Sommer",
              file_name = "map_precipitation_summer_points.png")

save_point_map(presence_plot_df,
              title = "Spring Precipitation and Presence Points",
              subtitle = "Farbskala zeigt Niederschlag im Frühling",
              color_var = "Niederschlag_Frühling",
              color_label = "Niederschlag Frühling",
              file_name = "map_precipitation_spring_points.png")

save_feature_pair_map(presence_plot_df,
                      var1 = "Hangausrichtung_SRTM1_puffer50m",
                      var2 = "Neigung_SRTM1_puffer50m",
                      label1 = "Hangausrichtung",
                      label2 = "Neigung",
                      file_name = "map_aspect_slope_points.png")

save_feature_pair_map(presence_plot_df,
                      var1 = "Temperatur_Jahr",
                      var2 = "Sonnenstunden_Jahr",
                      label1 = "Temperatur (Jahr)",
                      label2 = "Sonnenstunden (Jahr)",
                      file_name = "map_temperature_sunhours_points.png")

save_feature_pair_map(presence_plot_df,
                      var1 = "Temperatur_Sommer",
                      var2 = "Niederschlag_Sommer",
                      label1 = "Temperatur (Sommer)",
                      label2 = "Niederschlag (Sommer)",
                      file_name = "map_temperature_precipitation_summer_points.png")

# ── 4. Trainings-Datensatz zusammenbauen ──────────────────────────────────────
train_df <- bind_rows(
  select(presence_df, all_of(feature_cols), presence),
  select(absence_features, all_of(feature_cols), presence)
) %>%
  mutate(presence = factor(presence, levels = c(0, 1),
                           labels = c("absent", "present"))) %>%
  drop_na()

cat("Trainingsdaten:", nrow(train_df), "Zeilen |",
    sum(train_df$presence == "present"), "Presences |",
    sum(train_df$presence == "absent"),  "Absences\n")

# ── 5. Modell trainieren (Random Forest) ──────────────────────────────────────
set.seed(123)
ctrl <- trainControl(
  method          = "cv",
  number          = 5,
  classProbs      = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

rf_model <- train(
  presence ~ .,
  data      = train_df,
  method    = "ranger",
  trControl = ctrl,
  metric    = "ROC",
  tuneGrid  = expand.grid(
    mtry          = c(3, 5, 7),
    splitrule     = "gini",
    min.node.size = c(5, 20)
  ),
  importance = "impurity"
)

print(rf_model)
cat("\nBeste Parameter:\n"); print(rf_model$bestTune)
cat("\nVariable Importance:\n"); print(varImp(rf_model))

# ── 6. Vorhersage-Raster erstellen ────────────────────────────────────────────
res <- 0.01  # Rasterauflösung in Grad. Diese Einstellung steuert die Kachelgröße im Heatmap-Raster.

pred_grid <- expand.grid(
  lng = seq(bbox_bavaria["xmin"], bbox_bavaria["xmax"], by = res),
  lat = seq(bbox_bavaria["ymin"], bbox_bavaria["ymax"], by = res)
)

grid_sf    <- st_as_sf(pred_grid, coords = c("lng", "lat"), crs = 4326)
in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[, 1]
if (!is.logical(in_bavaria)) {
  in_bavaria <- as.logical(in_bavaria)
}
pred_grid  <- pred_grid[in_bavaria, ]

# Ensure grid_sf matches filtered pred_grid (avoid size mismatches later)
grid_sf <- grid_sf[in_bavaria, ]

cat(sprintf("in_bavaria class: %s\n", paste(class(in_bavaria), collapse=", ")))
cat(sprintf("pred_grid rows after Bavaria filter: %d\n", nrow(pred_grid)))
cat(sprintf("grid_sf rows after Bavaria filter: %d\n", nrow(grid_sf)))
cat(sprintf("in_bavaria length: %d\n", length(in_bavaria)))
cat(sprintf("sum(in_bavaria): %d\n", sum(in_bavaria)))

# Überprüfen: Wie viele Punkte liegen im Bayern-Polygon?
sum(in_bavaria)

# Für die Punkte im Bayern-Polygon: Feature-Werte aus den drei nächsten
# Präsenzpunkten ableiten, damit die Modellvorhersage räumlich variiert.
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
    presence_count_10km      = count_within_radius
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
# Optional: if environmental raster layers exist, prefer extracting real raster
# values at the prediction points. This avoids using presence-derived surrogates
# where possible. Search common folders (env, data, project root) for .tif files.
tif_files <- list.files(path = c("env", "data", "."), pattern = "\\.tif$",
                        recursive = TRUE, full.names = TRUE)
tif_files <- tif_files[file.exists(tif_files)]
if (length(tif_files) > 0) {
  cat("Found raster files, attempting terra::extract on:", paste0(tif_files, collapse = ", "), "\n")
  env_stack <- tryCatch(rast(tif_files), error = function(e) NULL)
  if (!is.null(env_stack)) {
    pred_grid_sf <- st_as_sf(pred_grid, coords = c("lng", "lat"), crs = 4326)
    env_vals <- terra::extract(env_stack, vect(pred_grid_sf))
    # terra::extract returns ID plus layers. Drop ID column
    if (ncol(env_vals) > 1) {
      env_vals2 <- as.data.frame(env_vals)[, -1, drop = FALSE]
      env_names <- names(env_vals2)
      common <- intersect(feature_cols, env_names)
      if (length(common) > 0) {
        cat("Using matching raster layers for features:", paste(common, collapse = ", "), "\n")
        pred_grid[, common] <- env_vals2[, common, drop = FALSE]
      } else if (ncol(env_vals2) == length(feature_cols)) {
        cat("Raster stack has same number of layers as feature_cols; assigning by position.\n")
        pred_grid[, feature_cols] <- env_vals2
      } else {
        cat("Raster layers didn't match feature names; keeping nearest-presence surrogates.\n")
      }
    }
  } else {
    cat("Could not read rasters with terra::rast; keeping nearest-presence surrogates.\n")
  }
} else {
  cat("No .tif rasters found in env/data/. — using nearest-presence surrogate features.\n")
}

# Debug: save pred_grid features and print summaries to help diagnose uniform predictions
debug_file <- "pred_grid_features_debug.csv"
write_csv(select(pred_grid, all_of(base_feature_cols)), debug_file)
cat("Saved prediction-grid features to:", debug_file, "\n")
cat("Feature summaries (prediction grid):\n")
for (col in feature_cols) {
  cat("--", col, "\n")
  print(summary(pred_grid[[col]]))
}

# Fast diagnostic: sample a subset of grid points and get predicted probabilities
diag_n <- min(5000, nrow(pred_grid))
set.seed(42)
diag_idx <- sample(nrow(pred_grid), diag_n)
diag_df <- pred_grid[diag_idx, , drop = FALSE]
diag_pred <- predict(rf_model, newdata = select(diag_df, all_of(feature_cols)), type = "prob")[, "present"]
cat(sprintf("Diagnostic sample predictions (%d points): min=%.4f mean=%.4f max=%.4f\n",
            length(diag_pred), min(diag_pred), mean(diag_pred), max(diag_pred)))
print(summary(diag_pred))
#pred_grid_sf <- st_as_sf(pred_grid, coords = c("lng","lat"), crs = 4326)
#env_stack <- terra::rast(c(
#  "height.tif",
#  "slope.tif",
#  "water.tif"
#))
#env_vals <- terra::extract(env_stack, vect(pred_grid_sf))
#pred_grid <- bind_cols(pred_grid, env_vals)

pred_prob              <- predict(rf_model,
                                  newdata = select(pred_grid, all_of(feature_cols)),
                                  type    = "prob")[, "present"]
pred_grid$prob_present <- pred_prob

# ── 7. terra-Raster bauen & speichern ─────────────────────────────────────────
r <- rast(
  xmin       = bbox_bavaria["xmin"], xmax = bbox_bavaria["xmax"],
  ymin       = bbox_bavaria["ymin"], ymax = bbox_bavaria["ymax"],
  resolution = res,
  crs        = "EPSG:4326"
)

cell_idx         <- cellFromXY(r, as.matrix(pred_grid[, c("lng", "lat")]))
r[cell_idx]      <- pred_grid$prob_present
names(r)         <- "fund_wahrscheinlichkeit"

# Überprüfen: Wie viele Rasterzellen wurden gefüllt? Daher die schwarzen flächen?
sum(is.na(cell_idx))
length(cell_idx)

save_name <- "heatmap_eisenzeit_bavaria_6"
save_name_tif <- paste(save_name, ".tif", sep = "")
save_name_png <- paste(save_name, ".png", sep = "")

writeRaster(r, save_name_tif, overwrite = TRUE)
cat(paste("GeoTIFF gespeichert: ", save_name_tif, "\n")) 

cat("Grid points:", nrow(pred_grid), "\n")
cat("NA cell indices:", sum(is.na(cell_idx)), "\n")
n_non_na <- global(r, "notNA", na.rm = TRUE)[1,1]
cat("Raster non-NA cells:", n_non_na, "\n")
cat("Bavaria points used:", sum(in_bavaria), "\n")

# ── 8. Visualisierung mit ggplot2 ─────────────────────────────────────────────
heatmap_df <- as.data.frame(r, xy = TRUE) %>%
  rename(lng = x, lat = y, prob = fund_wahrscheinlichkeit) %>%
  filter(!is.na(prob))

#bavaria_outline <- ne_states(country = "Germany", returnclass = "sf") %>%
#  filter(name == "Bayern")
# NEU (funktioniert ohne zusätzliche Pakete):
#bavaria_outline <- ne_states(country = "Germany", returnclass = "sf") %>%
#  filter(name == "Bayern") %>%
#  st_simplify(dTolerance = 0.01)

#bavaria_outline2 <- rnaturalearth::ne_download(
#  scale = 50, type = "states", category = "cultural",
#  returnclass = "sf"
#) %>%
#  dplyr::filter(name == "Bayern")

bavaria_outline <- germany_sf |>
  filter(name_de == "Bayern")


plot_subtitle <- sprintf(
  "Model: Random Forest | Präsenz: %d | Absenz: %d | Auflösung: %.3f°",
  sum(train_df$presence == "present"),
  sum(train_df$presence == "absent"),
  res
)

p <- ggplot() +
  geom_raster(data = heatmap_df, aes(x = lng, y = lat, fill = prob)) +
  geom_sf(data = bavaria_outline, fill = NA, color = "white", linewidth = 0.4) +
  geom_point(data = presence_df,
             aes(x = lng, y = lat),
             color = "yellow", size = 0.8, alpha = 0.6, shape = 16) +
  scale_fill_viridis_c(
    option = "magma",
    name   = "Auffindungs-\nwahrscheinlichkeit",
    limits = c(0, 1),
    breaks = seq(0, 1, 0.2),
    labels = scales::percent
  ) +
  coord_sf(
    xlim   = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim   = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE
  ) +
  labs(
    title    = "Predictive Modelling – Eisenzeitliche Fundstellen Bayern",
    subtitle = plot_subtitle,
    caption  = paste0("Datenbasis: VFPA Eisenzeit | Data Challenge SS2026 | Koo & Jungbeck",
                      " | Rasterzellen: ", nrow(pred_grid)),
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        legend.position = "right")

ggsave(save_name_png, plot = p,
       width = 12, height = 9, dpi = 300)
cat(paste("PNG gespeichert: ", save_name_png, "\n"))
print(p)

# ── 9. Modell-Evaluation ──────────────────────────────────────────────────────
cat("\n── Modell-Performance (5-fold CV) ──────────────────────\n")
cv_results <- rf_model$results %>%
  filter(mtry          == rf_model$bestTune$mtry,
         min.node.size == rf_model$bestTune$min.node.size)
cat(sprintf("  ROC-AUC:     %.3f\n", cv_results$ROC))
cat(sprintf("  Sensitivity: %.3f\n", cv_results$Sens))
cat(sprintf("  Specificity: %.3f\n", cv_results$Spec))

conf_mat <- confusionMatrix(
  rf_model$pred$pred,
  rf_model$pred$obs,
  positive = "present"
)
print(conf_mat)