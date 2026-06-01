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

presence_df <- df %>%
  select(lng, lat, all_of(feature_cols), presence) %>%
  drop_na()

cat("Verwertbare Fund-Datensätze:", nrow(presence_df), "\n")

# Reduziere die Anzahl der Präsenzpunkte für schnelleres Testen
presence_df <- presence_df %>%
  slice_sample(n = 2000)

cat("Nach Stichprobe (2000 Präsenzpunkte):", nrow(presence_df), "\n")

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

absence_pts    <- absence_pts[min_dist_m > 5000, ][1:n_absence, ]
absence_coords <- st_coordinates(absence_pts) %>%
  as.data.frame() %>%
  rename(lng = X, lat = Y)

print("Gefilterte Absence-Koordinaten (mind. 5 km von Präsenzpunkten entfernt):")
print(head(absence_coords))

absence_features <- absence_coords %>%
  mutate(
    Höhe_SRTM1_puffer50m            = mean(presence_df$Höhe_SRTM1_puffer50m,            na.rm = TRUE) + rnorm(n(), 0, 50),
    Neigung_SRTM1_puffer50m         = pmax(0, mean(presence_df$Neigung_SRTM1_puffer50m,  na.rm = TRUE) + rnorm(n(), 0, 2)),
    Hangausrichtung_SRTM1_puffer50m = runif(n(), 0, 360),
    Loess_1zu500k_puffer50m         = pmax(0, mean(presence_df$Loess_1zu500k_puffer50m,  na.rm = TRUE) + rnorm(n(), 0, 500)),
    Wasser_puffer50m                = pmax(0, mean(presence_df$Wasser_puffer50m,         na.rm = TRUE) + rnorm(n(), 0, 200)),
    Viewshed_km2                    = pmax(0, mean(presence_df$Viewshed_km2,             na.rm = TRUE) + rnorm(n(), 0, 5)),
    Umfeldanalyse_km2               = pmax(0, mean(presence_df$Umfeldanalyse_km2,        na.rm = TRUE) + rnorm(n(), 0, 5)),
    Reliefenergie                   = pmax(0, mean(presence_df$Reliefenergie,            na.rm = TRUE) + rnorm(n(), 0, 10)),
    Frosttage_Jahr                  = pmax(0, mean(presence_df$Frosttage_Jahr,           na.rm = TRUE) + rnorm(n(), 0, 10)),
    Niederschlag_Jahr               = pmax(0, mean(presence_df$Niederschlag_Jahr,        na.rm = TRUE) + rnorm(n(), 0, 50)),
    Sonnenstunden_Jahr              = pmax(0, mean(presence_df$Sonnenstunden_Jahr,       na.rm = TRUE) + rnorm(n(), 0, 50)),
    Temperatur_Jahr                 = mean(presence_df$Temperatur_Jahr,                  na.rm = TRUE) + rnorm(n(), 0, 0.5),
    presence = 0L
  )

print("Absence-Features (erstes paar Zeilen):")
print(head(absence_features))

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
    min.node.size = c(5, 10)
  ),
  importance = "impurity"
)

print(rf_model)
cat("\nBeste Parameter:\n"); print(rf_model$bestTune)
cat("\nVariable Importance:\n"); print(varImp(rf_model))

# ── 6. Vorhersage-Raster erstellen ────────────────────────────────────────────
res <- 0.02

pred_grid <- expand.grid(
  lng = seq(bbox_bavaria["xmin"], bbox_bavaria["xmax"], by = res),
  lat = seq(bbox_bavaria["ymin"], bbox_bavaria["ymax"], by = res)
)

grid_sf    <- st_as_sf(pred_grid, coords = c("lng", "lat"), crs = 4326)
in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[, 1]
pred_grid  <- pred_grid[in_bavaria, ]

# Überprüfen: Wie viele Punkte liegen im Bayern-Polygon?
sum(in_bavaria)

pred_grid <- pred_grid %>%
  mutate(
    Höhe_SRTM1_puffer50m            = mean(train_df$Höhe_SRTM1_puffer50m,            na.rm = TRUE),
    Neigung_SRTM1_puffer50m         = mean(train_df$Neigung_SRTM1_puffer50m,         na.rm = TRUE),
    Hangausrichtung_SRTM1_puffer50m = mean(train_df$Hangausrichtung_SRTM1_puffer50m, na.rm = TRUE),
    Loess_1zu500k_puffer50m         = mean(train_df$Loess_1zu500k_puffer50m,         na.rm = TRUE),
    Wasser_puffer50m                = mean(train_df$Wasser_puffer50m,                na.rm = TRUE),
    Viewshed_km2                    = mean(train_df$Viewshed_km2,                    na.rm = TRUE),
    Umfeldanalyse_km2               = mean(train_df$Umfeldanalyse_km2,               na.rm = TRUE),
    Reliefenergie                   = mean(train_df$Reliefenergie,                   na.rm = TRUE),
    Frosttage_Jahr                  = mean(train_df$Frosttage_Jahr,                  na.rm = TRUE),
    Niederschlag_Jahr               = mean(train_df$Niederschlag_Jahr,               na.rm = TRUE),
    Sonnenstunden_Jahr              = mean(train_df$Sonnenstunden_Jahr,              na.rm = TRUE),
    Temperatur_Jahr                 = mean(train_df$Temperatur_Jahr,                 na.rm = TRUE)
  )

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

save_name <- "heatmap_eisenzeit_bavaria_3"
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
    subtitle = "Modell: Random Forest | Gelbe Punkte = bekannte Fundorte",
    caption  = "Datenbasis: VFPA Eisenzeit, BLfD | Fender 2017 | Data Challenge SS2026",
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
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