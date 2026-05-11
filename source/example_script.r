# =============================================================================
# Predictive Modelling – Archäologische Fundstellen (Eisenzeit Bayern)
# Data Challenge SS 2026 – Goethe-Universität Frankfurt
#
# Datenbasis: ffm_vfpa_eisenzeit.sql  (Tabelle `4ffm`)
# Methode:    Random Forest (Binärklassifikation: Fund / kein Fund)
#             + Heatmap-Ausgabe via terra
#
# Voraussetzungen (einmalig installieren):
#   install.packages(c("RMySQL","DBI","terra","sf","ranger","ggplot2",
#                      "tidyverse","caret","viridis","rnaturalearth",
#                      "rnaturalearthdata"))
# =============================================================================

# ── 0. Pakete ─────────────────────────────────────────────────────────────────
library(DBI)
library(RMySQL)       # oder RMariaDB::MariaDB() falls bevorzugt
library(tidyverse)
library(sf)
library(terra)
library(ranger)       # schneller Random Forest
library(caret)
library(viridis)
library(rnaturalearth)
library(rnaturalearthdata)

# ── 1. Datenbankverbindung & Daten laden ──────────────────────────────────────
# Passe Host / Port / Credentials an deine lokale MySQL-Instanz an.
# Falls du den SQL-Dump direkt importiert hast (empfohlen):
#   mysql -u root -p < ffm_vfpa_eisenzeit.sql
#
# Alternativ: direkt aus CSV exportieren und read_csv() nutzen (siehe Abschnitt 1b).

con <- dbConnect(
  RMySQL::MySQL(),
  host     = "127.0.0.1",
  port     = 3306,
  dbname   = "vfpa_eisenzeit",   # oder "bayern_eisenzeit" je nach Import
  user     = "root",
  password = ""                  # ← eigenes Passwort eintragen
)

raw <- dbGetQuery(con, "SELECT * FROM `4ffm`")
dbDisconnect(con)

# ── 1b. Alternative: CSV (kein MySQL nötig) ───────────────────────────────────
# Exportiere die Tabelle z.B. mit:
#   SELECT * FROM `4ffm` INTO OUTFILE '/tmp/4ffm.csv' FIELDS TERMINATED BY ','
#             ENCLOSED BY '"' LINES TERMINATED BY '\n';
# raw <- read_csv("/tmp/4ffm.csv", locale = locale(decimal_mark = ","))

# ── 2. Daten bereinigen & Features aufbereiten ────────────────────────────────
# Die numerischen Spalten liegen als VARCHAR mit deutschem Komma (z.B. "1,23") vor.
# Konvertierung: Komma → Punkt, dann as.numeric().

fix_decimal <- function(x) as.numeric(gsub(",", ".", x))

df <- raw %>%
  filter(!is.na(lng_wgs84), !is.na(lat_wgs84)) %>%
  mutate(across(
    c(Höhe_SRTM1_puffer50m, Neigung_SRTM1_puffer50m, Hangausrichtung_SRTM1_puffer50m,
      Loess_1zu500k_puffer50m, Wasser_puffer50m, Umfeldanalyse_km2, Viewshed_km2,
      Reliefenergie, Frosttage_Jahr, Niederschlag_Jahr, Niederschlag_Sommer,
      `Niederschlag_Frühling`, Sonnenstunden_Jahr, Sonnenstunden_Frühling,
      Sonnenstunden_Sommer, Temperatur_Jahr, Temperatur_Sommer, `Temperatur_Frühling`),
    fix_decimal
  )) %>%
  mutate(
    lng = as.numeric(lng_wgs84),
    lat = as.numeric(lat_wgs84),
    # Zielvariable: jede Zeile ist ein echter Fund → presence = 1
    presence = 1L
  )

# Feature-Spalten (Prädiktoren)
feature_cols <- c(
  "Höhe_SRTM1_puffer50m",        # Höhenlage (m ü. NN)
  "Neigung_SRTM1_puffer50m",     # Hangneigung (%)
  "Hangausrichtung_SRTM1_puffer50m",  # Exposition (°)
  "Loess_1zu500k_puffer50m",     # Löss-Bedeckung (Bodengüte-Proxy, m²)
  "Wasser_puffer50m",            # Nähe zu Wasser (m²)
  "Viewshed_km2",                # Sichtbarkeit (km²)
  "Umfeldanalyse_km2",           # Nutzbare Fläche im Umfeld (km²)
  "Reliefenergie",               # Reliefenergie
  "Frosttage_Jahr",              # Frosttage pro Jahr
  "Niederschlag_Jahr",           # Jahresniederschlag (mm)
  "Sonnenstunden_Jahr",          # Sonnenstunden pro Jahr
  "Temperatur_Jahr"              # Jahresdurchschnittstemperatur (°C)
)

presence_df <- df %>%
  select(lng, lat, all_of(feature_cols), presence) %>%
  drop_na()

cat("Anzahl verwertbarer Fund-Datensätze:", nrow(presence_df), "\n")

# ── 3. Pseudo-Absence generieren ──────────────────────────────────────────────
# Für Binärklassifikation brauchen wir "Nicht-Fund"-Punkte (Pseudo-Absences).
# Strategie: zufällige Punkte innerhalb Bayerns (Bounding Box), die NICHT
# in der Nähe bekannter Fundstellen liegen (min. 5 km Abstand).

set.seed(42)

# Bounding Box Bayern (WGS84)
bbox_bavaria <- c(xmin = 8.9, xmax = 13.9, ymin = 47.2, ymax = 50.6)
n_absence    <- nrow(presence_df) * 2   # 2× so viele Absences wie Presences

# Bayerische Grenze laden (für räumliche Filterung)
bavaria_sf <- ne_states(country = "Germany", returnclass = "sf") %>%
  filter(name == "Bayern")

# Zufällige Punkte innerhalb Bayerns generieren
absence_pts <- st_sample(bavaria_sf, size = n_absence * 3) %>%
  st_as_sf() %>%
  st_set_crs(4326)

# Mindestabstand zu Fundstellen (5 km)
presence_sf <- st_as_sf(presence_df, coords = c("lng", "lat"), crs = 4326)
dist_matrix  <- st_distance(absence_pts, presence_sf)
min_dist_m   <- apply(dist_matrix, 1, min)

absence_pts  <- absence_pts[min_dist_m > 5000, ][1:n_absence, ]
absence_coords <- st_coordinates(absence_pts) %>% as.data.frame()
colnames(absence_coords) <- c("lng", "lat")

# Umweltvariablen für Absence-Punkte: Mittelwert + Rauschen (Placeholder).
# → ERSETZE DIESEN BLOCK durch terra::extract() aus echten Raster-Layern
#   (z.B. SRTM, WorldClim, Corine Land Cover), sobald diese vorliegen!
absence_features <- absence_coords %>%
  mutate(
    Höhe_SRTM1_puffer50m           = mean(presence_df$Höhe_SRTM1_puffer50m,        na.rm=TRUE) + rnorm(n(), 0, 50),
    Neigung_SRTM1_puffer50m        = pmax(0, mean(presence_df$Neigung_SRTM1_puffer50m, na.rm=TRUE) + rnorm(n(), 0, 2)),
    Hangausrichtung_SRTM1_puffer50m = runif(n(), 0, 360),
    Loess_1zu500k_puffer50m        = pmax(0, mean(presence_df$Loess_1zu500k_puffer50m, na.rm=TRUE) + rnorm(n(), 0, 500)),
    Wasser_puffer50m               = pmax(0, mean(presence_df$Wasser_puffer50m,       na.rm=TRUE) + rnorm(n(), 0, 200)),
    Viewshed_km2                   = pmax(0, mean(presence_df$Viewshed_km2,           na.rm=TRUE) + rnorm(n(), 0, 5)),
    Umfeldanalyse_km2              = pmax(0, mean(presence_df$Umfeldanalyse_km2,      na.rm=TRUE) + rnorm(n(), 0, 5)),
    Reliefenergie                  = pmax(0, mean(presence_df$Reliefenergie,          na.rm=TRUE) + rnorm(n(), 0, 10)),
    Frosttage_Jahr                 = pmax(0, mean(presence_df$Frosttage_Jahr,         na.rm=TRUE) + rnorm(n(), 0, 10)),
    Niederschlag_Jahr              = pmax(0, mean(presence_df$Niederschlag_Jahr,      na.rm=TRUE) + rnorm(n(), 0, 50)),
    Sonnenstunden_Jahr             = pmax(0, mean(presence_df$Sonnenstunden_Jahr,     na.rm=TRUE) + rnorm(n(), 0, 50)),
    Temperatur_Jahr                = mean(presence_df$Temperatur_Jahr,                na.rm=TRUE) + rnorm(n(), 0, 0.5),
    presence = 0L
  )

# ── 4. Trainings-Datensatz zusammenbauen ──────────────────────────────────────
train_df <- bind_rows(
  select(presence_df, all_of(feature_cols), presence),
  select(absence_features, all_of(feature_cols), presence)
) %>%
  mutate(presence = factor(presence, levels = c(0, 1),
                           labels = c("absent", "present"))) %>%
  drop_na()

cat("Trainingsdaten: ", nrow(train_df), "Zeilen  |",
    sum(train_df$presence == "present"), "Presences  |",
    sum(train_df$presence == "absent"),  "Absences\n")

# ── 5. Modell trainieren (Random Forest via caret + ranger) ───────────────────
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
    mtry            = c(3, 5, 7),
    splitrule       = "gini",
    min.node.size   = c(5, 10)
  ),
  importance = "impurity"
)

print(rf_model)
cat("\nBeste Parameter:\n")
print(rf_model$bestTune)

# Variable Importance
cat("\nVariable Importance:\n")
varImp(rf_model) %>% print()

# ── 6. Vorhersage-Raster erstellen ────────────────────────────────────────────
# Auflösung: 0.02° ≈ ~1.5 km (anpassen je nach Bedarf / Rechenzeit)
res <- 0.02

pred_grid <- expand.grid(
  lng = seq(bbox_bavaria["xmin"], bbox_bavaria["xmax"], by = res),
  lat = seq(bbox_bavaria["ymin"], bbox_bavaria["ymax"], by = res)
)

# Räumliche Filterung auf Bayern
grid_sf <- st_as_sf(pred_grid, coords = c("lng", "lat"), crs = 4326)
in_bavaria <- st_intersects(grid_sf, bavaria_sf, sparse = FALSE)[, 1]
pred_grid  <- pred_grid[in_bavaria, ]

# Umweltvariablen für das Raster:
# WICHTIG: Ersetze hier die Mittelwert-Fallbacks durch echte Raster-Extraktion!
#
#   srtm   <- rast("srtm_bavaria.tif")
#   wasser <- rast("distance_to_water.tif")
#   loess  <- rast("loess_500k.tif")
#   ...
#   grid_terra <- vect(pred_grid, geom=c("lng","lat"), crs="EPSG:4326")
#   pred_grid$Höhe_SRTM1_puffer50m <- extract(srtm, grid_terra)[,2]
#   ...  (analog für alle Feature-Spalten)
#
# Für dieses Demo-Skript: Feature-Werte durch globale Mittelwerte simulieren.
# Das ergibt eine räumlich UNIFORME Baseline – noch keine echte Heatmap!
# Ersetze die mutate()-Blöcke unten durch die Raster-Extraktion.

pred_grid <- pred_grid %>%
  mutate(
    Höhe_SRTM1_puffer50m           = mean(train_df$Höhe_SRTM1_puffer50m,           na.rm=TRUE),
    Neigung_SRTM1_puffer50m        = mean(train_df$Neigung_SRTM1_puffer50m,        na.rm=TRUE),
    Hangausrichtung_SRTM1_puffer50m = mean(train_df$Hangausrichtung_SRTM1_puffer50m, na.rm=TRUE),
    Loess_1zu500k_puffer50m        = mean(train_df$Loess_1zu500k_puffer50m,        na.rm=TRUE),
    Wasser_puffer50m               = mean(train_df$Wasser_puffer50m,               na.rm=TRUE),
    Viewshed_km2                   = mean(train_df$Viewshed_km2,                   na.rm=TRUE),
    Umfeldanalyse_km2              = mean(train_df$Umfeldanalyse_km2,              na.rm=TRUE),
    Reliefenergie                  = mean(train_df$Reliefenergie,                  na.rm=TRUE),
    Frosttage_Jahr                 = mean(train_df$Frosttage_Jahr,                 na.rm=TRUE),
    Niederschlag_Jahr              = mean(train_df$Niederschlag_Jahr,              na.rm=TRUE),
    Sonnenstunden_Jahr             = mean(train_df$Sonnenstunden_Jahr,             na.rm=TRUE),
    Temperatur_Jahr                = mean(train_df$Temperatur_Jahr,                na.rm=TRUE)
  )

# Wahrscheinlichkeits-Vorhersage
pred_prob <- predict(rf_model,
                     newdata = select(pred_grid, all_of(feature_cols)),
                     type    = "prob")[, "present"]
pred_grid$prob_present <- pred_prob

# ── 7. terra-Raster bauen & speichern ─────────────────────────────────────────
r <- rast(
  xmin = bbox_bavaria["xmin"], xmax = bbox_bavaria["xmax"],
  ymin = bbox_bavaria["ymin"], ymax = bbox_bavaria["ymax"],
  resolution = res,
  crs = "EPSG:4326"
)

# Koordinaten → Zell-Index
cell_idx <- cellFromXY(r, as.matrix(pred_grid[, c("lng", "lat")]))
r[cell_idx] <- pred_grid$prob_present

names(r) <- "fund_wahrscheinlichkeit"

# GeoTIFF exportieren (z.B. für QGIS / ArcGIS)
writeRaster(r, "heatmap_eisenzeit_bavaria.tif", overwrite = TRUE)
cat("GeoTIFF gespeichert: heatmap_eisenzeit_bavaria.tif\n")

# ── 8. Visualisierung mit ggplot2 ─────────────────────────────────────────────
# Raster → Data Frame für ggplot
heatmap_df <- as.data.frame(r, xy = TRUE) %>%
  rename(lng = x, lat = y, prob = fund_wahrscheinlichkeit) %>%
  filter(!is.na(prob))

# Bayerische Grenzen für Overlay
bavaria_outline <- ne_states(country = "Germany", returnclass = "sf") %>%
  filter(name == "Bayern")

p <- ggplot() +
  geom_raster(data = heatmap_df, aes(x = lng, y = lat, fill = prob)) +
  geom_sf(data = bavaria_outline, fill = NA, color = "white", linewidth = 0.4) +
  # Bekannte Fundstellen einzeichnen
  geom_point(data = presence_df,
             aes(x = lng, y = lat),
             color = "yellow", size = 0.8, alpha = 0.6, shape = 16) +
  scale_fill_viridis_c(
    option  = "magma",
    name    = "Auffindungs-\nwahrscheinlichkeit",
    limits  = c(0, 1),
    breaks  = seq(0, 1, 0.2),
    labels  = scales::percent
  ) +
  coord_sf(
    xlim = c(bbox_bavaria["xmin"], bbox_bavaria["xmax"]),
    ylim = c(bbox_bavaria["ymin"], bbox_bavaria["ymax"]),
    expand = FALSE
  ) +
  labs(
    title    = "Predictive Modelling – Eisenzeitliche Fundstellen Bayern",
    subtitle = "Modell: Random Forest | Gelbe Punkte = bekannte Fundorte (4ffm-Tabelle)",
    caption  = "Datenbasis: VFPA Eisenzeit, BLfD | Fender 2017 | Data Challenge SS2026",
    x = "Längengrad", y = "Breitengrad"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title    = element_text(face = "bold"),
    legend.position = "right"
  )

ggsave("heatmap_eisenzeit_bavaria.png", plot = p,
       width = 12, height = 9, dpi = 300)
cat("PNG gespeichert: heatmap_eisenzeit_bavaria.png\n")
print(p)

# ── 9. Modell-Evaluation (Kreuzvalidierung) ───────────────────────────────────
cat("\n── Modell-Performance (5-fold CV) ──────────────────────\n")
cv_results <- rf_model$results %>%
  filter(mtry == rf_model$bestTune$mtry,
         min.node.size == rf_model$bestTune$min.node.size)
cat(sprintf("  ROC-AUC:    %.3f\n", cv_results$ROC))
cat(sprintf("  Sensitivity:%.3f\n", cv_results$Sens))
cat(sprintf("  Specificity:%.3f\n", cv_results$Spec))

# Konfusionsmatrix auf Hold-out-Predictions
conf_mat <- confusionMatrix(
  rf_model$pred$pred,
  rf_model$pred$obs,
  positive = "present"
)
print(conf_mat)

# ── 10. Nächste Schritte / Hinweise ──────────────────────────────────────────
# 1. Echte Raster-Layer einbinden (Abschnitt 6):
#    - DEM/SRTM:        https://earthexplorer.usgs.gov/
#    - Weltklima:       WorldClim v2  (worldclim.org)
#    - Corine LC:       Copernicus Land Service
#    - Wasserabstand:   aus OpenStreetMap oder BKG-Daten ableiten
#    - Löss (MSQR):     BZB / BGR GeoViewer
#
# 2. Alternative Modelle testen:
#    - MaxENT (ENMeval-Paket) – klassisch in der Archäologie
#    - GLM / GAM (Binomial) – für interpretierbarere Koeffizienten
#    - BRT / XGBoost        – oft bessere AUC als RF
#
# 3. Räumliche Kreuzvalidierung (blockCV-Paket) statt zufälliger CV verwenden,
#    um spatial autocorrelation bias zu vermeiden.
#
# 4. Perioden-spezifische Teilmodelle (HA = Hallstatt, LT = Latène):
#    train_ha <- filter(df, grepl("^HA", periode))
#    train_lt <- filter(df, grepl("^LT", periode))
#
# 5. Heatmap exportieren → QGIS öffnen, Legende anpassen, Karte finalisieren.