############################################################
# Download von Umweltdaten für Bayern
############################################################

options(repos = c(CRAN = "https://cloud.r-project.org"))

# Pakete ---------------------------------------------------
if (!require("terra")) install.packages("terra")
if (!require("geodata")) install.packages("geodata")

library(terra)
library(geodata)

# Ordner anlegen -------------------------------------------
dir.create("data", showWarnings = FALSE)
dir.create("data/worldclim", showWarnings = FALSE)
dir.create("data/dem", showWarnings = FALSE)

############################################################
# WORLDCLIM
############################################################

# Auflösung:
# 10 = ca. 18 km
# 5  = ca. 9 km
# 2.5 = ca. 5 km
# 0.5 = ca. 1 km

resolution <- 0.5

cat("Lade WorldClim Temperatur...\n")

tavg <- worldclim_global(
  var = "tavg",
  res = resolution,
  path = "data/worldclim"
)

cat("Lade WorldClim Niederschlag...\n")

prec <- worldclim_global(
  var = "prec",
  res = resolution,
  path = "data/worldclim"
)

cat("Lade WorldClim Bioklima...\n")

bio <- worldclim_global(
  var = "bio",
  res = resolution,
  path = "data/worldclim"
)

############################################################
# DEM
############################################################

cat("Lade Höhenmodell...\n")

dem <- elevation_30s(
  country = "DEU",
  path = "data/dem"
)

############################################################
# Bayern ausschneiden
############################################################

# Bounding Box Bayern
bayern <- ext(
  8.9,
  13.9,
  47.2,
  50.6
)

tavg <- crop(tavg, bayern)
prec <- crop(prec, bayern)
bio  <- crop(bio,  bayern)
dem  <- crop(dem,  bayern)

############################################################
# Terrainvariablen
############################################################

terrain_layers <- terrain(
  dem,
  v = c("slope", "aspect"),
  unit = "degrees"
)

slope  <- terrain_layers[[1]]
aspect <- terrain_layers[[2]]

############################################################
# GeoTIFF speichern
############################################################

cat("Speichere GeoTIFFs...\n")

writeRaster(
  tavg,
  "data/worldclim/temperature_monthly.tif",
  overwrite = TRUE
)

writeRaster(
  prec,
  "data/worldclim/precipitation_monthly.tif",
  overwrite = TRUE
)

writeRaster(
  bio,
  "data/worldclim/bioclim_variables.tif",
  overwrite = TRUE
)

writeRaster(
  dem,
  "data/dem/dem_bayern.tif",
  overwrite = TRUE
)

writeRaster(
  slope,
  "data/dem/slope.tif",
  overwrite = TRUE
)

writeRaster(
  aspect,
  "data/dem/aspect.tif",
  overwrite = TRUE
)

cat("======================================\n")
cat("Alle GeoTIFFs wurden gespeichert.\n")
cat("======================================\n")