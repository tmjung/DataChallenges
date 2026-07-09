# Setup & Pakete
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c("terra", "geodata", "sf"))

library(terra)
library(geodata)
library(sf)

# Bayern Boundary
bavaria_all <- gadm("DEU", level = 1, path = "data")
print(names(bavaria_all))

bavaria <- bavaria_all[bavaria_all$NAME_1 == "Bayern", ]

bav_ext <- ext(bavaria)

#Temperatur + Niederschlag + DEM
dir.create("data/climate", recursive = TRUE, showWarnings = FALSE)
dir.create("data/dem", recursive = TRUE, showWarnings = FALSE)

# Temperatur (Jahresmittel)
tavg <- worldclim_global(
  var = "tavg",
  res = 0.5,
  path = "data/climate"
)

# Niederschlag
prec <- worldclim_global(
  var = "prec",
  res = 0.5,
  path = "data/climate"
)

# Crop auf Bayern
tavg <- crop(tavg, bav_ext)
prec <- crop(prec, bav_ext)

# DEM + Terrain
dem <- elevation_30s(country = "DEU", path = "data/dem")

dem <- crop(dem, bav_ext)

terrain_layers <- terrain(dem, c("slope", "aspect"), unit = "degrees")

slope <- terrain_layers[[1]]
aspect <- terrain_layers[[2]]

# GeoTIFF speichern
writeRaster(tavg, "data/climate/temperature.tif", overwrite = TRUE)
writeRaster(prec, "data/climate/precipitation.tif", overwrite = TRUE)

writeRaster(dem, "data/dem/dem.tif", overwrite = TRUE)
writeRaster(slope, "data/dem/slope.tif", overwrite = TRUE)
writeRaster(aspect, "data/dem/aspect.tif", overwrite = TRUE)