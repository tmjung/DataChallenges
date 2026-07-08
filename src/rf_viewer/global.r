################################################################################
# global.R
# Lädt alle Rasterdaten und stellt Hilfsfunktionen bereit
################################################################################

library(shiny)
library(leaflet)
library(terra)
library(sf)
library(viridis)
library(htmltools)

################################################################################
# Datenordner
################################################################################

DATA_DIR <- "data"

################################################################################
# Raster laden
################################################################################

load_raster <- function(filename) {

  path <- file.path(DATA_DIR, filename)

  if (!file.exists(path)) {

    warning(sprintf("Raster nicht gefunden: %s", path))

    return(NULL)

  }

  tryCatch(

    rast(path),

    error = function(e) {

      warning(sprintf("Fehler beim Laden von %s", filename))

      NULL

    }

  )

}

################################################################################
# Raster
################################################################################

prediction <- load_raster("rf_5_probability.tif")

precipitation <- load_raster("climate/precipitation.tif")

temperature <- load_raster("climate/temperature.tif")

dem <- load_raster("dem/dem.tif")

elevation <- load_raster("dem/DEU_elv_msk.tif")

print("\nFILE EXISTS: data/dem/DEU_elv_msk.tif\n")
print(file.exists("data/dem/DEU_elv_msk.tif"))

slope <- load_raster("dem/slope.tif")

aspect <- load_raster("dem/aspect.tif")

################################################################################
# Layerliste
################################################################################

layers <- list(

  "Fundwahrscheinlichkeit" = prediction,

  "Niederschlag" = precipitation,

  "Temperatur" = temperature,

  "Höhe (DEM)" = dem,

  "Elevation" = elevation,

  "Hangneigung" = slope,

  "Exposition" = aspect

)

################################################################################
# Bounding Box
################################################################################

if (!is.null(prediction)) {

  bbox <- ext(prediction)

} else {

  bbox <- ext(
    8.9,
    13.9,
    47.2,
    50.6
  )

}

################################################################################
# Farbpalette
################################################################################

palette_fun <- function(r) {

  colorNumeric(

    palette = viridis(256, option = "magma"),

    domain = values(r),

    na.color = "transparent"

  )

}

################################################################################
# Rasterwert auslesen
################################################################################

extract_value <- function(raster, lng, lat) {

  if (is.null(raster))
    return(NA)

  xy <- matrix(

    c(lng, lat),

    ncol = 2

  )

  value <- tryCatch({

    v <- terra::extract(raster, xy)

    if (is.null(v) || nrow(v) == 0)
      return(NA_real_)

    if (ncol(v) < 2)
      return(NA_real_)

    as.numeric(v[1,2])

  }, error=function(e){

    NA_real_

  })

  value

}

################################################################################
# Zahlen formatieren
################################################################################

fmt <- function(x, digits = 1, unit = "") {

  if (is.na(x))

    return("n/a")

  paste0(

    format(round(x, digits), nsmall = digits),

    unit

  )

}

################################################################################
# Alle Rasterwerte eines Punktes
################################################################################

extract_all_values <- function(lng, lat) {

  list(

    probability = extract_value(prediction, lng, lat),

    precipitation = extract_value(precipitation, lng, lat),

    temperature = extract_value(temperature, lng, lat),

    elevation = extract_value(elevation, lng, lat),

    dem = extract_value(dem, lng, lat),

    slope = extract_value(slope, lng, lat),

    aspect = extract_value(aspect, lng, lat)

  )

}

################################################################################
# Startinformationen
################################################################################

cat("\n")
cat("=========================================\n")
cat(" Random Forest Viewer\n")
cat("=========================================\n")
cat(sprintf("Raster geladen: %d\n", sum(!sapply(layers, is.null))))

for (n in names(layers)) {

  if (!is.null(layers[[n]])) {

    cat(sprintf("  ✓ %s\n", n))

  } else {

    cat(sprintf("  ✗ %s\n", n))

  }

}

cat("=========================================\n\n")
crs(prediction)
crs(precipitation)
crs(temperature)
crs(dem)
crs(slope)
crs(elevation)
cat("=========================================\n\n")
ext(prediction)
ext(precipitation)
ext(temperature)
ext(dem)
ext(slope)
ext(elevation)
cat("=========================================\n\n")

global(dem, "range", na.rm = TRUE)
global(slope, "range", na.rm = TRUE)
global(elevation, "range", na.rm = TRUE)

global(temperature, "range", na.rm = TRUE)
global(precipitation, "range", na.rm = TRUE)

cat("=========================================\n\n")
summary(values(dem))
summary(values(slope))
summary(values(elevation))