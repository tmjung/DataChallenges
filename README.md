# DataChallenges

Dieses Repository enthaelt den Code, die Eingangsdaten und zentrale Ergebnisdateien
fuer ein Data-Challenges-Projekt zur raeumlichen Modellierung archaeologischer
Fundwahrscheinlichkeiten in Bayern. Der Schwerpunkt liegt auf Iron-Age /
Eisenzeit-Fundstellen, Umweltpraediktoren und dem Vergleich mehrerer
Modellansaetze.

Die wichtigsten Modellfamilien sind:

- Random Forest mit Presence-/Absence-Punkten
- MaxEnt mit wiederholter Evaluation
- GLM als weiteres Vergleichsmodell
- Ensemble-Heatmaps aus Random Forest, MaxEnt und GLM

## Projektstruktur

```text
.
|-- data/                 Eingangsdaten, CSV/SQL und Rasterdaten
|-- output/               Versionierte Ergebnisdateien der Hauptmodelle
|-- praesentation/        Praesentationsdateien
|-- source/               Aeltere/alternative Hilfsskripte
|-- src/                  Aktive R-Skripte und Shiny-App
|-- sql_to_csv.py         Hilfsskript fuer SQL-Dump-nach-CSV-Konvertierung
|-- README.md
`-- .gitignore
```

Wichtige Unterordner:

- `src/`: alle aktiven Analyse- und Trainingsskripte
- `src/rf_viewer/`: Shiny-App zum Anzeigen der Ergebnisraster
- `data/climate/`: Klimadaten/Rasterdaten
- `data/dem/`: erwartete DEM-, Hoehen- und Slope-Raster
- `output/random_forest_5km_buffer/`: Random-Forest-Ergebnisse
- `output/maxent/` und `output/maxent_rf_grid/`: MaxEnt-Ergebnisse
- `output/glm/`: GLM-Ergebnisse
- `output/ensemble/` und `output/ensemble_glm/`: kombinierte Heatmaps und Diagnosen

## Daten

Zentrale Eingangsdaten:

- `data/ffm_vfpa_eisenzeit.csv`: Presence-Daten mit Koordinaten und
  Umweltattributen.
- `data/msqr.csv`: Muencheberger Soil Quality Rating (MSQR) als
  Bodenqualitaetsdaten.
- `data/export_*.csv`: Roh-/Exportdaten aus den zugrunde liegenden Quellen.
- `data/*.sql`: SQL-Dumps bzw. Quelldaten, aus denen CSV-Dateien erzeugt
  werden koennen.
- `data/climate/...`: WorldClim-/Klimadaten.

Mehrere Skripte erwarten ausserdem Raster unter:

```text
data/climate/precipitation.tif
data/climate/temperature.tif
data/dem/DEU_elv_msk.tif
data/dem/dem.tif
data/dem/slope.tif
```

Falls diese Dateien fehlen, muessen sie vor dem Ausfuehren der Trainingsskripte
erzeugt oder in diese Pfade gelegt werden.

## Voraussetzungen

### R

Die Analysen sind in R geschrieben. Benoetigte Pakete sind unter anderem:

```r
install.packages(c(
  "tidyverse",
  "sf",
  "terra",
  "ranger",
  "maxnet",
  "pROC",
  "ggplot2",
  "viridis",
  "rnaturalearth",
  "rnaturalearthdata",
  "shiny",
  "leaflet",
  "htmltools",
  "DBI",
  "RSQLite",
  "caret"
))
```

Einige Pakete wie `sf` und `terra` benoetigen je nach Betriebssystem
Systembibliotheken fuer Geodatenverarbeitung.

### Python

`sql_to_csv.py` nutzt nur Python-Standardbibliotheken:

- `re`
- `csv`
- `sys`

## Pfadlogik

Die wichtigsten Skripte in `src/` nutzen `src/project_paths.r`. Dadurch werden
Dateien relativ zum Projektroot gefunden, egal ob das Skript aus dem Root-Ordner
oder direkt aus `src/` gestartet wird.

Beispiele:

```r
project_path("data", "ffm_vfpa_eisenzeit.csv")
project_path("output", "glm", "glm_probability.tif")
```

## Typischer Workflow

Die Skripte koennen einzeln mit `Rscript` ausgefuehrt werden. Empfohlene
Reihenfolge fuer den aktuellen Hauptworkflow:

### 1. Random Forest trainieren

```bash
Rscript src/rf_test_4_absence_buffer.r
```

Erzeugt unter anderem:

- `output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit.tif`
- `output/random_forest_5km_buffer/random_forest_5km_buffer_fundwahrscheinlichkeit_heatmap.png`
- `output/random_forest_5km_buffer/random_forest_5km_buffer_trainingsdaten.csv`
- `output/random_forest_5km_buffer/random_forest_5km_buffer_modell.rds`

### 2. GLM trainieren

```bash
Rscript src/train_glm.r
```

Erzeugt unter anderem:

- `output/glm/glm_probability.tif`
- `output/glm/glm_probability_heatmap.png`
- `output/glm/glm_evaluation.csv`
- `output/glm/glm_coefficients.csv`

### 3. MaxEnt auf RF-Rastergrid trainieren

```bash
Rscript src/train_maxent_5.r
```

Dieses Skript nutzt bevorzugt das Random-Forest-Raster aus Schritt 1 als
Vorhersage-Grid.

Erzeugt unter anderem:

- `output/maxent_rf_grid/maxent_rf_grid_mittlere_standorteignung.tif`
- `output/maxent_rf_grid/maxent_rf_grid_unsicherheit_sd.tif`
- `output/maxent_rf_grid/maxent_rf_grid_auc_werte.csv`
- `output/maxent_rf_grid/maxent_rf_grid_prediction_grid.csv`

### 4. Ensemble aus Random Forest, MaxEnt und GLM berechnen

```bash
Rscript src/rf_me_glm_heatmap.r
```

Standardgewichte im Skript:

- Random Forest: `0.45`
- MaxEnt: `0.35`
- GLM: `0.20`

Erzeugt unter anderem:

- `output/ensemble_glm/rf_maxent_glm_weighted_mean.tif`
- `output/ensemble_glm/rf_maxent_glm_weighted_mean.png`
- `output/ensemble_glm/rf_maxent_glm_weighted_mean_presence_points.png`
- `output/ensemble_glm/rf_maxent_glm_agreement_sd.tif`
- `output/ensemble_glm/rf_maxent_glm_raster_diagnostics.csv`
- `output/ensemble_glm/rf_maxent_glm_evaluation.csv`

## Weitere wichtige Skripte

- `src/rf_test_5_msqr.r`: Random-Forest-Variante mit MSQR als zusaetzlichem
  Praediktor.
- `src/rf_roc_auc_verschiedene_werte.r`: Experimentiert mit verschiedenen
  Absence-Faktoren und Mindestdistanzen und gibt OOB Accuracy/ROC-AUC aus.
- `src/rf_me_mean_heatmap.r`: Ensemble aus Random Forest und MaxEnt ohne GLM.
- `src/train_maxent_4.r`: MaxEnt-Variante ohne externes Rastergrid.
- `src/create_tiff_terra_geodata.r` und
  `src/create_tiff_terra_geodata_less.r`: Hilfsskripte zum Erzeugen von
  Rasterdaten.
- `src/draw_points.r` und `src/draw_abs_and_pre.r`: Visualisierung von Punktdaten.
- `src/csv_reader.r` und `src/csv_creater.r`: Hilfsskripte fuer Datenimport und
  CSV-Erzeugung.

## Shiny Viewer

Der Viewer liegt in `src/rf_viewer/`.

Start:

```bash
Rscript src/rf_viewer/app.r
```

Oder innerhalb von R:

```r
shiny::runApp("src/rf_viewer")
```

Die App sucht im Projektroot nach den vorhandenen Output-Rastern und stellt sie
interaktiv dar.

## Ergebnisse und Versionierung

Der Ordner `output/` ist absichtlich versioniert. Die `.gitignore` enthaelt
globale Regeln fuer grosse Ergebnisformate wie `*.png`, `*.tif` und `*.rds`,
aber hebt diese Regeln fuer `output/` wieder auf:

```gitignore
!output/
!output/**
```

Damit werden die zentralen Modelloutputs im Repository mitgefuehrt, waehrend
andere generierte Raster/Bilder ausserhalb von `output/` weiterhin ignoriert
werden koennen.

## Reproduzierbarkeit

Viele Skripte setzen feste Seeds, zum Beispiel:

```r
set.seed(42)
```

Trotzdem koennen Ergebnisse leicht abweichen, wenn:

- andere Paketversionen installiert sind,
- externe Geodaten von `rnaturalearth` aktualisiert werden,
- Rasterdaten mit anderer Aufloesung oder Projektion verwendet werden,
- Zufallsstichproben in R-Versionen unterschiedlich behandelt werden.

Fuer reproduzierbare Vergleiche sollten die Eingangsdaten, Paketversionen und
Outputdateien gemeinsam dokumentiert werden.

## Hinweise zur Ausfuehrung

- Die rechenintensiven Skripte koennen je nach Rasteraufloesung lange laufen.
- Einige Skripte laden Bayern-Grenzen ueber `rnaturalearth`; dafuer kann beim
  ersten Lauf Netzwerkzugriff erforderlich sein.
- Die Ausgaben werden groesstenteils mit `overwrite = TRUE` geschrieben.
  Bestehende Dateien im jeweiligen Outputordner koennen also ueberschrieben
  werden.
- Falls ein Skript aus `src/` heraus gestartet wird, sollte `project_paths.r`
  automatisch den Projektroot finden.

## Kurzuebersicht der wichtigsten Outputs

```text
output/random_forest_5km_buffer/
|-- random_forest_5km_buffer_fundwahrscheinlichkeit.tif
|-- random_forest_5km_buffer_fundwahrscheinlichkeit_heatmap.png
|-- random_forest_5km_buffer_fundwahrscheinlichkeit_mit_presence_punkten.png
|-- random_forest_5km_buffer_presence_absence_punkte.png
|-- random_forest_5km_buffer_trainingsdaten.csv
|-- random_forest_5km_buffer_modell.rds
`-- random_forest_variable_importance.csv

output/glm/
|-- glm_probability.tif
|-- glm_probability_heatmap.png
|-- glm_fundwahrscheinlichkeit_mit_presence_punkten.png
|-- glm_presence_absence_points.png
|-- glm_evaluation.csv
`-- glm_coefficients.csv

output/ensemble_glm/
|-- rf_maxent_glm_weighted_mean.tif
|-- rf_maxent_glm_weighted_mean.png
|-- rf_maxent_glm_weighted_mean_presence_points.png
|-- rf_maxent_glm_agreement_sd.tif
|-- rf_maxent_glm_agreement_sd.png
|-- rf_maxent_glm_raster_diagnostics.csv
`-- rf_maxent_glm_evaluation.csv
```

## Lizenz und Quellen

Eine Lizenzdatei ist aktuell nicht im Repository enthalten. Vor einer
oeffentlichen Weitergabe sollte geklaert werden:

- unter welcher Lizenz der Code stehen soll,
- ob alle Datenquellen weiterverteilt werden duerfen,
- wie die verwendeten Datenquellen in Berichten oder Praesentationen zitiert
  werden muessen.
