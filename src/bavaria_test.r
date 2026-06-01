library(sf)
library(rnaturalearth)
library(dplyr)

states_raw <- ne_download(
  scale = 10,
  type = "states",
  category = "cultural",
  returnclass = "sf"
)

#unique(states_raw$admin)

states_raw %>%
  filter(admin == "Germany") %>%
  select(name, name_en, name_de, name_uk)

bavaria_outline <- states_raw |>
  filter(name_de == "Bayern")

print("Bayern-Umriss (erste paar Zeilen):")
print(head(bavaria_outline))