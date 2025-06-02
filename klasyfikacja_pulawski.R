# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT PUŁAWSKI


# Część 1: Pobieranie i przycinanie danych Sentinel-2

required_packages <- c("terra", "rstac", "sf")
#install.packages("terra", "rstac", "sf")
library("terra")
library("rstac")
library("sf")
options(timeout = 1200)

# Wczytanie granic powiatu puławskiego
setwd("C:/Users/Weronika/Desktop/zap_klasyfikacja_nadzorowana")
powiat <- vect("powiat_pulawski.shp")

# Przekonwertowanie układu współrzędnych na WGS84
powiat_sf <- st_as_sf(powiat)
powiat_wgs84_sf <- st_transform(powiat_sf, 4326)
powiat_wgs84 <- vect(powiat_wgs84_sf)

# Pobranie bounding box w WGS84
bbox_wgs84 <- ext(powiat_wgs84)
bbox_coords <- c(bbox_wgs84[1], bbox_wgs84[3], bbox_wgs84[2], bbox_wgs84[4])

# Margines bo bez marginesu nie chce się pobrać
margin <- 0.5
bbox_coords_buf <- c(bbox_coords[1] - margin,
                     bbox_coords[2] - margin,
                     bbox_coords[3] + margin,
                     bbox_coords[4] + margin)

# Wyszukiwanie scen Sentinel-2
stac_source = stac("https://earth-search.aws.element84.com/v1")
stac_source |>
  stac_search(
    collections = "sentinel-2-c1-l2a",
    bbox = bbox_coords_buf,
    datetime = "2024-03-01T00:00:00Z/2024-10-31T00:00:00Z") |>
  ext_query(`eo:cloud_cover` < 20) |>
  post_request() -> obrazy

print(paste("Znaleziono", length(obrazy$features), "scen"))



for(i in 1:length(obrazy$features)) {
  scena <- obrazy$features[[i]]
  
# info o scenie - żeby wybrać kafelek
  kafelek <- scena$properties$`mgrs:utm_zone`
  latitude_band <- scena$properties$`mgrs:latitude_band`
  grid_square <- scena$properties$`mgrs:grid_square`
  kafelek_pelny <- paste0(kafelek, latitude_band, grid_square)
  
  data <- substr(scena$properties$datetime, 1, 10)
  zachmurzenie <- round(scena$properties$`eo:cloud_cover`, 2)
  id <- scena$id
  
  print(paste("Scena", i, "- Kafelek:", kafelek_pelny))
  print(paste("  Data:", data))
  print(paste("  Zachmurzenie:", zachmurzenie, "%"))
  print(paste("  ID:", id))
  print("---")
}


kanaly = c("blue", "green", "red", "nir")

# 34UEB
obrazy |>
  items_filter(id == "S2A_T34UEB_20241018T095240_L2A") |>
  assets_download(asset_names = kanaly, overwrite = TRUE)

# 34UFB
obrazy |>
  items_filter(id == "S2B_T34UFB_20241020T094031_L2A") |>
  assets_download(asset_names = kanaly, overwrite = TRUE)

# 34UFC
obrazy |>
  items_filter(id == "S2B_T34UFC_20241017T093032_L2A") |>
  assets_download(asset_names = kanaly, overwrite = TRUE)

# 34UEC
obrazy |>
  items_filter(id == "S2A_T34UEC_20241018T095240_L2A") |>
  assets_download(asset_names = kanaly, overwrite = TRUE)


# Wczytanie pobranych rastrów
rastry_sentinel <- list.files("sentinel-2-c1-l2a", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

