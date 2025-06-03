# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT WŁODAWSKI

# Wymagane pakiety
library("terra")
library("rstac")
library("sf")
options(timeout = 1200)

# wczytanie i przygotowanie granic powiatu
setwd("C:/Users/Weronika/Desktop/zap_klasyfikacja_nadzorowana")
powiat <- vect("powiat_wlodawski.shp")
powiat_wgs84 <- project(powiat, "EPSG:4326")
powiat_utm <- project(powiat, "EPSG:32634")

# wyszukiwanie obrazów Sentinel-2
bbox_wgs84 <- ext(powiat_wgs84)
margin <- 0.1
bbox_buf <- c(bbox_wgs84[1] - margin, bbox_wgs84[3] - margin,
              bbox_wgs84[2] + margin, bbox_wgs84[4] + margin)

obrazy <- stac("https://earth-search.aws.element84.com/v1") |>
  stac_search(
    collections = "sentinel-2-c1-l2a",
    bbox = bbox_buf,
    datetime = "2023-05-01T00:00:00Z/2023-09-30T00:00:00Z"
  ) |>
  ext_query(`eo:cloud_cover` < 10) |>
  post_request() |>
  items_fetch()

# definicja kanałów i kafli
kanaly <- c("blue", "green", "red", "nir")
tile_ids <- c("S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFB_N05.09",
              "S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFC_N05.09")

# pobieranie kafli do folderu sentinel
base_dir <- "sentinel"
dir.create(base_dir, showWarnings = FALSE)

for (tile in tile_ids) {
  tile_dir <- file.path(base_dir, tile)
  dir.create(tile_dir, showWarnings = FALSE)
  
  urls <- obrazy |>
    items_filter(properties$`s2:tile_id` == tile) |>
    assets_select(asset_names = kanaly) |>
    assets_url()
  
  for (i in seq_along(urls)) {
    download.file(urls[i], destfile = file.path(tile_dir, basename(urls[i])), mode = "wb")
  }
}

# wczytanie
rastry <- list.files("sentinel", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# mozaikowanie
kafle <- lapply(unique(dirname(rastry)), function(dir) {
  r <- rast(list.files(dir, pattern = "\\.tif$", full.names = TRUE))
  names(r) <- kanaly
  return(r)
})

# docięcie, maskowanie i ustawianie wartosci spoza zakresu na NA
r_final <- mask(crop(mosaic(sprc(kafle)), powiat_utm), powiat_utm) |>
  clamp(lower = 0, upper = 1, values = FALSE)

# Wizualizacja
plotRGB(r_final, r = 3, g = 2, b = 1, stretch = "lin")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)