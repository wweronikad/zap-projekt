# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT WŁODAWSKI

# Wymagane pakiety
library("terra")
library("rstac")
library("sf")
options(timeout = 1200)

# wczytanie i przygotowanie granic powiatu włodawskiego
setwd("D:/Nowy folder/studia/zap_projekt/zap_klasyfikacja_nadzorowana")
powiat <- vect("powiat_wlodawski.shp")
powiat_wgs84 <- project(powiat, "EPSG:4326")
powiat_utm <- project(powiat, "EPSG:32634")

############## wyszukiwanie obrazów Sentinel-2 z 2023 roku do analizy ###################
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

print(paste("Znaleziono", length(obrazy$features), "obrazów dla 2023 roku"))

# Wyświetlenie podstawowych właściwości
if (length(obrazy$features) > 0) {
  obrazy_info <- data.frame(
    nr = 1:length(obrazy$features),
    tile_id = sapply(obrazy$features, function(x) x$properties$`s2:tile_id`),
    data = sapply(obrazy$features, function(x) substr(x$properties$datetime, 1, 10)),
    chmury = round(sapply(obrazy$features, function(x) x$properties$`eo:cloud_cover`), 1),
    platforma = sapply(obrazy$features, function(x) x$properties$platform)
  )
  print(obrazy_info)
}

# definicja kanałów i kafli
kanaly <- c("blue", "green", "red", "nir")
tile_ids <- c("S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFB_N05.09",
              "S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFC_N05.09")

# pobieranie kafli do folderu sentinel_2023
base_dir <- "sentinel_2023"
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
rastry <- list.files("sentinel_2023", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# mozaikowanie
kafle <- lapply(unique(dirname(rastry)), function(dir) {
  r <- rast(list.files(dir, pattern = "\\.tif$", full.names = TRUE))
  names(r) <- kanaly
  return(r)
})

# docięcie, maskowanie i ustawianie wartosci spoza zakresu na NA
r_final <- mask(crop(mosaic(sprc(kafle)), powiat_utm), powiat_utm) |>
  clamp(lower = 0, upper = 1, values = FALSE)

r_2023 <- r_final

# Wizualizacja
plotRGB(r_final, r = 3, g = 2, b = 1, stretch = "lin")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)



########### wyszukiwanie obrazów Sentinel-2 z 2018 roku ##############

obrazy2 <- stac("https://earth-search.aws.element84.com/v1") |>
  stac_search(
    collections = "sentinel-2-c1-l2a",
    bbox = bbox_buf,
    datetime = "2018-05-01T00:00:00Z/2018-09-30T00:00:00Z"
  ) |>
  ext_query(`eo:cloud_cover` < 10) |>
  post_request() |>
  items_fetch()

# Sprawdzenie czy znaleziono jakieś obrazy
print(paste("Znaleziono", length(obrazy2$features), "obrazów dla 2018 roku"))

# Wyświetlenie podstawowych właściwości
if (length(obrazy2$features) > 0) {
  obrazy_info <- data.frame(
    nr = 1:length(obrazy2$features),
    tile_id = sapply(obrazy2$features, function(x) x$properties$`s2:tile_id`),
    data = sapply(obrazy2$features, function(x) substr(x$properties$datetime, 1, 10)),
    chmury = round(sapply(obrazy2$features, function(x) x$properties$`eo:cloud_cover`), 1),
    platforma = sapply(obrazy2$features, function(x) x$properties$platform)
  )
  print(obrazy_info)
}

# definicja kanałów i kafli
kanaly <- c("blue", "green", "red", "nir")
tile_ids <- c("S2B_OPER_MSI_L2A_TL_S2RP_20230711T000257_A007641_T34UFC_N05.00",
              "S2B_OPER_MSI_L2A_TL_S2RP_20230708T134857_A007598_T34UFB_N05.00")

# pobieranie kafli do folderu sentinel_2018
base_dir <- "sentinel_2018"
dir.create(base_dir, showWarnings = FALSE)

for (tile in tile_ids) {
  tile_dir <- file.path(base_dir, tile)
  dir.create(tile_dir, showWarnings = FALSE)
  
  urls <- obrazy2 |>
    items_filter(properties$`s2:tile_id` == tile) |>
    assets_select(asset_names = kanaly) |>
    assets_url()
  
  for (i in seq_along(urls)) {
    download.file(urls[i], destfile = file.path(tile_dir, basename(urls[i])), mode = "wb")
  }
}

# wczytanie
rastry <- list.files("sentinel_2018", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)

# mozaikowanie
kafle <- lapply(unique(dirname(rastry)), function(dir) {
  r <- rast(list.files(dir, pattern = "\\.tif$", full.names = TRUE))
  names(r) <- kanaly
  return(r)
})

# docięcie, maskowanie i ustawianie wartosci spoza zakresu na NA
r_final <- mask(crop(mosaic(sprc(kafle)), powiat_utm), powiat_utm) |>
  clamp(lower = 0, upper = 1, values = FALSE)

r_2018 <- r_final

par(mfrow = c(1, 2))
plotRGB(r_2023, r = 3, g = 2, b = 1, stretch = "lin", main = "2023")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

plotRGB(r_2018, r = 3, g = 2, b = 1, stretch = "lin", main = "2018") 
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

par(mfrow = c(1, 1))


########### OBLICZANIE NVDI DLA OBSZARU POWIATU ######################

# granice analizy z powiatu
valores_2023 <- extract(r_2023, powiat_utm)
valores_2018 <- extract(r_2018, powiat_utm)

print(paste("Liczba wyciągniętych pikseli (2023):", nrow(valores_2023)))
print(paste("Liczba wyciągniętych pikseli (2018):", nrow(valores_2018)))

# usunięcie kolumny ID (bo był z nią problem)
valores_2023_vals <- valores_2023[, -1]
valores_2018_vals <- valores_2018[, -1]

# czyszczenie pikseli z NA
valid_rows <- complete.cases(valores_2023_vals) & complete.cases(valores_2018_vals)
valores_2023_clean <- valores_2023_vals[valid_rows, ]
valores_2018_clean <- valores_2018_vals[valid_rows, ]

print(paste("Liczba prawidłowych pikseli po czyszczeniu:", nrow(valores_2023_clean)))

# obliczenie NDVI
ndvi_2023_clean <- (valores_2023_clean$nir - valores_2023_clean$red) / 
  (valores_2023_clean$nir + valores_2023_clean$red)
ndvi_2018_clean <- (valores_2018_clean$nir - valores_2018_clean$red) / 
  (valores_2018_clean$nir + valores_2018_clean$red)

# różnica NDVI
ndvi_diff_clean <- ndvi_2023_clean - ndvi_2018_clean

print("NDVI 2023:")
print(summary(ndvi_2023_clean))

print("NDVI 2018:")
print(summary(ndvi_2018_clean))

print("Różnica NDVI:")
print(summary(ndvi_diff_clean))




# WIZUALIZACJA RÓŻNIC NDVI

# Histogram różnic NDVI
par(mfrow = c(1, 1))

# Wybierz losowo tylko 1% pikseli do wizualizacji (ok.125k punktów)
set.seed(123)  # powtarzalnosc
sample_size <- round(length(ndvi_diff_clean) * 0.01)  # 1% danych
sample_indices <- sample(1:length(ndvi_diff_clean), sample_size)

# Próbkowane dane
ndvi_diff_sample <- ndvi_diff_clean[sample_indices]
ndvi_2023_sample <- ndvi_2023_clean[sample_indices]
ndvi_2018_sample <- ndvi_2018_clean[sample_indices]

print(paste("Użyto", sample_size, "pikseli z", length(ndvi_diff_clean), "dostępnych"))

# WYKRESY
par(mfrow = c(1, 2))

# 1. Histogram
hist(ndvi_diff_sample, breaks = 50, 
     main = "Histogram różnic NDVI (próbka 1%)", 
     xlab = "Różnica NDVI", 
     col = "lightblue")
abline(v = 0, col = "red", lwd = 2, lty = 2)
abline(v = mean(ndvi_diff_sample), col = "green", lwd = 2)

# 3. Density plot
plot(density(ndvi_diff_sample), 
     main = "Rozkład różnic NDVI",
     col = "darkgreen", lwd = 2)
polygon(density(ndvi_diff_sample), col = "lightgreen")
abline(v = 0, col = "red", lwd = 2, lty = 2)

par(mfrow = c(1, 1))


# Dodatkowe analizy
print("ANALIZA KIERUNKU ZMIAN NDVI")
wzrost <- sum(ndvi_diff_clean > 0)
spadek <- sum(ndvi_diff_clean < 0)
brak_zmian <- sum(ndvi_diff_clean == 0)

print(paste("Piksele ze wzrostem NDVI:", wzrost, "(",round(wzrost/length(ndvi_diff_clean)*100,1),"%)"))
print(paste("Piksele ze spadkiem NDVI:", spadek, "(",round(spadek/length(ndvi_diff_clean)*100,1),"%)"))
print(paste("Piksele bez zmian:", brak_zmian, "(",round(brak_zmian/length(ndvi_diff_clean)*100,1),"%)"))

# Percentyle różnic
print("PERCENTYLE RÓŻNIC NDVI")
percentyle <- quantile(ndvi_diff_clean, probs = c(0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99))
print(round(percentyle, 4))

# Znaczące zmiany
print("ZNACZĄCE ZMIANY")
duzy_wzrost <- sum(ndvi_diff_clean > 0.2)
duzy_spadek <- sum(ndvi_diff_clean < -0.2)
print(paste("Duży wzrost NDVI (>0.2):", duzy_wzrost, "pikseli"))
print(paste("Duży spadek NDVI (<-0.2):", duzy_spadek, "pikseli"))

############# MAPA ZMIAN NDVI ############################################

ndvi_2018 <- (r_2018$nir - r_2018$red) / (r_2018$nir + r_2018$red)
ndvi_2023 <- (r_2023$nir - r_2023$red) / (r_2023$nir + r_2023$red)

par(mfrow = c(1, 1))
plot(ndvi_2023 - ndvi_2018, main = "Zmiana", col = hcl.colors(100, "RdBu", rev = FALSE))
