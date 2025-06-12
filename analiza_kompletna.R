# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT WŁODAWSKI
# ANALIZA PORÓWNAWCZA 2018 vs 2023

# Wymagane pakiety
library("terra")
library("rstac")
library("sf")
library("caret")         # do podziału na zbiór treningowy i testowy
library("randomForest")  # klasyfikator Random Forest
library("RColorBrewer")  # ładne kolory do wizualizacji
options(timeout = 1200)

# Ustawienia folderu
setwd("D:/Nowy folder/studia/zap_projekt/zap_klasyfikacja_nadzorowana")

# 1. WCZYTANIE I PRZYGOTOWANIE GRANIC POWIATU
print("1. Wczytywanie granic powiatu...")
powiat <- vect("powiat_wlodawski.shp")
powiat_wgs84 <- project(powiat, "EPSG:4326")
powiat_utm <- project(powiat, "EPSG:32634")

# 2. PRZYGOTOWANIE BBOX DLA WYSZUKIWANIA OBRAZÓW
bbox_wgs84 <- ext(powiat_wgs84)
margin <- 0.1
bbox_buf <- c(bbox_wgs84[1] - margin, bbox_wgs84[3] - margin,
              bbox_wgs84[2] + margin, bbox_wgs84[4] + margin)

kanaly <- c("blue", "green", "red", "nir")

#======================================================================
#### 3. POBIERANIE I PRZETWARZANIE DANYCH SENTINEL-2 DLA 2023 ROKU ####
#=======================================================================

print("2. Pobieranie danych Sentinel-2 dla 2023...")

obrazy_2023 <- stac("https://earth-search.aws.element84.com/v1") |>
  stac_search(
    collections = "sentinel-2-c1-l2a",
    bbox = bbox_buf,
    datetime = "2023-05-01T00:00:00Z/2023-09-30T00:00:00Z"
  ) |>
  ext_query(`eo:cloud_cover` < 10) |>
  post_request() |>
  items_fetch()

print(paste("Znaleziono", length(obrazy_2023$features), "obrazów dla 2023 roku"))

# wyświetlenie podstawowych właściwości, tak aby można było wybrać odpowiednie kafle
if (length(obrazy_2023$features) > 0) {
  obrazy_info <- data.frame(
    nr = 1:length(obrazy_2023$features),
    tile_id = sapply(obrazy_2023$features, function(x) x$properties$`s2:tile_id`),
    data = sapply(obrazy_2023$features, function(x) substr(x$properties$datetime, 1, 10)),
    chmury = round(sapply(obrazy_2023$features, function(x) x$properties$`eo:cloud_cover`), 1),
    platforma = sapply(obrazy_2023$features, function(x) x$properties$platform)
  )
  print(obrazy_info)
}

# Pobieranie wybranych kafli 2023
tile_ids_2023 <- c("S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFB_N05.09",
                   "S2A_OPER_MSI_L2A_TL_2APS_20230928T135957_A043176_T34UFC_N05.09")

base_dir_2023 <- "sentinel_2023"
dir.create(base_dir_2023, showWarnings = FALSE)

# sprawdzenie czy wszystkie kafle są już pobrane, jeżeli tak to nie pobiera ponownie
all_tiles_downloaded_2023 <- TRUE
for (tile in tile_ids_2023) {
  tile_dir <- file.path(base_dir_2023, tile)
  if (!dir.exists(tile_dir)) {
    all_tiles_downloaded_2023 <- FALSE
    break
  }
  # sprawdzenie czy wszystkie kanały są pobrane dla tego kafla
  expected_files <- length(kanaly)
  actual_files <- length(list.files(tile_dir, pattern = "\\.tif$"))
  if (actual_files < expected_files) {
    all_tiles_downloaded_2023 <- FALSE
    break
  }
}

if (all_tiles_downloaded_2023) {
  print("Wszystkie kafle 2023 są już pobrane. Pomijam pobieranie.")
} else {
  print("Pobieranie kafli 2023...")
  for (tile in tile_ids_2023) {
    tile_dir <- file.path(base_dir_2023, tile)
    dir.create(tile_dir, showWarnings = FALSE)
    
    urls <- obrazy_2023 |>
      items_filter(properties$`s2:tile_id` == tile) |>
      assets_select(asset_names = kanaly) |>
      assets_url()
    
    for (i in seq_along(urls)) {
      file_path <- file.path(tile_dir, basename(urls[i]))
      if (!file.exists(file_path)) {
        print(paste("Pobieranie:", basename(file_path)))
        download.file(urls[i], destfile = file_path, mode = "wb")
      } else {
        print(paste("Plik już istnieje:", basename(file_path)))
      }
    }
  }
  print("Pobieranie kafli 2023 zakończone.")
}

# mozaikowanie 2023 - sprawdzenie czy plik już istnieje, jeżeli nie to mozaikuje
mozaika_2023_path <- "wyniki/mozaika_2023.tif"
dir.create("wyniki", showWarnings = FALSE)

if (file.exists(mozaika_2023_path)) {
  print("Wczytywanie istniejącej mozaiki 2023...")
  r_2023 <- rast(mozaika_2023_path)
} else {
  print("Tworzenie mozaiki 2023...")
  rastry_2023 <- list.files("sentinel_2023", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  kafle_2023 <- lapply(unique(dirname(rastry_2023)), function(dir) {
    r <- rast(list.files(dir, pattern = "\\.tif$", full.names = TRUE))
    names(r) <- kanaly
    return(r)
  })
  
  r_2023 <- mask(crop(mosaic(sprc(kafle_2023)), powiat_utm), powiat_utm) |>
    clamp(lower = 0, upper = 1, values = FALSE)
  
# Zapisanie mozaiki
  writeRaster(r_2023, mozaika_2023_path, overwrite = TRUE)
  print("Mozaika 2023 zapisana do pliku.")
}

#======================================================================
#### 3. POBIERANIE I PRZETWARZANIE DANYCH SENTINEL-2 DLA 2018 ROKU ####
#======================================================================

print("3. Pobieranie danych Sentinel-2 dla 2018...")

obrazy_2018 <- stac("https://earth-search.aws.element84.com/v1") |>
  stac_search(
    collections = "sentinel-2-c1-l2a",
    bbox = bbox_buf,
    datetime = "2018-05-01T00:00:00Z/2018-09-30T00:00:00Z"
  ) |>
  ext_query(`eo:cloud_cover` < 10) |>
  post_request() |>
  items_fetch()

print(paste("Znaleziono", length(obrazy_2018$features), "obrazów dla 2018 roku"))

# Pobieranie kafli 2018
tile_ids_2018 <- c("S2B_OPER_MSI_L2A_TL_S2RP_20230711T000257_A007641_T34UFC_N05.00",
                   "S2B_OPER_MSI_L2A_TL_S2RP_20230708T134857_A007598_T34UFB_N05.00")

base_dir_2018 <- "sentinel_2018"
dir.create(base_dir_2018, showWarnings = FALSE)

# sprawdzenie czy wszystkie kafle są już pobrane, jeżeli nie to pobiera
all_tiles_downloaded_2018 <- TRUE
for (tile in tile_ids_2018) {
  tile_dir <- file.path(base_dir_2018, tile)
  if (!dir.exists(tile_dir)) {
    all_tiles_downloaded_2018 <- FALSE
    break
  }
  # sprawdzenie czy wszystkie kanały są pobrane dla tego kafla
  expected_files <- length(kanaly)
  actual_files <- length(list.files(tile_dir, pattern = "\\.tif$"))
  if (actual_files < expected_files) {
    all_tiles_downloaded_2018 <- FALSE
    break
  }
}

if (all_tiles_downloaded_2018) {
  print("Wszystkie kafle 2018 są już pobrane. Pomijam pobieranie.")
} else {
  print("Pobieranie kafli 2018...")
  for (tile in tile_ids_2018) {
    tile_dir <- file.path(base_dir_2018, tile)
    dir.create(tile_dir, showWarnings = FALSE)
    
    urls <- obrazy_2018 |>
      items_filter(properties$`s2:tile_id` == tile) |>
      assets_select(asset_names = kanaly) |>
      assets_url()
    
    for (i in seq_along(urls)) {
      file_path <- file.path(tile_dir, basename(urls[i]))
      if (!file.exists(file_path)) {
        print(paste("Pobieranie:", basename(file_path)))
        download.file(urls[i], destfile = file_path, mode = "wb")
      } else {
        print(paste("Plik już istnieje:", basename(file_path)))
      }
    }
  }
  print("Pobieranie kafli 2018 zakończone.")
}

# mozaikowanie 2018 - sprawdzenie czy plik już istnieje, jeżeli nie to mozaikuje
mozaika_2018_path <- "wyniki/mozaika_2018.tif"

if (file.exists(mozaika_2018_path)) {
  print("Wczytywanie istniejącej mozaiki 2018...")
  r_2018 <- rast(mozaika_2018_path)
} else {
  print("Tworzenie mozaiki 2018...")
  rastry_2018 <- list.files("sentinel_2018", pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
  kafle_2018 <- lapply(unique(dirname(rastry_2018)), function(dir) {
    r <- rast(list.files(dir, pattern = "\\.tif$", full.names = TRUE))
    names(r) <- kanaly
    return(r)
  })
  
  r_2018 <- mask(crop(mosaic(sprc(kafle_2018)), powiat_utm), powiat_utm) |>
    clamp(lower = 0, upper = 1, values = FALSE)
  
  # Zapisanie mozaiki
  writeRaster(r_2018, mozaika_2018_path, overwrite = TRUE)
  print("Mozaika 2018 zapisana do pliku.")
}

#=======================================================
#### 4. PRZYGOTOWANIE DANYCH REFERENCYJNYCH (CLC) ####
#=======================================================
# dane zostały pobrane ze strony https://land.copernicus.eu/en/products/corine-land-cover

print("4. Przygotowywanie danych referencyjnych CLC")

clc <- vect("CLC18_PL/clc18_PL.shp") # już pobrane dane clc18_PL.shp, znajdujące sie w folderze CLC18_PL
clc_utm <- crop(project(clc, crs(r_2023)), powiat_utm)

# wybór klas pokrycia terenu
selected_classes <- c("112", "211", "311", "511", "231")  # zabudowa, pola, lasy, wody, pastwiska
clc_sel <- clc_utm[clc_utm$CODE_18 %in% selected_classes, ]
clc_sel$class <- as.factor(clc_sel$CODE_18)

# generowanie punktów treningowych
set.seed(42)
points_train <- spatSample(clc_sel, size = 600, method = "random")
points_train$class <- as.factor(points_train$CODE_18)

#=============================================================================
#### 5. KLASYFIKACJA POKRYCIA TERENU DLA 2023 ROKU ####
#=============================================================================
print("5. Klasyfikacja pokrycia terenu dla 2023...")

# ekstrakcja wartości spektralnych dla punktów treningowych
train_data_2023 <- extract(r_2023, points_train, df = TRUE)
train_data_2023$class <- points_train$class
train_data_2023 <- na.omit(train_data_2023)

# podział na zbiór treningowy i testowy
set.seed(123)
split_index <- createDataPartition(train_data_2023$class, p = 0.7, list = FALSE)
train_set_2023 <- train_data_2023[split_index, ]
test_set_2023 <- train_data_2023[-split_index, ]

# sprawdzenie czy model i klasyfikacja już istnieją
model_2023_path <- "wyniki/model_rf_2023.rds"
classified_2023_path <- "wyniki/klasyfikacja_2023.tif"

if (file.exists(model_2023_path) && file.exists(classified_2023_path)) {
  print("Wczytywanie istniejącego modelu i klasyfikacji 2023...")
  model_rf_2023 <- readRDS(model_2023_path)
  classified_2023 <- rast(classified_2023_path)
} else {
  print("Trening modelu Random Forest dla 2023...")
  # trening modelu Random Forest dla 2023
  model_rf_2023 <- randomForest(as.factor(class) ~ blue + green + red + nir,
                                data = train_set_2023,
                                ntree = 200,
                                importance = TRUE)
  
  # zapisanie modelu
  saveRDS(model_rf_2023, model_2023_path)
  
  # klasyfikacja całego obrazu 2023
  print("Klasyfikacja obrazu 2023...")
  classified_2023 <- predict(r_2023, model_rf_2023, type = "response")
  
  # zapisanie klasyfikacji
  writeRaster(classified_2023, classified_2023_path, overwrite = TRUE)
  print("Model i klasyfikacja 2023 zapisane do plików.")
}

# ocena skuteczności dla 2023
pred_test_2023 <- predict(model_rf_2023, test_set_2023)
conf_matrix_2023 <- confusionMatrix(pred_test_2023, test_set_2023$class)

print("WYNIKI KLASYFIKACJI 2023")
print(conf_matrix_2023)

# macierz pomyłek i najczęstsze błędy
print(conf_matrix_2023$table)
print(conf_matrix_2023$byClass[, c("Precision", "Recall", "F1")])

#=============================================================================
# 6. KLASYFIKACJA POKRYCIA TERENU DLA 2018 ROKU
#=============================================================================
print("6. Klasyfikacja pokrycia terenu dla 2018...")

# ekstrakcja wartości spektralnych dla punktów treningowych (te same punkty!)
train_data_2018 <- extract(r_2018, points_train, df = TRUE)
train_data_2018$class <- points_train$class
train_data_2018 <- na.omit(train_data_2018)

# podział na zbiór treningowy i testowy (ten sam podział co dla 2023)
train_set_2018 <- train_data_2018[split_index, ]
test_set_2018 <- train_data_2018[-split_index, ]

# sprawdzenie czy model i klasyfikacja już istnieją
model_2018_path <- "wyniki/model_rf_2018.rds"
classified_2018_path <- "wyniki/klasyfikacja_2018.tif"

if (file.exists(model_2018_path) && file.exists(classified_2018_path)) {
  print("Wczytywanie istniejącego modelu i klasyfikacji 2018...")
  model_rf_2018 <- readRDS(model_2018_path)
  classified_2018 <- rast(classified_2018_path)
} else {
  print("Trening modelu Random Forest dla 2018...")
  # trening modelu Random Forest dla 2018
  model_rf_2018 <- randomForest(as.factor(class) ~ blue + green + red + nir,
                                data = train_set_2018,
                                ntree = 200,
                                importance = TRUE)
  
  # zapisanie modelu
  saveRDS(model_rf_2018, model_2018_path)
  
  # klasyfikacja całego obrazu 2018
  print("Klasyfikacja obrazu 2018...")
  classified_2018 <- predict(r_2018, model_rf_2018, type = "response")
  
  # zapisanie klasyfikacji
  writeRaster(classified_2018, classified_2018_path, overwrite = TRUE)
  print("Model i klasyfikacja 2018 zapisane do plików.")
}

# ocena skuteczności dla 2018
pred_test_2018 <- predict(model_rf_2018, test_set_2018)
conf_matrix_2018 <- confusionMatrix(pred_test_2018, test_set_2018$class)

print("WYNIKI KLASYFIKACJI 2018")
print(conf_matrix_2018)

# macierz pomyłek i najczęstsze błędy
print(conf_matrix_2018$table)
print(conf_matrix_2018$byClass[, c("Precision", "Recall", "F1")])

#=============================================================================
# 7. ANALIZA ZMIAN KLASYFIKACJI 2018 vs 2023
#=============================================================================
print("7. Analiza zmian klasyfikacji...")

# tworzenie mapy zmian klasyfikacji
change_matrix <- cbind(values(classified_2018), values(classified_2023))
change_matrix <- change_matrix[complete.cases(change_matrix), ]

# macierz przejść między klasami
transition_matrix <- table(
  "2018" = change_matrix[,1], 
  "2023" = change_matrix[,2]
)

print("MACIERZ PRZEJŚĆ MIĘDZY KLASAMI")
print(transition_matrix)

# obliczanie zmian powierzchni dla każdej klasy
freq_2018 <- freq(classified_2018)
freq_2023 <- freq(classified_2023)

change_stats <- data.frame(
  class = selected_classes,
  area_2018 = freq_2018[match(selected_classes, freq_2018[,"value"]), "count"],
  area_2023 = freq_2023[match(selected_classes, freq_2023[,"value"]), "count"]
)
change_stats$area_2018[is.na(change_stats$area_2018)] <- 0
change_stats$area_2023[is.na(change_stats$area_2023)] <- 0
change_stats$change_pixels <- change_stats$area_2023 - change_stats$area_2018
change_stats$change_percent <- round(100 * change_stats$change_pixels / change_stats$area_2018, 2)

print("ZMIANY POWIERZCHNI KLAS POKRYCIA TERENU")
print(change_stats)

#=============================================================================
### 8. ANALIZA NDVI I ZMIAN NDVI ####
#=============================================================================
print("8. Analiza NDVI")

# sprawdzenie czy rastry NDVI już istnieją
ndvi_2018_path <- "wyniki/ndvi_2018.tif"
ndvi_2023_path <- "wyniki/ndvi_2023.tif"
ndvi_diff_path <- "wyniki/ndvi_zmiana.tif"

if (file.exists(ndvi_2018_path) && file.exists(ndvi_2023_path) && file.exists(ndvi_diff_path)) {
  print("Wczytywanie istniejących rastrów NDVI")
  ndvi_2018 <- rast(ndvi_2018_path)
  ndvi_2023 <- rast(ndvi_2023_path)
  ndvi_diff <- rast(ndvi_diff_path)
} else {
  print("Obliczanie NDVI dla obu lat...")
  # obliczenie NDVI dla 2018 i 2023
  ndvi_2018 <- (r_2018$nir - r_2018$red) / (r_2018$nir + r_2018$red)
  ndvi_2023 <- (r_2023$nir - r_2023$red) / (r_2023$nir + r_2023$red)
  
  # różnica NDVI
  ndvi_diff <- ndvi_2023 - ndvi_2018
  
  # zapisanie rastrów NDVI
  writeRaster(ndvi_2018, ndvi_2018_path, overwrite = TRUE)
  writeRaster(ndvi_2023, ndvi_2023_path, overwrite = TRUE)
  writeRaster(ndvi_diff, ndvi_diff_path, overwrite = TRUE)
  print("Rastry NDVI zapisane do plików.")
}

# ekstrakcja wartości NDVI dla całego powiatu
valores_2023 <- extract(r_2023, powiat_utm)
valores_2018 <- extract(r_2018, powiat_utm)

# czyszczenie danych
valores_2023_vals <- valores_2023[, -1]
valores_2018_vals <- valores_2018[, -1]
valid_rows <- complete.cases(valores_2023_vals) & complete.cases(valores_2018_vals)
valores_2023_clean <- valores_2023_vals[valid_rows, ]
valores_2018_clean <- valores_2018_vals[valid_rows, ]

# obliczenie NDVI dla pikseli
ndvi_2023_clean <- (valores_2023_clean$nir - valores_2023_clean$red) / 
  (valores_2023_clean$nir + valores_2023_clean$red)
ndvi_2018_clean <- (valores_2018_clean$nir - valores_2018_clean$red) / 
  (valores_2018_clean$nir + valores_2018_clean$red)
ndvi_diff_clean <- ndvi_2023_clean - ndvi_2018_clean

print("STATYSTYKI NDVI")
print("NDVI 2018:")
print(summary(ndvi_2018_clean))
print("NDVI 2023:")
print(summary(ndvi_2023_clean))
print("Różnica NDVI:")
print(summary(ndvi_diff_clean))

#=============================================================================
### 9. WIZUALIZACJE ####
# =============================================================================
print("9. Tworzenie wizualizacji")

# kolory dla klas
colors <- brewer.pal(length(selected_classes), "Set1")
class_labels <- c("Zabudowa", "Pola uprawne", "Lasy", "Wody", "Pastwiska")

# WIZUALIZACJA 1: Porównanie obrazów RGB
par(mfrow = c(1, 2))
plotRGB(r_2018, r = 3, g = 2, b = 1, stretch = "lin", main = "Sentinel-2 RGB 2018")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

plotRGB(r_2023, r = 3, g = 2, b = 1, stretch = "lin", main = "Sentinel-2 RGB 2023")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

# WIZUALIZACJA 2: Porównanie klasyfikacji
plot(classified_2018, col = colors, main = "Klasyfikacja pokrycia terenu 2018")
legend("topright", legend = class_labels, fill = colors, title = "Klasy", cex = 0.8)

plot(classified_2023, col = colors, main = "Klasyfikacja pokrycia terenu 2023")
legend("topright", legend = class_labels, fill = colors, title = "Klasy", cex = 0.8)


# WIZUALIZACJA 3: porównanie NDVI
plot(ndvi_2018, main = "NDVI 2018", col = hcl.colors(100, "Greens"))
plot(ndvi_2023, main = "NDVI 2023", col = hcl.colors(100, "Greens"))

# różnica NDVI
par(mfrow = c(1, 1))
plot(ndvi_diff, main = "Zmiana NDVI (2023-2018)", col = hcl.colors(100, "RdBu", rev = FALSE))

# Histogram zmian NDVI
hist(ndvi_diff_clean, breaks = 50, main = "Histogram zmian NDVI", 
     xlab = "Różnica NDVI", col = "lightblue")
abline(v = 0, col = "red", lwd = 2, lty = 2)

#=============================================================================
# 10. PORÓWNANIE Z BDOT10K
#=============================================================================
print("10. Porównanie z danymi BDOT10k")

# Wczytanie warstw BDOT10k
ptlz <- vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTLZ_A.shp")   # lasy
ptwp <- vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTWP_A.shp")   # wody
pttr <- vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTTR_A.shp")   # pastwiska
ptzb <- vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTZB_A.shp")   # zabudowa

# Przekształcenie do układu współrzędnych
crs_utm <- crs(r_2023)
ptlz <- project(ptlz, crs_utm)
ptwp <- project(ptwp, crs_utm)
pttr <- project(pttr, crs_utm)
ptzb <- project(ptzb, crs_utm)

# Dodanie kodów klas
ptlz$KOD_KLASY <- 1  # las
ptwp$KOD_KLASY <- 2  # wody
pttr$KOD_KLASY <- 3  # trawiasta
ptzb$KOD_KLASY <- 4  # zabudowa

# Połączenie warstw
bdot_all <- rbind(ptlz, ptwp, pttr, ptzb)

# Sprawdzenie czy raster BDOT już istnieje
bdot_raster_path <- "wyniki/bdot_raster.tif"

if (file.exists(bdot_raster_path)) {
  print("Wczytywanie istniejącego rastra BDOT...")
  bdot_crop <- rast(bdot_raster_path)
} else {
  print("Rasteryzacja danych BDOT10k...")
  # Rasteryzacja
  bdot_raster <- rasterize(bdot_all, r_2023[[1]], field = "KOD_KLASY")
  bdot_crop <- mask(crop(bdot_raster, powiat_utm), powiat_utm)
  
  # Zapisanie rastra BDOT
  writeRaster(bdot_crop, bdot_raster_path, overwrite = TRUE)
  print("Raster BDOT zapisany do pliku.")
}

# Porównanie z klasyfikacją 2023
compare_2023 <- cbind(values(classified_2023), values(bdot_crop))
compare_2023 <- compare_2023[complete.cases(compare_2023), ]
conf_matrix_bdot_2023 <- table(klasyfikacja = compare_2023[,1], bdot = compare_2023[,2])

print("PORÓWNANIE KLASYFIKACJI 2023 Z BDOT10K")
print(conf_matrix_bdot_2023)

# =============================================================================
# 11 PORÓWNANIE ZGODNOŚCI KLASYFIKACJI Z CLC DLA OBU LAT
# =============================================================================
print("11. PORÓWNANIE KLASYFIKACJI Z CLC")

# Przygotowanie rastra CLC
clc_raster_path <- "wyniki/clc_raster.tif"
if (file.exists(clc_raster_path)) {
  clc_raster_crop <- rast(clc_raster_path)
} else {
  clc_raster_crop <- mask(crop(rasterize(clc_sel, r_2023[[1]], field = "CODE_18"), powiat_utm), powiat_utm)
  writeRaster(clc_raster_crop, clc_raster_path, overwrite = TRUE)
}

# funkcja do obliczania zgodności
calculate_agreement <- function(classified, clc_raster, year) {
  compare_data <- cbind(values(classified), values(clc_raster))
  compare_data <- compare_data[complete.cases(compare_data), ]
  
  if (nrow(compare_data) > 0) {
    agreement <- sum(compare_data[,1] == compare_data[,2])
    total <- nrow(compare_data)
    accuracy <- round(100 * agreement / total, 2)
    
    print(paste("Zgodność", year, "z CLC:", accuracy, "% (", agreement, "/", total, "pikseli)"))
    return(accuracy)
  }
  return(NA)
}

# obliczenie zgodności dla obu lat
accuracy_2018_clc <- calculate_agreement(classified_2018, clc_raster_crop, "2018")
accuracy_2023_clc <- calculate_agreement(classified_2023, clc_raster_crop, "2023")

# porównanie zgodności
if (!is.na(accuracy_2018_clc) && !is.na(accuracy_2023_clc)) {
  diff_clc <- accuracy_2023_clc - accuracy_2018_clc
  print(paste("Różnica zgodności (2023-2018):", round(diff_clc, 2), "p.p."))
  if (diff_clc > 0) print("✅ Poprawa jakości klasyfikacji")
  else if (diff_clc < 0) print("⚠️ Pogorszenie jakości klasyfikacji")
  else print("➡️ Brak zmian w jakości klasyfikacji")
}

# WIZUALIZACJA PORÓWNANIA Z BDOT10K
print("Tworzenie wizualizacji porównania z BDOT10k...")

# przygotowanie BDOT do wizualizacji
bdot_viz <- bdot_crop
bdot_viz <- as.factor(bdot_viz)
levels(bdot_viz) <- data.frame(value = 1:4,
                               label = c("Lasy", "Wody", "Pola i pastwiska", "Zabudowa"))

# Kolory dla BDOT (te same co dla klasyfikacji)
kolory_bdot <- c("darkgreen", "blue", "gold", "red")

# wizualizacja porównania
par(mfrow = c(1, 2))

# 1. Klasyfikacja 2023
plot(classified_2023, col = colors, main = "Klasyfikacja 2023 (Random Forest)")
legend("topright", legend = class_labels, fill = colors, title = "Klasy RF", cex = 0.7)

# 2. BDOT10k
plot(bdot_viz, col = kolory_bdot, main = "Pokrycie terenu BDOT10k")
legend("topright", legend = c("Lasy", "Wody", "Pola i pastwiska", "Zabudowa"), 
       fill = kolory_bdot, title = "Klasy BDOT", cex = 0.7)

# 3. Mapa różnic (gdzie klasyfikacja się różni od BDOT)
# Tworzenie mapy zgodności/różnic
classification_vals <- values(classified_2023)
bdot_vals <- values(bdot_crop)

# Mapowanie klas do wspólnych kategorii (uproszczone)
# RF: 112->1(zabudowa), 211->2(pola), 311->3(lasy), 511->4(wody), 231->2(pastwiska)
# BDOT: 1->3(lasy), 2->4(wody), 3->2(pola/pastwiska), 4->1(zabudowa)
rf_mapped <- classification_vals
rf_mapped[classification_vals == 112] <- 1  # zabudowa
rf_mapped[classification_vals == 211] <- 2  # pola
rf_mapped[classification_vals == 311] <- 3  # lasy  
rf_mapped[classification_vals == 511] <- 4  # wody
rf_mapped[classification_vals == 231] <- 2  # pastwiska -> pola

bdot_mapped <- bdot_vals
bdot_mapped[bdot_vals == 1] <- 3  # lasy
bdot_mapped[bdot_vals == 2] <- 4  # wody
bdot_mapped[bdot_vals == 3] <- 2  # pola/pastwiska
bdot_mapped[bdot_vals == 4] <- 1  # zabudowa

# mapa zgodności (1 = zgodne, 0 = różne)
agreement_vals <- ifelse(rf_mapped == bdot_mapped, 1, 0)
agreement_raster <- classified_2023
values(agreement_raster) <- agreement_vals

par(mfrow = c(1, 1))
plot(agreement_raster, col = c("red", "green"), main = "Zgodność klasyfikacji z BDOT")
legend("topright", legend = c("Różnice", "Zgodność"), fill = c("red", "green"), 
       title = "Porównanie", cex = 0.7)


# obliczenie procentu zgodności
agreement_percent <- round(100 * sum(agreement_vals == 1, na.rm = TRUE) / sum(!is.na(agreement_vals)), 2)
print(paste("Ogólna zgodność klasyfikacji z BDOT10k:", agreement_percent, "%"))

# =============================================================================
# 12. PODSUMOWANIE
# =============================================================================

print("PODSUMOWANIE ANALIZY")
print("Pobrano i przetworzono dane Sentinel-2 dla 2018 i 2023")
print("Wykonano klasyfikację pokrycia terenu dla obu lat")
print("Przeanalizowano zmiany między latami")
print("Obliczono NDVI i jego zmiany")
print("Zidentyfikowano hotspoty zmian")
print("Porównano z danymi referencyjnymi BDOT10k")
print("Zapisano wszystkie wyniki w folderze 'wyniki/'")

print("GŁÓWNE WNIOSKI")
cat("Dokładność klasyfikacji 2018:", round(conf_matrix_2018$overall["Accuracy"], 3), "\n")
cat("Dokładność klasyfikacji 2023:", round(conf_matrix_2023$overall["Accuracy"], 3), "\n")
cat("Zgodność z BDOT10k:", agreement_percent, "%\n")
if (exists("accuracy_2018_clc") && exists("accuracy_2023_clc")) {
  cat("Zgodność 2018 z CLC:", accuracy_2018_clc, "%\n")
  cat("Zgodność 2023 z CLC:", accuracy_2023_clc, "%\n")
}
cat("Średnia zmiana NDVI:", round(mean(ndvi_diff_clean), 4), "\n")
cat("Piksele ze wzrostem NDVI:", sum(ndvi_diff_clean > 0), "(", round(100*sum(ndvi_diff_clean > 0)/length(ndvi_diff_clean), 1), "%)\n")
cat("Piksele ze spadkiem NDVI:", sum(ndvi_diff_clean < 0), "(", round(100*sum(ndvi_diff_clean < 0)/length(ndvi_diff_clean), 1), "%)\n")
cat("Hotspoty wzrostu NDVI (>", prog_wzrost, "):", piksele_wzrost, "pikseli\n")
cat("Hotspoty spadku NDVI (<", prog_spadek, "):", piksele_spadek, "pikseli\n")

# 1 = błędny w obu, 2 = poprawny tylko 2018, 3 = poprawny tylko 2023, 4 = poprawny w obu
change_map <- rep(NA, length(correct_2018))
idx <- which(!is.na(correct_2018) & !is.na(correct_2023))
change_map[idx] <- 1 + correct_2018[idx] + 2 * correct_2023[idx]

change_raster <- classified_2018
values(change_raster) <- change_map
change_raster <- as.factor(change_raster)
change_labs <- c("Błędny w obu", "Poprawny tylko 2018", "Poprawny tylko 2023", "Poprawny w obu")
levels(change_raster) <- data.frame(ID=1:4, Poprawnosc=change_labs)
change_cols <- c("grey80", "orange", "blue", "green")

par(mfrow = c(1,2))
plot(change_raster, col = change_cols, main = "Zmiana poprawności klasyfikacji", axes=FALSE, legend=TRUE)
plot(ndvi_diff, main = "Zmiana NDVI (2023-2018)", col = hcl.colors(100, "RdBu", rev = FALSE))
