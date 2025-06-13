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


# WIZUALIZACJA 1: Porównanie obrazów RGB
par(mfrow = c(1, 2))
plotRGB(r_2018, r = 3, g = 2, b = 1, stretch = "lin", main = "Sentinel-2 RGB 2018")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

plotRGB(r_2023, r = 3, g = 2, b = 1, stretch = "lin", main = "Sentinel-2 RGB 2023")
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)

#RP
############TRENING NA PODSTAWIE CLC I PORÓWNANIE Z BDOT10K#########################



# 1. WCZYTANIE I PRZYGOTOWANIE DANYCH REFERENCYJNYCH (np. CLC)
# Zakładamy, że warstwa CLC jest już pobrana i dotyczy całej Polski lub regionu
clc = vect("CLC18_PL/clc18_PL.shp")  # zamień ścieżkę na odpowiednią
clc_utm = crop(project(clc, crs(r_2023)), powiat_utm)  # przekształcenie i docięcie do powiatu

# 2. WYBÓR KILKU KLAS I PRZYGOTOWANIE PRÓBEK
# Wybieramy tylko kilka klas (np. lasy, pola uprawne, zabudowa)
# Załóżmy, że kolumna 'CODE_18' zawiera identyfikatory klas CLC
selected_classes = c("112", "211", "311", "511", "231")  # przykładowo: 112 - zabudowa, 511-wody, 211 - pola, 311 - lasy, 231- "pastwiska"
clc_sel = clc_utm[clc_utm$CODE_18 %in% selected_classes, ]
clc_sel$class = as.factor(clc_sel$CODE_18)

# 3. GENEROWANIE PUNKTÓW TRENINGOWYCH
set.seed(42)
points_train = spatSample(clc_sel, size = 600, method = "random")
points_train$class = as.factor(points_train$CODE_18)

# 4. EKSTRAKCJA WARTOŚCI SPEKTRALNYCH DLA PUNKTÓW
train_data = extract(r_2023, points_train, df = TRUE)
train_data$class = points_train$class
train_data = na.omit(train_data)  # usunięcie punktów z brakującymi wartościami

# 5. PODZIAŁ NA ZBIÓR TRENINGOWY I TESTOWY
set.seed(123)
split_index = createDataPartition(train_data$class, p = 0.7, list = FALSE)
train_set = train_data[split_index, ]
test_set = train_data[-split_index, ]

# 6. TRENING MODELU RANDOM FOREST
model_rf = randomForest(as.factor(class) ~ blue + green + red + nir,
                        data = train_set,
                        ntree = 200,
                        importance = TRUE)

# 7. KLASYFIKACJA CAŁEGO OBRAZU
classified_2023 <- predict(r_2023, model_rf, type = "response")

# 8. WIZUALIZACJA WYNIKU KLASYFIKACJI
par(mfrow = c(1, 3))

colors <- brewer.pal(length(selected_classes), "Set1")
plot(classified_2023, col = colors, main = "Wynik klasyfikacji pokrycia terenu (RF) - 2023")
legend("topright", legend = selected_classes, fill = colors, title = "Klasy CLC")

# KLASYFIKACJA DLA 2018 ROKU
print("Klasyfikacja dla roku 2018...")
classified_2018 <- predict(r_2018, model_rf, type = "response")

# WIZUALIZACJA WYNIKU KLASYFIKACJI DLA 2018
plot(classified_2018, col = colors, main = "Wynik klasyfikacji pokrycia terenu (RF) - 2018")
legend("topright", legend = selected_classes, fill = colors, title = "Klasy CLC")

# 9. OCENA SKUTECZNOŚCI NA ZBIORZE TESTOWYM
pred_test <- predict(model_rf, test_set)
conf_matrix <- confusionMatrix(pred_test, test_set$class)
print(conf_matrix)

# 10. MACIERZ POMYŁEK I NAJCZĘSTSZE BŁĘDY
print(conf_matrix$table)
print(conf_matrix$byClass[, c("Precision", "Recall", "F1")])

# 11. PORÓWNANIE Z MAPĄ REFERENCYJNĄ (jeśli masz CLC lub S2GLC w rastrowej formie)
library(terra)


ptlz = vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTLZ_A.shp")   # lasy
ptwp = vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTWP_A.shp")   # wody
pttr = vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTTR_A.shp")   # pastwiska,trawiasta
ptzb = vect("bdot10k/bdot10k_0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTZB_A.shp")   # zabudowa

#  Przekształcenie do układu współrzędnych klasyfikacji
crs_utm = crs(r_2023)
ptlz = project(ptlz, crs_utm)
ptwp = project(ptwp, crs_utm)
pttr = project(pttr, crs_utm)
ptzb = project(ptzb, crs_utm)
#  Dodajanie kolumn KOD_KLASY (uproszczone klasy tematyczne) 

ptlz$KOD_KLASY = 1  # las
ptwp$KOD_KLASY = 2  # wody
pttr$KOD_KLASY = 3  # trawiasta, grunty rolne
ptzb$KOD_KLASY = 4 #zabudowy
#  Połączenie warstwy w jedną ===
bdot_all = rbind(ptlz, ptwp, pttr, ptzb)

#  Rasteryzacja na podstawie pola KOD_KLASY
bdot_raster = rasterize(bdot_all, r_2023[[1]], field = "KOD_KLASY")

#  Przycięcie i zamaska do granic powiatu 
bdot_crop = crop(bdot_raster, powiat_utm)
bdot_crop = mask(bdot_crop, powiat_utm)

#  Porównanie z klasyfikacją (np. z randomForest lub inną metodą) 
# Zakładamy, że wynik klasyfikacji to raster "classified"
compare <- cbind(values(classified_2023), values(bdot_crop))
conf_matrix <- table(klasyfikacja = compare[,1], bdot = compare[,2])
print(conf_matrix)

# Kolory dla klas
kolory = c("darkgreen", "blue", "gold", "red")

bdot_crop = as.factor(bdot_crop)
levels(bdot_crop) = data.frame(value = 1:4,
                               label = c("Lasy", "Wody", "Pola i pastwiska", "Zabudowa"))

# Użyj kolorów i automatycznej legendy
plot(bdot_crop, col = kolory, main = "Pokrycie terenu BDOT")

# 12. UDZIAŁ PROCENTOWY POSZCZEGÓLNYCH KLAS
freq_table_2023 <- freq(classified_2023)
freq_df_2023 <- data.frame(class = freq_table_2023[, "value"],
                           count = freq_table_2023[, "count"])
freq_df_2023$percent <- round(100 * freq_df_2023$count / sum(freq_df_2023$count), 2)
print("Udział procentowy klas w 2023:")
print(freq_df_2023)

# UDZIAŁ PROCENTOWY POSZCZEGÓLNYCH KLAS DLA 2018
print("Udział procentowy klas w 2018:")
freq_table_2018 <- freq(classified_2018)
freq_df_2018 <- data.frame(class = freq_table_2018[, "value"],
                           count = freq_table_2018[, "count"])
freq_df_2018$percent <- round(100 * freq_df_2018$count / sum(freq_df_2018$count), 2)
print(freq_df_2018)
