# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT WŁODAWSKI

# Wymagane pakiety
library("terra")
library("rstac")
library("sf")
library("caret")         # do podziału na zbiór treningowy i testowy
library("randomForest")  # klasyfikator Random Forest
library("RColorBrewer")  # ładne kolory do wizualizacji
options(timeout = 1200)

#======================================================================
#### KONFIGURACJA ŚCIEŻEK ####
#======================================================================
main_dir <- "D:/Nowy folder/studia/zap_projekt/zap_klasyfikacja_nadzorowana"
setwd(main_dir)

# Ścieżki do danych wejściowych
powiat_path <- file.path(main_dir, "powiat_wlodawski.shp")
clc_path <- file.path(main_dir, "CLC18_PL", "clc18_PL.shp")
bdot_dir <- file.path(main_dir, "bdot10k", "bdot10k_0619_SHP")
ptlz_path <- file.path(bdot_dir, "PL.PZGiK.3700.BDOT10k.0619__OT_PTLZ_A.shp")
ptwp_path <- file.path(bdot_dir, "PL.PZGiK.3700.BDOT10k.0619__OT_PTWP_A.shp")
pttr_path <- file.path(bdot_dir, "PL.PZGiK.3700.BDOT10k.0619__OT_PTTR_A.shp")
ptzb_path <- file.path(bdot_dir, "PL.PZGiK.3700.BDOT10k.0619__OT_PTZB_A.shp")

# Ścieżki do danych wyjściowych i tymczasowych
wyniki_dir <- file.path(main_dir, "wyniki")
sentinel_2023_dir <- file.path(main_dir, "sentinel_2023")
mozaika_2023_path <- file.path(wyniki_dir, "mozaika_2023.tif")

#======================================================================

# Ustawienia folderu
# setwd("D:/Nowy folder/studia/zap_projekt/zap_klasyfikacja_nadzorowana")

# 1. WCZYTANIE I PRZYGOTOWANIE GRANIC POWIATU
print("1. Wczytywanie granic powiatu...")
powiat <- vect(powiat_path)
powiat_wgs84 <- project(powiat, "EPSG:4326")
powiat_utm <- project(powiat, "EPSG:32634")

# 2. PRZYGOTOWANIE BBOX DLA WYSZUKIWANIA OBRAZÓW
bbox_wgs84 <- ext(powiat_wgs84)
margin <- 0.1
bbox_buf <- c(bbox_wgs84[1] - margin, bbox_wgs84[3] - margin,
              bbox_wgs84[2] + margin, bbox_wgs84[4] + margin)

kanaly <- c("blue", "green", "red", "nir")

#======================================================================
#### 1. POBIERANIE I PRZETWARZANIE DANYCH SENTINEL-2 DLA 2023 ROKU ####
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
tile_ids_2023 <- c("S2A_OPER_MSI_L2A_TL_2APS_20230829T133355_A042747_T34UFB_N05.09",
                   "S2A_OPER_MSI_L2A_TL_2APS_20230829T133355_A042747_T34UFC_N05.09")

base_dir_2023 <- sentinel_2023_dir
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
  print("Kafle 2023 są już pobrane.")
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
dir.create(wyniki_dir, showWarnings = FALSE)

if (file.exists(mozaika_2023_path)) {
  print("Wczytywanie istniejącej mozaiki 2023...")
  r_2023 <- rast(mozaika_2023_path)
} else {
  print("Tworzenie mozaiki 2023...")
  rastry_2023 <- list.files(sentinel_2023_dir, pattern = "\\.tif$", full.names = TRUE, recursive = TRUE)
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

# WIZUALIZACJA 1: Obraz RGB

dev.new() # Otwarcie nowego, czystego okna graficznego
plotRGB(r_2023, r = 3, g = 2, b = 1, stretch = "lin", main = "Sentinel-2 RGB 2023", mar=c(2,2,4,1))
plot(powiat_utm, add = TRUE, col = NA, border = "red", lwd = 2)
#RP
#======================================================================
#### 1. TRENING NA PODSTAWIE CLC I PORÓWNANIE Z BDOT10K ####
#=======================================================================


# 1. WCZYTANIE I PRZYGOTOWANIE DANYCH REFERENCYJNYCH (np. CLC)
# Zakładamy, że warstwa CLC jest już pobrana i dotyczy całej Polski lub regionu
clc = vect(clc_path)  # zamień ścieżkę na odpowiednią
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
classified = predict(r_2023, model_rf, type = "response")

# 8. WIZUALIZACJA WYNIKU KLASYFIKACJI
par(mfrow = c(1, 2))
colors = brewer.pal(length(selected_classes), "Set1")
plot(classified, col = colors, main = "Wynik klasyfikacji pokrycia terenu (RF)")
legend("topright", legend = selected_classes, fill = colors, title = "Klasy CLC")

# 9. OCENA SKUTECZNOŚCI NA ZBIORZE TESTOWYM
pred_test = predict(model_rf, test_set)
conf_matrix = confusionMatrix(pred_test, test_set$class)
print(conf_matrix)

# 10. MACIERZ POMYŁEK I NAJCZĘSTSZE BŁĘDY
print(conf_matrix$table)
print(conf_matrix$byClass[, c("Precision", "Recall", "F1")])

# 11. PORÓWNANIE Z MAPĄ REFERENCYJNĄ (jeśli masz CLC lub S2GLC w rastrowej formie)
library(terra)


ptlz = vect(ptlz_path)   # lasy
ptwp = vect(ptwp_path)   # wody
pttr = vect(pttr_path)   # pastwiska,trawiasta
ptzb = vect(ptzb_path)   # zabudowa

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

# Porównanie z klasyfikacją
compare = cbind(values(classified), values(bdot_crop))
conf_matrix_bdot = table(klasyfikacja = compare[,1], bdot = compare[,2])
print("Macierz pomyłek Klasyfikacja vs BDOT:")
print(conf_matrix_bdot)

# 11. WIZUALIZACJA PORÓWNAWCZA
# Definicja czytelnych etykiet dla legend
class_labels_rf <- c("Zabudowa", "Pola uprawne", "Pastwiska", "Lasy", "Wody")
class_labels_bdot <- c("Lasy", "Wody", "Pola i pastwiska", "Zabudowa")

# Ustawienia kolorów
colors_rf <- brewer.pal(length(selected_classes), "Set1")
colors_bdot <- c("darkgreen", "blue", "gold", "red")

# Ustawienie okna na dwa wykresy z odpowiednimi marginesami
# mar = c(dół, lewo, góra, prawo); oma = zewnętrzny margines
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1), oma = c(0, 0, 2, 0))

# Wykres 1: Wynik klasyfikacji Random Forest
plot(classified, col = colors_rf, main = "Klasyfikacja RF (2023)", legend = FALSE)
legend("topright", legend = class_labels_rf, fill = colors_rf, title = "Klasy", bg = "white", cex = 0.8, inset = 0.02)

# Wykres 2: Mapa referencyjna BDOT
bdot_crop_factor = as.factor(bdot_crop)
levels(bdot_crop_factor) = data.frame(value = 1:4, label = class_labels_bdot)
plot(bdot_crop_factor, col = colors_bdot, main = "Pokrycie terenu BDOT", legend = FALSE)
legend("topright", legend = class_labels_bdot, fill = colors_bdot, title = "Klasy BDOT", bg = "white", cex = 0.8, inset = 0.02)

# Resetowanie ustawień graficznych jest ważne, aby uniknąć wpływu na kolejne wykresy.
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


# 12. UDZIAŁ PROCENTOWY POSZCZEGÓLNYCH KLAS
freq_table = freq(classified)
freq_df = data.frame(class = freq_table[, "value"],
                     count = freq_table[, "count"])
freq_df$percent = round(100 * freq_df$count / sum(freq_df$count), 2)
print(freq_df)



# 12. WIZUALIZACJA Z UJEDNOLICONYMI KOLORAMI (KLASYFIKACJA RF vs BDOT)

# Stworzenie macierzy do reklasyfikacji w celu ujednolicenia klas
# Mapujemy wartości rastra klasyfikacji (1-5) na kody klas BDOT (1-4)
# Wartości RF (z sorted(selected_classes)): 1="112", 2="211", 3="231", 4="311", 5="511"
# Mapowanie na BDOT: Zabudowa(4), Pola(3), Pastwiska(3), Lasy(1), Wody(2)
reclass_matrix <- matrix(c(1, 4,  # Zabudowa
                           2, 3,  # Pola uprawne
                           3, 3,  # Pastwiska
                           4, 1,  # Lasy
                           5, 2), # Wody
                         ncol = 2, byrow = TRUE)

# Reklasyfikacja rastra RF
classified_reclass <- classify(classified, reclass_matrix, others = NA)

# Ustawienie okna na dwa wykresy
par(mfrow = c(1, 2), mar = c(2, 2, 3, 1), oma = c(0, 0, 2, 0))

# Wykres 1: Reklasyfikowany wynik RF z kolorami BDOT
plot(classified_reclass, col = colors_bdot, main = "Klasyfikacja RF (kolory BDOT)", legend = FALSE)
legend("topright", legend = class_labels_bdot, fill = colors_bdot, title = "Klasy", bg = "white", cex = 0.8, inset = 0.02)

# Wykres 2: Oryginalna mapa BDOT
plot(bdot_crop_factor, col = colors_bdot, main = "Pokrycie terenu BDOT", legend = FALSE)
legend("topright", legend = class_labels_bdot, fill = colors_bdot, title = "Klasy BDOT", bg = "white", cex = 0.8, inset = 0.02)

# Resetowanie ustawień na końcu jest dobrą praktyką
par(mfrow = c(1, 1), mar = c(5.1, 4.1, 4.1, 2.1), oma = c(0, 0, 0, 0))


# 13. UDZIAŁ PROCENTOWY POSZCZEGÓLNYCH KLAS
freq_table = freq(classified)
freq_df = data.frame(class = freq_table[, "value"],
                     count = freq_table[, "count"])
freq_df$percent = round(100 * freq_df$count / sum(freq_df$count), 2)
print(freq_df)