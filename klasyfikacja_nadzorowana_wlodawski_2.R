# KLASYFIKACJA NADZOROWANA POKRYCIA TERENU - POWIAT WŁODAWSKI

# Wymagane pakiety
library("terra")
library("rstac")
library("sf")
library("caret")         # do podziału na zbiór treningowy i testowy
library("randomForest")  # klasyfikator Random Forest
library("RColorBrewer")  # ładne kolory do wizualizacji
options(timeout = 1200)

# wczytanie i przygotowanie granic powiatu
setwd("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie")
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



#RP
############TRENING NA PODSTAWIE CLC I PORÓWNANIE Z BDOT10K#########################



# 1. WCZYTANIE I PRZYGOTOWANIE DANYCH REFERENCYJNYCH (np. CLC)
# Zakładamy, że warstwa CLC jest już pobrana i dotyczy całej Polski lub regionu
clc = vect("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie/clc18_PL.shp")  # zamień ścieżkę na odpowiednią
clc_utm = crop(project(clc, crs(r_final)), powiat_utm)  # przekształcenie i docięcie do powiatu

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
train_data = extract(r_final, points_train, df = TRUE)
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
classified = predict(r_final, model_rf, type = "response")

# 8. WIZUALIZACJA WYNIKU KLASYFIKACJI
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


ptlz = vect("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie/0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTLZ_A.shp")   # lasy
ptwp = vect("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie/0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTWP_A.shp")   # wody
pttr = vect("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie/0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTTR_A.shp")   # pastwiska,trawiasta
ptzb = vect("E:/geoinfIIstopien/isemestr/zaawansowane_analizy_przestrzenne/zaliczenie/0619_SHP/PL.PZGiK.3700.BDOT10k.0619__OT_PTZB_A.shp")   # zabudowa

#  Przekształcenie do układu współrzędnych klasyfikacji
crs_utm = crs(r_final)
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
bdot_raster = rasterize(bdot_all, r_final[[1]], field = "KOD_KLASY")

#  Przycięcie i zamaska do granic powiatu 
bdot_crop = crop(bdot_raster, powiat_utm)
bdot_crop = mask(bdot_crop, powiat_utm)

#  Porównanie z klasyfikacją (np. z randomForest lub inną metodą) 
# Zakładamy, że wynik klasyfikacji to raster "classified"
compare = cbind(values(classified), values(bdot_crop))
conf_matrix = table(klasyfikacja = compare[,1], bdot = compare[,2])
print(conf_matrix)

# Kolory dla klas
kolory = c("darkgreen", "blue", "gold", "red")

bdot_crop = as.factor(bdot_crop)
levels(bdot_crop) = data.frame(value = 1:4,
                                label = c("Lasy", "Wody", "Pola i pastwiska", "Zabudowa"))

# Użyj kolorów i automatycznej legendy
plot(bdot_crop, col = kolory, main = "Pokrycie terenu BDOT")



# 12. UDZIAŁ PROCENTOWY POSZCZEGÓLNYCH KLAS
freq_table = freq(classified)
freq_df = data.frame(class = freq_table[, "value"],
                      count = freq_table[, "count"])
freq_df$percent = round(100 * freq_df$count / sum(freq_df$count), 2)
print(freq_df)




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