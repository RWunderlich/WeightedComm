source(file = "main.R")
require(RStoolbox)
require(virtualspecies)

Forest_list <- list.files(path = "Layers", pattern = "Earth", full.names = T)[1:4]
Forest_stack <- stack(x = Forest_list)
plot(Forest_stack)

EE_1_Forest <- calc(x = Forest_stack, fun = sum)
EE_1_Forest@data@names <- "EE_1_Forest"
plot(EE_1_Forest)
rm(Forest_list, Forest_stack)

EE_2_Shrub <- raster(x = "Layers/EarthEnv_05_Shrubs.tif")
EE_2_Shrub@data@names <- "EE_2_Shrub"

EE_3_Herbac <- raster(x = "Layers/EarthEnv_06_Herbaceous.tif")
EE_3_Herbac@data@names <- "EE_3_Herbac"

EE_4_Cultiv <- raster(x = "Layers/EarthEnv_07_Cultivated.tif")
EE_4_Cultiv@data@names <- "EE_4_Cultiv"

EE_5_Flood <- raster(x = "Layers/EarthEnv_08_Flooded.tif")
EE_5_Flood@data@names <- "EE_5_Flood"

EE_6_Water <- raster(x = "Layers/EarthEnv_12_Water.tif")
EE_6_Water@data@names <- "EE_6_Water"

EE_7_Built <- raster(x = "Layers/EarthEnv_09_Built.tif")
EE_7_Built@data@names <- "EE_7_Built"

HeatLoad <- raster(x = list.files(path = "Layers", pattern = "tif", full.names = T)[17])
HeatLoad@data@names <- "HeatLoad"

layer_list <- list.files(path = "Layers", pattern = "tif", full.names = T)[1:6]

my_stack <- stack(layer_list, EE_1_Forest, EE_2_Shrub, EE_3_Herb, EE_4_Cultiv, EE_5_Flood, EE_6_Water, EE_7_Built, HeatLoad)

my_PCA <- rasterPCA(img = my_stack, nSamples = NULL, nComp = 10, spca = 1, maskCheck = T)

my_PCs <- stack(my_PCA$map)[[1:8]]

my_PCtable <- data.frame("min" = rep(NA, 8), "med" = rep(NA, 8), "max" = rep(NA, 8))

for (j in (1:3)) {
  for (i in (1:8)) {
    if (j == 1) {
      my_PCtable[i, j] <- minValue(my_PCs[[i]])
    }
    if (j == 2) {
      my_PCtable[i, j] <- median(values(my_PCs[[i]]), na.rm = T)
    }
    if (j == 3) {
      my_PCtable[i, j] <- maxValue(my_PCs[[i]])
    }
  }
}

my_VS_mean_SD <- array(data = rep(NA, 80), dim = c(5, 8, 2), dimnames = list(c("VS1", "VS2", "VS3", "VS4", "VS5"),
                                                                             c("PC1", "PC2", "PC3","PC4", "PC5","PC6", "PC7", "PC8"),
                                                                             c("mean", "sd")))
set.seed(2607)
for (i in (1:8)) {
  my_VS_mean_SD[1,i,1] <- my_PCtable[i,2]
  my_VS_mean_SD[2,i,1] <- my_PCtable[i,2]
  my_VS_mean_SD[1,i,2] <- 0.5
  my_VS_mean_SD[2,i,2] <- 2.5
  for (j in (3:5)) {
    my_VS_mean_SD[j,i,1] <- sample(x = (my_PCtable[i,1]:my_PCtable[i,3]), size = 1)
    my_VS_mean_SD[j,i,2] <- runif(n = 1, min = 0.5, max = 2.5)
  }
}

my_VS_parameters_1 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,1,1], sd = my_VS_mean_SD[1,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,2,1], sd = my_VS_mean_SD[1,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,3,1], sd = my_VS_mean_SD[1,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,4,1], sd = my_VS_mean_SD[1,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,5,1], sd = my_VS_mean_SD[1,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,6,1], sd = my_VS_mean_SD[1,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,7,1], sd = my_VS_mean_SD[1,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,8,1], sd = my_VS_mean_SD[1,8,2]))

my_VS_1 <- generateSpFromFun(raster.stack = my_PCs, scale = T
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_1,
                             plot = TRUE)

my_VS_parameters_2 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,1,1], sd = my_VS_mean_SD[2,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,2,1], sd = my_VS_mean_SD[2,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,3,1], sd = my_VS_mean_SD[2,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,4,1], sd = my_VS_mean_SD[2,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,5,1], sd = my_VS_mean_SD[2,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,6,1], sd = my_VS_mean_SD[2,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,7,1], sd = my_VS_mean_SD[2,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,8,1], sd = my_VS_mean_SD[2,8,2]))

my_VS_2 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_2,
                             plot = TRUE)

my_VS_parameters_3 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,1,1], sd = my_VS_mean_SD[3,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,2,1], sd = my_VS_mean_SD[3,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,3,1], sd = my_VS_mean_SD[3,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,4,1], sd = my_VS_mean_SD[3,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,5,1], sd = my_VS_mean_SD[3,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,6,1], sd = my_VS_mean_SD[3,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,7,1], sd = my_VS_mean_SD[3,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,8,1], sd = my_VS_mean_SD[3,8,2]))

my_VS_3 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_3,
                             plot = TRUE)

my_VS_parameters_4 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,1,1], sd = my_VS_mean_SD[4,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,2,1], sd = my_VS_mean_SD[4,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,3,1], sd = my_VS_mean_SD[4,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,4,1], sd = my_VS_mean_SD[4,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,5,1], sd = my_VS_mean_SD[4,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,6,1], sd = my_VS_mean_SD[4,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,7,1], sd = my_VS_mean_SD[4,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,8,1], sd = my_VS_mean_SD[4,8,2]))

my_VS_4 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_4,
                             plot = TRUE)

my_VS_parameters_5 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,1,1], sd = my_VS_mean_SD[5,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,2,1], sd = my_VS_mean_SD[5,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,3,1], sd = my_VS_mean_SD[5,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,4,1], sd = my_VS_mean_SD[5,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,5,1], sd = my_VS_mean_SD[5,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,6,1], sd = my_VS_mean_SD[5,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,7,1], sd = my_VS_mean_SD[5,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,8,1], sd = my_VS_mean_SD[5,8,2]))

my_VS_5 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_5,
                             plot = TRUE)

set.seed(2607)
my_PA_1 <- convertToPA(x = my_VS_1, PA.method = "probability", species.prevalence = 0.10, alpha = -0.01, plot = T)
my_occurrence_data_PD_1_po <- sampleOccurrences(x = my_PA_1, n = 1000, type = "presence only", detection.probability = 1, plot = T)
my_occurrence_data_ID_1_po <- sampleOccurrences(x = my_PA_1, n = 1000, type = "presence only", detection.probability = 0.75, plot = T)
my_occurrence_data_PD_1 <- sampleOccurrences(x = my_PA_1, n = 5000, type = "presence-absence", detection.probability = 1, plot = T)
my_occurrence_data_ID_1 <- sampleOccurrences(x = my_PA_1, n = 5000, type = "presence-absence", detection.probability = 0.75, plot = T)


set.seed(2607)
my_PA_2 <- convertToPA(x = my_VS_2, PA.method = "probability", species.prevalence = 0.25, alpha = -0.01, plot = T)
my_occurrence_data_PD_2_po <- sampleOccurrences(x = my_PA_2, n = 1000, type = "presence only", detection.probability = 1, plot = T)
my_occurrence_data_ID_2_po <- sampleOccurrences(x = my_PA_2, n = 1000, type = "presence only", detection.probability = 0.75, plot = T)
my_occurrence_data_PD_2 <- sampleOccurrences(x = my_PA_2, n = 5000, type = "presence-absence", detection.probability = 1, plot = T)
my_occurrence_data_ID_2 <- sampleOccurrences(x = my_PA_2, n = 5000, type = "presence-absence", detection.probability = 0.75, plot = T)

set.seed(2607)
my_PA_3 <- convertToPA(x = my_VS_3, PA.method = "probability", species.prevalence = 0.175, alpha = -0.01, plot = T)
my_occurrence_data_PD_3_po <- sampleOccurrences(x = my_PA_3, n = 1000, type = "presence only", detection.probability = 1, plot = T)
my_occurrence_data_ID_3_po <- sampleOccurrences(x = my_PA_3, n = 1000, type = "presence only", detection.probability = 0.75, plot = T)
my_occurrence_data_PD_3 <- sampleOccurrences(x = my_PA_3, n = 5000, type = "presence-absence", detection.probability = 1, plot = T)
my_occurrence_data_ID_3 <- sampleOccurrences(x = my_PA_3, n = 5000, type = "presence-absence", detection.probability = 0.75, plot = T)

set.seed(2607)
my_PA_4 <- convertToPA(x = my_VS_4, PA.method = "probability", species.prevalence = 0.175, alpha = -0.01, plot = T)
my_occurrence_data_PD_4_po <- sampleOccurrences(x = my_PA_4, n = 1000, type = "presence only", detection.probability = 1, plot = T)
my_occurrence_data_ID_4_po <- sampleOccurrences(x = my_PA_4, n = 1000, type = "presence only", detection.probability = 0.75, plot = T)
my_occurrence_data_PD_4 <- sampleOccurrences(x = my_PA_4, n = 5000, type = "presence-absence", detection.probability = 1, plot = T)
my_occurrence_data_ID_4 <- sampleOccurrences(x = my_PA_4, n = 5000, type = "presence-absence", detection.probability = 0.75, plot = T)

set.seed(2607)
my_PA_5 <- convertToPA(x = my_VS_5, PA.method = "probability", species.prevalence = 0.175, alpha = -0.01, plot = T)
my_occurrence_data_PD_5_po <- sampleOccurrences(x = my_PA_5, n = 1000, type = "presence only", detection.probability = 1, plot = T)
my_occurrence_data_ID_5_po <- sampleOccurrences(x = my_PA_5, n = 1000, type = "presence only", detection.probability = 0.75, plot = T)
my_occurrence_data_PD_5 <- sampleOccurrences(x = my_PA_5, n = 5000, type = "presence-absence", detection.probability = 1, plot = T)
my_occurrence_data_ID_5 <- sampleOccurrences(x = my_PA_5, n = 5000, type = "presence-absence", detection.probability = 0.75, plot = T)
