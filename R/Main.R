source(file = "main.R")
require(RStoolbox)
require(virtualspecies)
require(dplyr)

(Forest_list <- list.files(path = "Layers", pattern = "Earth", full.names = T)[1:4])
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

(layer_list <- list.files(path = "Layers", pattern = "tif", full.names = T)[c(1:6,17:18)])

my_stack <- stack(layer_list, EE_1_Forest, EE_2_Shrub, EE_3_Herbac, EE_4_Cultiv, EE_5_Flood, EE_6_Water, EE_7_Built)

my_stack[my_stack <= -99] <- NA # make sure no mistakes in NA values

my_PCA <- rasterPCA(img = my_stack, nSamples = NULL, nComp = 10, spca = 1, maskCheck = T)

my_PCs <- stack(my_PCA$map)[[1:9]]

for (i in (1:nlayers(my_PCs))) {
  writeRaster(x = my_PCs[[i]], filename = paste0("Layers/PCA/PC", as.character(i), ".tif"), format = "GTiff", overwrite = T)
}

my_PCtable <- data.frame("min" = rep(NA, 9), "med" = rep(NA, 9), "max" = rep(NA, 9))

for (j in (1:3)) {
  for (i in (1:9)) {
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

my_VS_mean_SD <- array(data = rep(NA, 80), dim = c(13, 9, 2), dimnames = list(c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07",
                                                                                "VS08", "VS09", "VS10", "VS11", "VS12", "VS13"),
                                                                              c("PC1", "PC2", "PC3","PC4", "PC5","PC6", "PC7", "PC8", "PC9"),
                                                                              c("mean", "sd")))
set.seed(2607)
for (i in (1:9)) {
  my_VS_mean_SD[1,i,1] <- my_PCtable[i,2]
  my_VS_mean_SD[2,i,1] <- my_PCtable[i,2]
  my_VS_mean_SD[1,i,2] <- 0.5
  my_VS_mean_SD[2,i,2] <- 2.5
  for (j in (3:13)) {
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
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,8,1], sd = my_VS_mean_SD[1,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[1,9,1], sd = my_VS_mean_SD[1,9,2]))

my_VS_1 <- generateSpFromFun(raster.stack = my_PCs,
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
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,8,1], sd = my_VS_mean_SD[2,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[2,9,1], sd = my_VS_mean_SD[2,9,2]))

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
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,8,1], sd = my_VS_mean_SD[3,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[3,9,1], sd = my_VS_mean_SD[3,9,2]))

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
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,8,1], sd = my_VS_mean_SD[4,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[4,9,1], sd = my_VS_mean_SD[4,9,2]))

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
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,8,1], sd = my_VS_mean_SD[5,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[5,9,1], sd = my_VS_mean_SD[5,9,2]))

my_VS_5 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_5,
                             plot = TRUE)

my_VS_parameters_6 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,1,1], sd = my_VS_mean_SD[6,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,2,1], sd = my_VS_mean_SD[6,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,3,1], sd = my_VS_mean_SD[6,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,4,1], sd = my_VS_mean_SD[6,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,5,1], sd = my_VS_mean_SD[6,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,6,1], sd = my_VS_mean_SD[6,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,7,1], sd = my_VS_mean_SD[6,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,8,1], sd = my_VS_mean_SD[6,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[6,9,1], sd = my_VS_mean_SD[6,9,2]))

my_VS_6 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_6,
                             plot = TRUE)

my_VS_parameters_7 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,1,1], sd = my_VS_mean_SD[7,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,2,1], sd = my_VS_mean_SD[7,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,3,1], sd = my_VS_mean_SD[7,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,4,1], sd = my_VS_mean_SD[7,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,5,1], sd = my_VS_mean_SD[7,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,6,1], sd = my_VS_mean_SD[7,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,7,1], sd = my_VS_mean_SD[7,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,8,1], sd = my_VS_mean_SD[7,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[7,9,1], sd = my_VS_mean_SD[7,9,2]))

my_VS_7 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_7,
                             plot = TRUE)

my_VS_parameters_8 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,1,1], sd = my_VS_mean_SD[8,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,2,1], sd = my_VS_mean_SD[8,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,3,1], sd = my_VS_mean_SD[8,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,4,1], sd = my_VS_mean_SD[8,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,5,1], sd = my_VS_mean_SD[8,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,6,1], sd = my_VS_mean_SD[8,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,7,1], sd = my_VS_mean_SD[8,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,8,1], sd = my_VS_mean_SD[8,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[8,9,1], sd = my_VS_mean_SD[8,9,2]))

my_VS_8 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_8,
                             plot = TRUE)

my_VS_parameters_9 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,1,1], sd = my_VS_mean_SD[9,1,2]),
                                      PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,2,1], sd = my_VS_mean_SD[9,2,2]),
                                      PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,3,1], sd = my_VS_mean_SD[9,3,2]),
                                      PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,4,1], sd = my_VS_mean_SD[9,4,2]),
                                      PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,5,1], sd = my_VS_mean_SD[9,5,2]),
                                      PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,6,1], sd = my_VS_mean_SD[9,6,2]),
                                      PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,7,1], sd = my_VS_mean_SD[9,7,2]),
                                      PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,8,1], sd = my_VS_mean_SD[9,8,2]),
                                      PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[9,9,1], sd = my_VS_mean_SD[9,9,2]))

my_VS_9 <- generateSpFromFun(raster.stack = my_PCs,
                             species.type = "additive",
                             rescale = T,
                             parameters = my_VS_parameters_9,
                             plot = TRUE)

my_VS_parameters_10 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,1,1], sd = my_VS_mean_SD[10,1,2]),
                                       PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,2,1], sd = my_VS_mean_SD[10,2,2]),
                                       PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,3,1], sd = my_VS_mean_SD[10,3,2]),
                                       PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,4,1], sd = my_VS_mean_SD[10,4,2]),
                                       PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,5,1], sd = my_VS_mean_SD[10,5,2]),
                                       PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,6,1], sd = my_VS_mean_SD[10,6,2]),
                                       PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,7,1], sd = my_VS_mean_SD[10,7,2]),
                                       PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,8,1], sd = my_VS_mean_SD[10,8,2]),
                                       PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[10,9,1], sd = my_VS_mean_SD[10,9,2]))

my_VS_10 <- generateSpFromFun(raster.stack = my_PCs,
                              species.type = "additive",
                              rescale = T,
                              parameters = my_VS_parameters_10,
                              plot = TRUE)

my_VS_parameters_11 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,1,1], sd = my_VS_mean_SD[11,1,2]),
                                       PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,2,1], sd = my_VS_mean_SD[11,2,2]),
                                       PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,3,1], sd = my_VS_mean_SD[11,3,2]),
                                       PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,4,1], sd = my_VS_mean_SD[11,4,2]),
                                       PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,5,1], sd = my_VS_mean_SD[11,5,2]),
                                       PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,6,1], sd = my_VS_mean_SD[11,6,2]),
                                       PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,7,1], sd = my_VS_mean_SD[11,7,2]),
                                       PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,8,1], sd = my_VS_mean_SD[11,8,2]),
                                       PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[11,9,1], sd = my_VS_mean_SD[11,9,2]))

my_VS_11 <- generateSpFromFun(raster.stack = my_PCs,
                              species.type = "additive",
                              rescale = T,
                              parameters = my_VS_parameters_11,
                              plot = TRUE)

my_VS_parameters_12 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,1,1], sd = my_VS_mean_SD[12,1,2]),
                                       PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,2,1], sd = my_VS_mean_SD[12,2,2]),
                                       PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,3,1], sd = my_VS_mean_SD[12,3,2]),
                                       PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,4,1], sd = my_VS_mean_SD[12,4,2]),
                                       PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,5,1], sd = my_VS_mean_SD[12,5,2]),
                                       PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,6,1], sd = my_VS_mean_SD[12,6,2]),
                                       PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,7,1], sd = my_VS_mean_SD[12,7,2]),
                                       PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,8,1], sd = my_VS_mean_SD[12,8,2]),
                                       PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[12,9,1], sd = my_VS_mean_SD[12,9,2]))

my_VS_12 <- generateSpFromFun(raster.stack = my_PCs,
                              species.type = "additive",
                              rescale = T,
                              parameters = my_VS_parameters_12,
                              plot = TRUE)

my_VS_parameters_13 <- formatFunctions(PC1 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,1,1], sd = my_VS_mean_SD[13,1,2]),
                                       PC2 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,2,1], sd = my_VS_mean_SD[13,2,2]),
                                       PC3 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,3,1], sd = my_VS_mean_SD[13,3,2]),
                                       PC4 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,4,1], sd = my_VS_mean_SD[13,4,2]),
                                       PC5 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,5,1], sd = my_VS_mean_SD[13,5,2]),
                                       PC6 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,6,1], sd = my_VS_mean_SD[13,6,2]),
                                       PC7 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,7,1], sd = my_VS_mean_SD[13,7,2]),
                                       PC8 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,8,1], sd = my_VS_mean_SD[13,8,2]),
                                       PC9 = c(fun = 'dnorm', mean = my_VS_mean_SD[13,9,1], sd = my_VS_mean_SD[13,9,2]))

my_VS_13 <- generateSpFromFun(raster.stack = my_PCs,
                              species.type = "additive",
                              rescale = T,
                              parameters = my_VS_parameters_13,
                              plot = TRUE)

##########################
# Generate true distribution
# and sample occurrence data
##########################

set.seed(2607)
my_PA_1 <- convertToPA(x = my_VS_1, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_1$pa.raster, filename = "VS/Truth/my_PA_1.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_1_po <- sampleOccurrences(x = my_PA_1, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_1_po$sample.points[,1:2], file = "VS/Locations/PD_1_po.csv", row.names = F)
write.csv(x = PD_1_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/1/1_po.csv", row.names = F)     #1
write.csv(x = PD_1_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/1/2_po.csv", row.names = F)   #2
write.csv(x = PD_1_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/1/3_po.csv", row.names = F)   #3
write.csv(x = PD_1_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/1/4_po.csv", row.names = F)   #4
write.csv(x = PD_1_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/1/5_po.csv", row.names = F)   #5
write.csv(x = PD_1_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/1/6_po.csv", row.names = F)   #6
write.csv(x = PD_1_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/1/7_po.csv", row.names = F)   #7
write.csv(x = PD_1_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/1/8_po.csv", row.names = F)   #8
write.csv(x = PD_1_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/1/9_po.csv", row.names = F)   #9
write.csv(x = PD_1_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/1/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_2 <- convertToPA(x = my_VS_2, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_2$pa.raster, filename = "VS/Truth/my_PA_2.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_2_po <- sampleOccurrences(x = my_PA_2, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_2_po$sample.points[,1:2], file = "VS/Locations/PD_2_po.csv", row.names = F)
write.csv(x = PD_2_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/2/1_po.csv", row.names = F)     #1
write.csv(x = PD_2_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/2/2_po.csv", row.names = F)   #2
write.csv(x = PD_2_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/2/3_po.csv", row.names = F)   #3
write.csv(x = PD_2_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/2/4_po.csv", row.names = F)   #4
write.csv(x = PD_2_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/2/5_po.csv", row.names = F)   #5
write.csv(x = PD_2_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/2/6_po.csv", row.names = F)   #6
write.csv(x = PD_2_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/2/7_po.csv", row.names = F)   #7
write.csv(x = PD_2_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/2/8_po.csv", row.names = F)   #8
write.csv(x = PD_2_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/2/9_po.csv", row.names = F)   #9
write.csv(x = PD_2_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/2/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_3 <- convertToPA(x = my_VS_3, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_3$pa.raster, filename = "VS/Truth/my_PA_3.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_3_po <- sampleOccurrences(x = my_PA_3, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_3_po$sample.points[,1:2], file = "VS/Locations/PD_3_po.csv", row.names = F)
write.csv(x = PD_3_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/3/1_po.csv", row.names = F)     #1
write.csv(x = PD_3_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/3/2_po.csv", row.names = F)   #2
write.csv(x = PD_3_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/3/3_po.csv", row.names = F)   #3
write.csv(x = PD_3_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/3/4_po.csv", row.names = F)   #4
write.csv(x = PD_3_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/3/5_po.csv", row.names = F)   #5
write.csv(x = PD_3_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/3/6_po.csv", row.names = F)   #6
write.csv(x = PD_3_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/3/7_po.csv", row.names = F)   #7
write.csv(x = PD_3_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/3/8_po.csv", row.names = F)   #8
write.csv(x = PD_3_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/3/9_po.csv", row.names = F)   #9
write.csv(x = PD_3_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/3/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_4 <- convertToPA(x = my_VS_4, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_4$pa.raster, filename = "VS/Truth/my_PA_4.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_4_po <- sampleOccurrences(x = my_PA_4, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_4_po$sample.points[,1:2], file = "VS/Locations/PD_4_po.csv", row.names = F)
write.csv(x = PD_4_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/4/1_po.csv", row.names = F)     #1
write.csv(x = PD_4_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/4/2_po.csv", row.names = F)   #2
write.csv(x = PD_4_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/4/3_po.csv", row.names = F)   #3
write.csv(x = PD_4_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/4/4_po.csv", row.names = F)   #4
write.csv(x = PD_4_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/4/5_po.csv", row.names = F)   #5
write.csv(x = PD_4_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/4/6_po.csv", row.names = F)   #6
write.csv(x = PD_4_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/4/7_po.csv", row.names = F)   #7
write.csv(x = PD_4_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/4/8_po.csv", row.names = F)   #8
write.csv(x = PD_4_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/4/9_po.csv", row.names = F)   #9
write.csv(x = PD_4_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/4/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_5 <- convertToPA(x = my_VS_5, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_5$pa.raster, filename = "VS/Truth/my_PA_5.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_5_po <- sampleOccurrences(x = my_PA_5, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_5_po$sample.points[,1:2], file = "VS/Locations/PD_5_po.csv", row.names = F)
write.csv(x = PD_5_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/5/1_po.csv", row.names = F)     #1
write.csv(x = PD_5_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/5/2_po.csv", row.names = F)   #2
write.csv(x = PD_5_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/5/3_po.csv", row.names = F)   #3
write.csv(x = PD_5_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/5/4_po.csv", row.names = F)   #4
write.csv(x = PD_5_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/5/5_po.csv", row.names = F)   #5
write.csv(x = PD_5_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/5/6_po.csv", row.names = F)   #6
write.csv(x = PD_5_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/5/7_po.csv", row.names = F)   #7
write.csv(x = PD_5_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/5/8_po.csv", row.names = F)   #8
write.csv(x = PD_5_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/5/9_po.csv", row.names = F)   #9
write.csv(x = PD_5_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/5/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_6 <- convertToPA(x = my_VS_6, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_6$pa.raster, filename = "VS/Truth/my_PA_6.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_6_po <- sampleOccurrences(x = my_PA_6, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_6_po$sample.points[,1:2], file = "VS/Locations/PD_6_po.csv", row.names = F)
write.csv(x = PD_6_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/6/1_po.csv", row.names = F)     #1
write.csv(x = PD_6_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/6/2_po.csv", row.names = F)   #2
write.csv(x = PD_6_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/6/3_po.csv", row.names = F)   #3
write.csv(x = PD_6_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/6/4_po.csv", row.names = F)   #4
write.csv(x = PD_6_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/6/5_po.csv", row.names = F)   #5
write.csv(x = PD_6_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/6/6_po.csv", row.names = F)   #6
write.csv(x = PD_6_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/6/7_po.csv", row.names = F)   #7
write.csv(x = PD_6_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/6/8_po.csv", row.names = F)   #8
write.csv(x = PD_6_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/6/9_po.csv", row.names = F)   #9
write.csv(x = PD_6_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/6/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_7 <- convertToPA(x = my_VS_7, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_7$pa.raster, filename = "VS/Truth/my_PA_7.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_7_po <- sampleOccurrences(x = my_PA_7, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_7_po$sample.points[,1:2], file = "VS/Locations/PD_7_po.csv", row.names = F)
write.csv(x = PD_7_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/7/1_po.csv", row.names = F)     #1
write.csv(x = PD_7_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/7/2_po.csv", row.names = F)   #2
write.csv(x = PD_7_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/7/3_po.csv", row.names = F)   #3
write.csv(x = PD_7_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/7/4_po.csv", row.names = F)   #4
write.csv(x = PD_7_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/7/5_po.csv", row.names = F)   #5
write.csv(x = PD_7_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/7/6_po.csv", row.names = F)   #6
write.csv(x = PD_7_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/7/7_po.csv", row.names = F)   #7
write.csv(x = PD_7_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/7/8_po.csv", row.names = F)   #8
write.csv(x = PD_7_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/7/9_po.csv", row.names = F)   #9
write.csv(x = PD_7_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/7/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_8 <- convertToPA(x = my_VS_8, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_8$pa.raster, filename = "VS/Truth/my_PA_8.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_8_po <- sampleOccurrences(x = my_PA_8, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_8_po$sample.points[,1:2], file = "VS/Locations/PD_8_po.csv", row.names = F)
write.csv(x = PD_8_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/8/1_po.csv", row.names = F)     #1
write.csv(x = PD_8_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/8/2_po.csv", row.names = F)   #2
write.csv(x = PD_8_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/8/3_po.csv", row.names = F)   #3
write.csv(x = PD_8_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/8/4_po.csv", row.names = F)   #4
write.csv(x = PD_8_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/8/5_po.csv", row.names = F)   #5
write.csv(x = PD_8_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/8/6_po.csv", row.names = F)   #6
write.csv(x = PD_8_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/8/7_po.csv", row.names = F)   #7
write.csv(x = PD_8_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/8/8_po.csv", row.names = F)   #8
write.csv(x = PD_8_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/8/9_po.csv", row.names = F)   #9
write.csv(x = PD_8_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/8/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_9 <- convertToPA(x = my_VS_9, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_9$pa.raster, filename = "VS/Truth/my_PA_9.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_9_po <- sampleOccurrences(x = my_PA_9, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_9_po$sample.points[,1:2], file = "VS/Locations/PD_9_po.csv", row.names = F)
write.csv(x = PD_9_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/9/1_po.csv", row.names = F)     #1
write.csv(x = PD_9_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/9/2_po.csv", row.names = F)   #2
write.csv(x = PD_9_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/9/3_po.csv", row.names = F)   #3
write.csv(x = PD_9_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/9/4_po.csv", row.names = F)   #4
write.csv(x = PD_9_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/9/5_po.csv", row.names = F)   #5
write.csv(x = PD_9_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/9/6_po.csv", row.names = F)   #6
write.csv(x = PD_9_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/9/7_po.csv", row.names = F)   #7
write.csv(x = PD_9_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/9/8_po.csv", row.names = F)   #8
write.csv(x = PD_9_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/9/9_po.csv", row.names = F)   #9
write.csv(x = PD_9_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/9/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_10 <- convertToPA(x = my_VS_10, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_10$pa.raster, filename = "VS/Truth/my_PA_10.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_10_po <- sampleOccurrences(x = my_PA_10, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_10_po$sample.points[,1:2], file = "VS/Locations/PD_10_po.csv", row.names = F)
write.csv(x = PD_10_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/10/1_po.csv", row.names = F)     #1
write.csv(x = PD_10_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/10/2_po.csv", row.names = F)   #2
write.csv(x = PD_10_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/10/3_po.csv", row.names = F)   #3
write.csv(x = PD_10_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/10/4_po.csv", row.names = F)   #4
write.csv(x = PD_10_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/10/5_po.csv", row.names = F)   #5
write.csv(x = PD_10_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/10/6_po.csv", row.names = F)   #6
write.csv(x = PD_10_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/10/7_po.csv", row.names = F)   #7
write.csv(x = PD_10_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/10/8_po.csv", row.names = F)   #8
write.csv(x = PD_10_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/10/9_po.csv", row.names = F)   #9
write.csv(x = PD_10_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/10/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_11 <- convertToPA(x = my_VS_11, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_11$pa.raster, filename = "VS/Truth/my_PA_11.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_11_po <- sampleOccurrences(x = my_PA_11, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_11_po$sample.points[,1:2], file = "VS/Locations/PD_11_po.csv", row.names = F)
write.csv(x = PD_11_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/11/1_po.csv", row.names = F)     #1
write.csv(x = PD_11_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/11/2_po.csv", row.names = F)   #2
write.csv(x = PD_11_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/11/3_po.csv", row.names = F)   #3
write.csv(x = PD_11_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/11/4_po.csv", row.names = F)   #4
write.csv(x = PD_11_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/11/5_po.csv", row.names = F)   #5
write.csv(x = PD_11_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/11/6_po.csv", row.names = F)   #6
write.csv(x = PD_11_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/11/7_po.csv", row.names = F)   #7
write.csv(x = PD_11_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/11/8_po.csv", row.names = F)   #8
write.csv(x = PD_11_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/11/9_po.csv", row.names = F)   #9
write.csv(x = PD_11_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/11/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_12 <- convertToPA(x = my_VS_12, PA.method = "probability", species.prevalence = 0.10, alpha = -0.05, plot = T)
writeRaster(x = my_PA_12$pa.raster, filename = "VS/Truth/my_PA_12.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_12_po <- sampleOccurrences(x = my_PA_12, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_12_po$sample.points[,1:2], file = "VS/Locations/PD_12_po.csv", row.names = F)
write.csv(x = PD_12_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/12/1_po.csv", row.names = F)     #1
write.csv(x = PD_12_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/12/2_po.csv", row.names = F)   #2
write.csv(x = PD_12_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/12/3_po.csv", row.names = F)   #3
write.csv(x = PD_12_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/12/4_po.csv", row.names = F)   #4
write.csv(x = PD_12_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/12/5_po.csv", row.names = F)   #5
write.csv(x = PD_12_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/12/6_po.csv", row.names = F)   #6
write.csv(x = PD_12_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/12/7_po.csv", row.names = F)   #7
write.csv(x = PD_12_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/12/8_po.csv", row.names = F)   #8
write.csv(x = PD_12_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/12/9_po.csv", row.names = F)   #9
write.csv(x = PD_12_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/12/10_po.csv", row.names = F) #10

set.seed(2607)
my_PA_13 <- convertToPA(x = my_VS_13, PA.method = "probability", species.prevalence = 0.50, alpha = -0.05, plot = T)
writeRaster(x = my_PA_13$pa.raster, filename = "VS/Truth/my_PA_13.tif", format = "GTiff", overwrite = T)
# perfect detection
PD_13_po <- sampleOccurrences(x = my_PA_13, n = 1000, type = "presence only", detection.probability = 1, plot = T)
write.csv(x = PD_13_po$sample.points[,1:2], file = "VS/Locations/PD_13_po.csv", row.names = F)
write.csv(x = PD_13_po$sample.points[1:100,1:2], file = "VS/Locations/subsamples/13/1_po.csv", row.names = F)     #1
write.csv(x = PD_13_po$sample.points[101:200,1:2], file = "VS/Locations/subsamples/13/2_po.csv", row.names = F)   #2
write.csv(x = PD_13_po$sample.points[201:300,1:2], file = "VS/Locations/subsamples/13/3_po.csv", row.names = F)   #3
write.csv(x = PD_13_po$sample.points[301:400,1:2], file = "VS/Locations/subsamples/13/4_po.csv", row.names = F)   #4
write.csv(x = PD_13_po$sample.points[401:500,1:2], file = "VS/Locations/subsamples/13/5_po.csv", row.names = F)   #5
write.csv(x = PD_13_po$sample.points[501:600,1:2], file = "VS/Locations/subsamples/13/6_po.csv", row.names = F)   #6
write.csv(x = PD_13_po$sample.points[601:700,1:2], file = "VS/Locations/subsamples/13/7_po.csv", row.names = F)   #7
write.csv(x = PD_13_po$sample.points[701:800,1:2], file = "VS/Locations/subsamples/13/8_po.csv", row.names = F)   #8
write.csv(x = PD_13_po$sample.points[801:900,1:2], file = "VS/Locations/subsamples/13/9_po.csv", row.names = F)   #9
write.csv(x = PD_13_po$sample.points[901:1000,1:2], file = "VS/Locations/subsamples/13/10_po.csv", row.names = F) #10



##########################
# Calibrate models and
# predict across study area
##########################

#### DEBUG maxnet functions; first bugged: need to remove "..." from glmnet arguments (8th); also need to increase maxit def is 1e+06 -> try 1e+07
###  rest not loaded for unknown reasons
#' @import stats
#' @export
maxnet <-
  function(p, data, f=maxnet.formula(p, data), regmult=1.0,
           regfun=maxnet.default.regularization, ...)
  {
    if (anyNA(data)) stop("NA values in data table. Please remove them and rerun.")
    mm <- model.matrix(f, data)
    reg <- regfun(p,mm) * regmult
    weights <- p+(1-p)*100
    glmnet::glmnet.control(pmin=1.0e-8, fdev=0, devmax = 0.99)
    model <- glmnet::glmnet(x=mm, y=as.factor(p), family="binomial", standardize=F,
                            penalty.factor=reg, lambda=10^(seq(4,0,length.out=200))*sum(reg)/length(reg)*sum(p)/sum(weights),
                            weights=weights, maxit = 1e+8, thresh = 1e-5)   # "..." CAUSES AN ERROR!!! thresh reduced to allow convergence
    # Also increasing maxit from 1e+05 to 1e+08 
    # Also changing convergence threshold from 1e-07 to 1e-05, the same as in the recent java version
    # Also changing convergence threshold from 1e-07 to 1e-05, the same as in the recent java version
    class(model) <- c("maxnet", class(model))
    if (length(model$beta) < 200) stop("Error: glmnet failed to complete regularization path")
    bb <- model$beta[,200]
    model$betas <- bb[bb!=0]
    model$alpha <- 0
    rr <- predict.maxnet(model, data[p==0, , drop = FALSE], type="exponent", clamp=F)
    raw <- rr / sum(rr)
    model$entropy <- -sum(raw * log(raw))
    model$alpha <- -log(sum(rr))
    model$penalty.factor <- reg
    model$featuremins <- apply(mm, 2, min)
    model$featuremaxs <- apply(mm, 2, max)
    vv <- (sapply(data, class)!="factor")
    model$varmin <- apply(data[,vv, drop = FALSE], 2, min)
    model$varmax <- apply(data[,vv, drop = FALSE], 2, max)
    means <- apply(data[p==1,vv, drop = FALSE], 2, mean)
    majorities <- sapply(names(data)[!vv],
                         function(n) which.max(table(data[p==1,n, drop = FALSE])))
    names(majorities) <- names(data)[!vv]
    model$samplemeans <- unlist(c(means, majorities))
    model$levels <- lapply(data, levels)
    model
  }
#' 
#' 
#' #' @export
predict.maxnet <-
  function(object, newdata, clamp=F, type=c("link","exponential","cloglog","logistic"), ...) # our predict function DOES NOT clamp!
  {
    if (clamp) {
      for (v in intersect(names(object$varmax), names(newdata))) {
        newdata[,v] <- pmin(pmax(newdata[,v], object$varmin[v]), object$varmax[v])
      }
    }
    terms <- sub("hinge\\((.*)\\):(.*):(.*)$", "hingeval(\\1,\\2,\\3)", names(object$betas))
    terms <- sub("categorical\\((.*)\\):(.*)$", "categoricalval(\\1,\\2)", terms)
    terms <- sub("thresholds\\((.*)\\):(.*)$", "thresholdval(\\1,\\2)", terms)
    f <- formula(paste("~", paste(terms, collapse=" + "), "-1"))
    mm <- model.matrix(f, data.frame(newdata))
    if (clamp) mm <- t(pmin(pmax(t(mm), object$featuremins[names(object$betas)]),
                            object$featuremaxs[names(object$betas)]))
    link <- (mm %*% object$betas) + object$alpha
    type <- match.arg(type)
    if (type=="link") return(link)
    if (type=="exponential") return(exp(link))
    if (type=="cloglog") return(1-exp(0-exp(object$entropy+link)))
    if (type=="logistic") return(1/(1+exp(-object$entropy-link)))
  }
#' 
hingeval <-
  function(x, min, max)
  {
    pmin(1, pmax(0, (x-min)/(max-min)))
  }

require(maxnet)

## for whatever reason not found in namespace...
thresholdval <-
  function(x, knot) 
  {
    ifelse(x >= knot, 1, 0)
  }

#### DEBUG END

# put all species names in a vector
species <- c("Amphiesma_stol", "Boiga_krae", "Cyclophiops_majo", "Elaphe_cari",
             "Eutropis_long", "Japalura_swin", "Lycodon_rufo", "Lycodon_ruhs",
             "Oligodon_form", "Oreocryptophis_porp", "Ptyas_muco", "Trimeresurus_stej",
             "Xenochrophis_pisc")


num_bg <- 10000 ######Changed to 5000 for bugfixing only

NAcells <- which((is.na(my_PCs[[1]][]) == T))

VS_TRAIN <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

ES_TRAIN <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

for (i in (1:13)) {
  for (j in (1:10)) {
    VS_TRAIN[[i]][[j]] <- read.csv(file = paste0("VS/Locations/subsamples/", as.character(i), "/", as.character(j), "_po.csv"))
    ES_TRAIN[[i]][[j]] <- read.csv(file = paste0("ES/Locations/subsamples/occurrences_", species[i], "_", as.character(j), ".csv"))[,c(2,3)]
  }
}


###################
####  First VS  ###
###################

presence_matrix <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
presence_cn <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
samplefrom <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
bg_sample <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

maxnet_data <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

maxnet_p <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_models_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_models_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

lq_predictions_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
lq_predictions_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())


set.seed(2607)
for (i in 1:13) {
  for (j in 1:10) {
    
    presence_matrix[[i]][[j]] <- matrix(c(VS_TRAIN[[i]][[j]][,1], VS_TRAIN[[i]][[j]][,2]), ncol = 2)
    
    presence_cn[[i]][[j]] <- cellFromXY(object = my_PCs, xy = presence_matrix[[i]][[j]]) # occurrence cell numbers
    
    samplefrom[[i]][[j]] <- (1:ncell(my_PCs[[1]]))[!is.element((1:ncell(my_PCs[[1]]))[], c(NAcells, unlist(presence_cn[[i]])))]
    
    bg_sample[[i]][[j]] <- sample(x = samplefrom[[i]][[j]], size = num_bg, replace = F)
    
    maxnet_data[[i]][[j]] <- data.frame("presence" = c(rep(x = 1, times = NROW(VS_TRAIN[[i]][[j]])), rep(x = 0, times = num_bg)),
                                        "PC1" = c(as.vector(raster::extract(x = my_PCs[[1]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[1]], y = bg_sample[[i]][[j]]))),
                                        "PC2" = c(as.vector(raster::extract(x = my_PCs[[2]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[2]], y = bg_sample[[i]][[j]]))),
                                        "PC3" = c(as.vector(raster::extract(x = my_PCs[[3]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[3]], y = bg_sample[[i]][[j]]))),
                                        "PC4" = c(as.vector(raster::extract(x = my_PCs[[4]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[4]], y = bg_sample[[i]][[j]]))),
                                        "PC5" = c(as.vector(raster::extract(x = my_PCs[[5]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[5]], y = bg_sample[[i]][[j]]))),
                                        "PC6" = c(as.vector(raster::extract(x = my_PCs[[6]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[6]], y = bg_sample[[i]][[j]]))),
                                        "PC7" = c(as.vector(raster::extract(x = my_PCs[[7]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[7]], y = bg_sample[[i]][[j]]))),
                                        "PC8" = c(as.vector(raster::extract(x = my_PCs[[8]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[8]], y = bg_sample[[i]][[j]]))),
                                        "PC9" = c(as.vector(raster::extract(x = my_PCs[[9]], y = presence_matrix[[i]][[j]])),
                                                  as.vector(raster::extract(x = my_PCs[[9]], y = bg_sample[[i]][[j]])))
    )
    
    
    
    maxnet_p[[i]][[j]] <- maxnet_data[[i]][[j]]$presence
    
    lq_models_F9_0.5[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -1],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -1],
                                                            classes = "lq"),
                                         regmult = 0.5)
    
    lq_models_F9_1.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -1],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -1],
                                                            classes = "lq"),
                                         regmult = 1.0)
    
    lq_models_F9_4.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -1],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -1],
                                                            classes = "lq"),
                                         regmult = 4.0)
    
    lq_models_R7_0.5[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                                            classes = "lq"),
                                         regmult = 0.5)
    
    lq_models_R7_1.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                                            classes = "lq"),
                                         regmult = 1.0)
    
    lq_models_R7_4.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,9:10)],
                                                            classes = "lq"),
                                         regmult = 4.0)
    
    lq_models_R4_0.5[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                                            classes = "lq"),
                                         regmult = 0.5)
    
    lq_models_R4_1.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                                            classes = "lq"),
                                         regmult = 1.0)
    
    lq_models_R4_4.0[[i]][[j]] <- maxnet(p = maxnet_p[[i]][[j]], type = "cloglog", data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                         f = maxnet.formula(p = maxnet_p[[i]][[j]],
                                                            data = maxnet_data[[i]][[j]][, -c(1,6:10)],
                                                            classes = "lq"),
                                         regmult = 4.0)
    
    lq_predictions_F9_0.5[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                               model =  lq_models_F9_0.5[[i]][[j]],
                                               filename = paste0("VS/Models/F9/M0.5/", as.character(i), "/F9_05_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_F9_1.0[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                               model =  lq_models_F9_1.0[[i]][[j]],
                                               filename = paste0("VS/Models/F9/M1/", as.character(i), "/F9_10_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_F9_4.0[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                               model =  lq_models_F9_4.0[[i]][[j]],
                                               filename = paste0("VS/Models/F9/M4/", as.character(i), "/F9_40_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    
    lq_predictions_R7_0.5[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                               model =  lq_models_R7_0.5[[i]][[j]],
                                               filename = paste0("VS/Models/R7/M0.5/", as.character(i), "/R7_05_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_R7_1.0[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                               model =  lq_models_R7_1.0[[i]][[j]],
                                               filename = paste0("VS/Models/R7/M1/", as.character(i), "/R7_10_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_R7_4.0[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                               model =  lq_models_R7_4.0[[i]][[j]],
                                               filename = paste0("VS/Models/R7/M4/", as.character(i), "/R7_40_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_R4_0.5[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                               model =  lq_models_R4_0.5[[i]][[j]],
                                               filename = paste0("VS/Models/R4/M0.5/", as.character(i), "/R4_05_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_R4_1.0[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                               model =  lq_models_R4_1.0[[i]][[j]],
                                               filename = paste0("VS/Models/R4/M1/", as.character(i), "/R4_10_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    lq_predictions_R4_4.0[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                               model =  lq_models_R4_4.0[[i]][[j]],
                                               filename = paste0("VS/Models/R4/M4/", as.character(i), "/R4_40_", as.character(j), ".tif"),
                                               format = "GTiff", overwrite = T)
    
    if (j == 10) {
      mean_lq_predictions_F9_0.5 <- calc(x = stack(lq_predictions_F9_0.5[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/F9/M0.5/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_F9_1.0 <- calc(x = stack(lq_predictions_F9_1.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/F9/M1/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_F9_4.0 <- calc(x = stack(lq_predictions_F9_4.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/F9/M4/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      
      mean_lq_predictions_R7_0.5 <- calc(x = stack(lq_predictions_R7_0.5[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R7/M0.5/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_R7_1.0 <- calc(x = stack(lq_predictions_R7_1.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R7/M1/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_R7_4.0 <- calc(x = stack(lq_predictions_R7_4.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R7/M4/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      
      mean_lq_predictions_R4_0.5 <- calc(x = stack(lq_predictions_R4_0.5[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R4/M0.5/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_R4_1.0 <- calc(x = stack(lq_predictions_R4_1.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R4/M1/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
      mean_lq_predictions_R4_4.0 <- calc(x = stack(lq_predictions_R4_4.0[[i]]), fun = mean, na.rm = T,
                                         filename = paste0("VS/Models/R4/M4/", as.character(i), "/mean.tif"),
                                         format = "GTiff", overwrite = T)
    }
    
  }
}

###################
####  Evaluate  ###
###################
require(ENMTools)

SchoenersD <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                         "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

for (i in (1:13)) {
  SchoenersD[i, 2] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 3] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/F9/M1/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 4] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/F9/M4/", as.character(i), "/mean.tif")))$D
  
  SchoenersD[i, 5] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 6] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/R7/M1/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 7] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/R7/M4/", as.character(i), "/mean.tif")))$D
  
  SchoenersD[i, 8] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 9] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                               y = raster(x = paste0("VS/Models/R4/M1/", as.character(i), "/mean.tif")))$D
  SchoenersD[i, 10] <- ENMTools::raster.overlap(x = raster(x = paste0("VS/Truth/my_PA_", as.character(i), ".tif")),
                                                y = raster(x = paste0("VS/Models/R4/M4/", as.character(i), "/mean.tif")))$D
}

#######
# AUC
#######

AUCvs <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                    "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  AUCtemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    AUCtemp[j, 1] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 2] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 3] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 4] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 5] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 6] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 7] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 8] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 9] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    
    if (j == 10) {
      AUCvs[i, 2] <- mean(AUCtemp[, 1])
      AUCvs[i, 3] <- mean(AUCtemp[, 2])
      AUCvs[i, 4] <- mean(AUCtemp[, 3])
      AUCvs[i, 5] <- mean(AUCtemp[, 4])
      AUCvs[i, 6] <- mean(AUCtemp[, 5])
      AUCvs[i, 7] <- mean(AUCtemp[, 6])
      AUCvs[i, 8] <- mean(AUCtemp[, 7])
      AUCvs[i, 9] <- mean(AUCtemp[, 8])
      AUCvs[i, 10] <- mean(AUCtemp[, 9])
    }
  }
}

#######
# maxSSS
#######

maxSSSvs <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                       "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  maxSSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    maxSSStemp[j, 1] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 2] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 3] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 4] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 5] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 6] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 7] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 8] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 9] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    
    if (j == 10) {
      maxSSSvs[i, 2] <- mean(maxSSStemp[, 1])
      maxSSSvs[i, 3] <- mean(maxSSStemp[, 2])
      maxSSSvs[i, 4] <- mean(maxSSStemp[, 3])
      maxSSSvs[i, 5] <- mean(maxSSStemp[, 4])
      maxSSSvs[i, 6] <- mean(maxSSStemp[, 5])
      maxSSSvs[i, 7] <- mean(maxSSStemp[, 6])
      maxSSSvs[i, 8] <- mean(maxSSStemp[, 7])
      maxSSSvs[i, 9] <- mean(maxSSStemp[, 8])
      maxSSSvs[i, 10] <- mean(maxSSStemp[, 9])
    }
  }
}

#########################
#    compute weights    #
#########################

WEIGHTS_C <- raster::stack(my_PCs[[1]])

WEIGHTS_O <- list(raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]),
                  raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]),
                  raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]))


mean_radius_VS <- c()
mean_radius_VS[c(1,3,4,5,9,11,12)] <- 400 # operational scale rare species == radius
mean_radius_VS[c(2,6,7,8,10,13)] <- 400 # operational scale common species
# i want to change scale of common species to 500 and also change their number of cell in neighbours from 9 to 5

mean_weight_C <- list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
mean_weight_O <- list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())

PRESENCES_XY <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
PRESENCES_CN <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

SCE <- c("F9/M0.5/", "F9/M1/", "F9/M4/",
         "R7/M0.5/", "R7/M1/", "R7/M4/",
         "R4/M0.5/", "R4/M1/", "R4/M4/")

for (i in (1:13)) { #13
  PRESENCES_XY[[i]] <- read.csv(file = paste0("VS/Locations/PD_", as.character(i), "_po.csv"))
  PRESENCES_CN[[i]] <- unique(cellFromXY(object = my_PCs, xy = PRESENCES_XY[[i]]))
  set.seed(2607)
  bg <- sample(x = unique(unlist(bg_sample[[i]])), size = 10000, replace = F)
  WEIGHTS_C[[i]] <- distanceFromPoints(object = my_PCs[[1]], xy = PRESENCES_XY[[i]])
  WEIGHTS_C[[i]] <- WEIGHTS_C[[i]] / (mean_radius_VS[i] * 2)
  WEIGHTS_C[[i]][WEIGHTS_C[[i]] < 1] <-  1
  WEIGHTS_C[[i]] <- (1 - (1/WEIGHTS_C[[i]]))
  WEIGHTS_C <- raster::mask(x = WEIGHTS_C, mask = my_PCs[[1]])
  
  for (j in (1:9)) { #9
    tempprediction <- raster(x = paste0("VS/Models/", SCE[j], as.character(i), "/mean.tif"))
    
    tempweightC <- c()
    tempweightO <- c()
    
    WEIGHTS_O[[i]][[j]] <- focal(x = tempprediction, w = foc_weights, fun = my_foc_fun)
    WEIGHTS_O[[i]][[j]] <- WEIGHTS_O[[i]][[j]] - minValue(WEIGHTS_O[[i]][[j]])
    WEIGHTS_O[[i]][[j]] <- WEIGHTS_O[[i]][[j]] * 1/maxValue(WEIGHTS_O[[i]][[j]]) # scaling to 0..1
    
    COMMISSIONS_cells <- raster::extract(x = tempprediction, y = xyFromCell(object = tempprediction, cell = bg), cellnumbers = T)
    COMMISSIONS_cells <- COMMISSIONS_cells[COMMISSIONS_cells[,2] > maxSSSvs[i,j], 1]
    
    OMISSIONS_cells <- raster::extract(x = tempprediction, y = xyFromCell(object = tempprediction, cell = PRESENCES_CN[[i]]), cellnumbers = T)
    OMISSIONS_cells <- OMISSIONS_cells[OMISSIONS_cells[,2] < maxSSSvs[i,j], 1]
    
    if (length(COMMISSIONS_cells) != 0) {
      tempweightC <- raster::extract(x = WEIGHTS_C[[i]], y = COMMISSIONS_cells)
      (mean_weight_C[[i]][j] <- mean(tempweightC))
    }
    
    else {
      (mean_weight_C[[i]][j] <- 1) # weight irrelevant if there are no such errors, so we just set it to 1 to avoid NAs
    }
    
    if (length(OMISSIONS_cells) != 0) {
      tempweightO <- raster::extract(x = WEIGHTS_O[[i]][[j]], y = OMISSIONS_cells)
      (mean_weight_O[[i]][j] <- mean(tempweightO))
    }
    
    else {
      (mean_weight_O[[i]][j] <- 1) # weight irrelevant if there are no such errors, so we just set it to 1 to avoid NAs
    }
  }
  WEIGHTS_O[[i]] <- raster::mask(x = WEIGHTS_O[[i]], mask = my_PCs[[1]])
}

#######
# TSS
#######

TSSvs <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                    "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSvs[i, 2] <- mean(TSStemp[, 1])
      TSSvs[i, 3] <- mean(TSStemp[, 2])
      TSSvs[i, 4] <- mean(TSStemp[, 3])
      TSSvs[i, 5] <- mean(TSStemp[, 4])
      TSSvs[i, 6] <- mean(TSStemp[, 5])
      TSSvs[i, 7] <- mean(TSStemp[, 6])
      TSSvs[i, 8] <- mean(TSStemp[, 7])
      TSSvs[i, 9] <- mean(TSStemp[, 8])
      TSSvs[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS C
#######

TSSvsC <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                    "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSvsC[i, 2] <- mean(TSStemp[, 1])
      TSSvsC[i, 3] <- mean(TSStemp[, 2])
      TSSvsC[i, 4] <- mean(TSStemp[, 3])
      TSSvsC[i, 5] <- mean(TSStemp[, 4])
      TSSvsC[i, 6] <- mean(TSStemp[, 5])
      TSSvsC[i, 7] <- mean(TSStemp[, 6])
      TSSvsC[i, 8] <- mean(TSStemp[, 7])
      TSSvsC[i, 9] <- mean(TSStemp[, 8])
      TSSvsC[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS O
#######

TSSvsO <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSvsO[i, 2] <- mean(TSStemp[, 1])
      TSSvsO[i, 3] <- mean(TSStemp[, 2])
      TSSvsO[i, 4] <- mean(TSStemp[, 3])
      TSSvsO[i, 5] <- mean(TSStemp[, 4])
      TSSvsO[i, 6] <- mean(TSStemp[, 5])
      TSSvsO[i, 7] <- mean(TSStemp[, 6])
      TSSvsO[i, 8] <- mean(TSStemp[, 7])
      TSSvsO[i, 9] <- mean(TSStemp[, 8])
      TSSvsO[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS B
#######

TSSvsB <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][1])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][2])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][3])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][4])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][5])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][6])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][7])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][8])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][9])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSvsB[i, 2] <- mean(TSStemp[, 1])
      TSSvsB[i, 3] <- mean(TSStemp[, 2])
      TSSvsB[i, 4] <- mean(TSStemp[, 3])
      TSSvsB[i, 5] <- mean(TSStemp[, 4])
      TSSvsB[i, 6] <- mean(TSStemp[, 5])
      TSSvsB[i, 7] <- mean(TSStemp[, 6])
      TSSvsB[i, 8] <- mean(TSStemp[, 7])
      TSSvsB[i, 9] <- mean(TSStemp[, 8])
      TSSvsB[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# SEDI
#######

SEDIvs <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIvs[i, 2] <- mean(SEDItemp[, 1])
      SEDIvs[i, 3] <- mean(SEDItemp[, 2])
      SEDIvs[i, 4] <- mean(SEDItemp[, 3])
      SEDIvs[i, 5] <- mean(SEDItemp[, 4])
      SEDIvs[i, 6] <- mean(SEDItemp[, 5])
      SEDIvs[i, 7] <- mean(SEDItemp[, 6])
      SEDIvs[i, 8] <- mean(SEDItemp[, 7])
      SEDIvs[i, 9] <- mean(SEDItemp[, 8])
      SEDIvs[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI C
#######

SEDIvsC <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIvsC[i, 2] <- mean(SEDItemp[, 1])
      SEDIvsC[i, 3] <- mean(SEDItemp[, 2])
      SEDIvsC[i, 4] <- mean(SEDItemp[, 3])
      SEDIvsC[i, 5] <- mean(SEDItemp[, 4])
      SEDIvsC[i, 6] <- mean(SEDItemp[, 5])
      SEDIvsC[i, 7] <- mean(SEDItemp[, 6])
      SEDIvsC[i, 8] <- mean(SEDItemp[, 7])
      SEDIvsC[i, 9] <- mean(SEDItemp[, 8])
      SEDIvsC[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI O
#######

SEDIvsO <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                      "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIvsO[i, 2] <- mean(SEDItemp[, 1])
      SEDIvsO[i, 3] <- mean(SEDItemp[, 2])
      SEDIvsO[i, 4] <- mean(SEDItemp[, 3])
      SEDIvsO[i, 5] <- mean(SEDItemp[, 4])
      SEDIvsO[i, 6] <- mean(SEDItemp[, 5])
      SEDIvsO[i, 7] <- mean(SEDItemp[, 6])
      SEDIvsO[i, 8] <- mean(SEDItemp[, 7])
      SEDIvsO[i, 9] <- mean(SEDItemp[, 8])
      SEDIvsO[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI B
#######

SEDIvsB <- data.frame("VS" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                      "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][1])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][2])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][3])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][4])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][5])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][6])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][7])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][8])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("VS/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSvs[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_O[[i]][9])
    M[2,1] <- round(M[2,1] * mean_weight_C[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIvsB[i, 2] <- mean(SEDItemp[, 1])
      SEDIvsB[i, 3] <- mean(SEDItemp[, 2])
      SEDIvsB[i, 4] <- mean(SEDItemp[, 3])
      SEDIvsB[i, 5] <- mean(SEDItemp[, 4])
      SEDIvsB[i, 6] <- mean(SEDItemp[, 5])
      SEDIvsB[i, 7] <- mean(SEDItemp[, 6])
      SEDIvsB[i, 8] <- mean(SEDItemp[, 7])
      SEDIvsB[i, 9] <- mean(SEDItemp[, 8])
      SEDIvsB[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#########################
# SUMMARY of VS RESULTS #
#########################

VS_results_NO <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07",
                                                                     "VS08", "VS09", "VS10", "VS11", "VS12", "VS13"),
                                                                   c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                   c("TSS", "SEDI")))

VS_results_C <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07",
                                                                       "VS08", "VS09", "VS10", "VS11", "VS12", "VS13"),
                                                                     c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                     c("TSS", "SEDI")))

VS_results_O <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07",
                                                                       "VS08", "VS09", "VS10", "VS11", "VS12", "VS13"),
                                                                     c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                     c("TSS", "SEDI")))

VS_results_B <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("VS01", "VS02", "VS03", "VS04", "VS05", "VS06", "VS07",
                                                                       "VS08", "VS09", "VS10", "VS11", "VS12", "VS13"),
                                                                     c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                     c("TSS", "SEDI")))

VS_results_NO[,,1] <- as.matrix(TSSvs)[,2:10]
VS_results_NO[,,2] <- as.matrix(SEDIvs)[,2:10]

VS_results_C[,,1] <- as.matrix(TSSvsC)[,2:10]
VS_results_C[,,2] <- as.matrix(SEDIvsC)[,2:10]

VS_results_O[,,1] <- as.matrix(TSSvsO)[,2:10]
VS_results_O[,,2] <- as.matrix(SEDIvsO)[,2:10]

VS_results_B[,,1] <- as.matrix(TSSvsB)[,2:10]
VS_results_B[,,2] <- as.matrix(SEDIvsB)[,2:10]



######################################
###################
####  Then ES   ###
###################
######################################

Epresence_matrix <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Epresence_cn <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Esamplefrom <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Ebg_sample <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Emaxnet_data <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Emaxnet_p <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_models_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_models_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_F9_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_F9_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_F9_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_R7_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R7_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R7_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

Elqp_predictions_R4_0.5 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R4_1.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
Elqp_predictions_R4_4.0 <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

set.seed(2607)
for (i in 1:13) {
  for (j in 1:10) {
    
    Epresence_matrix[[i]][[j]] <- matrix(c(ES_TRAIN[[i]][[j]][,1], ES_TRAIN[[i]][[j]][,2]), ncol = 2)
    
    Epresence_cn[[i]][[j]] <- cellFromXY(object = my_PCs, xy = Epresence_matrix[[i]][[j]]) # occurrence cell numbers
    
    Esamplefrom[[i]][[j]] <- (1:ncell(my_PCs[[1]]))[!is.element((1:ncell(my_PCs[[1]]))[], c(NAcells, unique(unlist(Epresence_cn[[i]]))))]
    
    Ebg_sample[[i]][[j]] <- sample(x = Esamplefrom[[i]][[j]], size = num_bg, replace = F)
    
    Emaxnet_data[[i]][[j]] <- data.frame("presence" = c(rep(x = 1, times = NROW(ES_TRAIN[[i]][[j]])), rep(x = 0, times = num_bg)),
                                         "PC1" = c(as.vector(raster::extract(x = my_PCs[[1]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[1]], y = Ebg_sample[[i]][[j]]))),
                                         "PC2" = c(as.vector(raster::extract(x = my_PCs[[2]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[2]], y = Ebg_sample[[i]][[j]]))),
                                         "PC3" = c(as.vector(raster::extract(x = my_PCs[[3]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[3]], y = Ebg_sample[[i]][[j]]))),
                                         "PC4" = c(as.vector(raster::extract(x = my_PCs[[4]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[4]], y = Ebg_sample[[i]][[j]]))),
                                         "PC5" = c(as.vector(raster::extract(x = my_PCs[[5]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[5]], y = Ebg_sample[[i]][[j]]))),
                                         "PC6" = c(as.vector(raster::extract(x = my_PCs[[6]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[6]], y = Ebg_sample[[i]][[j]]))),
                                         "PC7" = c(as.vector(raster::extract(x = my_PCs[[7]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[7]], y = Ebg_sample[[i]][[j]]))),
                                         "PC8" = c(as.vector(raster::extract(x = my_PCs[[8]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[8]], y = Ebg_sample[[i]][[j]]))),
                                         "PC9" = c(as.vector(raster::extract(x = my_PCs[[9]], y = Epresence_matrix[[i]][[j]])),
                                                   as.vector(raster::extract(x = my_PCs[[9]], y = Ebg_sample[[i]][[j]])))
    )
    
    
    
    Emaxnet_p[[i]][[j]] <- Emaxnet_data[[i]][[j]]$presence
    
    Elqp_models_F9_0.5[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -1],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -1],
                                                              classes = "lqp"),
                                           regmult = 0.5)
    
    Elqp_models_F9_1.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -1],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -1],
                                                              classes = "lqp"),
                                           regmult = 1.0)
    
    Elqp_models_F9_4.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -1],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -1],
                                                              classes = "lqp"),
                                           regmult = 4.0)
    
    Elqp_models_R7_0.5[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                                              classes = "lqp"),
                                           regmult = 0.5)
    
    Elqp_models_R7_1.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                                              classes = "lqp"),
                                           regmult = 1.0)
    
    Elqp_models_R7_4.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,9:10)],
                                                              classes = "lqp"),
                                           regmult = 4.0)
    
    Elqp_models_R4_0.5[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                                              classes = "lqp"),
                                           regmult = 0.5)
    
    Elqp_models_R4_1.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                                              classes = "lqp"),
                                           regmult = 1.0)
    
    Elqp_models_R4_4.0[[i]][[j]] <- maxnet(p = Emaxnet_p[[i]][[j]], type = "cloglog", data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                           f = maxnet.formula(p = Emaxnet_p[[i]][[j]],
                                                              data = Emaxnet_data[[i]][[j]][, -c(1,6:10)],
                                                              classes = "lqp"),
                                           regmult = 4.0)
    
    Elqp_predictions_F9_0.5[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                                 model =  Elqp_models_F9_0.5[[i]][[j]],
                                                 filename = paste0("ES/Models/F9/M0.5/", as.character(i), "/F9_05_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_F9_1.0[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                                 model =  Elqp_models_F9_1.0[[i]][[j]],
                                                 filename = paste0("ES/Models/F9/M1/", as.character(i), "/F9_10_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_F9_4.0[[i]][[j]] <- predict(object = my_PCs, type = "cloglog",
                                                 model =  Elqp_models_F9_4.0[[i]][[j]],
                                                 filename = paste0("ES/Models/F9/M4/", as.character(i), "/F9_40_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    
    Elqp_predictions_R7_0.5[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                                 model =  Elqp_models_R7_0.5[[i]][[j]],
                                                 filename = paste0("ES/Models/R7/M0.5/", as.character(i), "/R7_05_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_R7_1.0[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                                 model =  Elqp_models_R7_1.0[[i]][[j]],
                                                 filename = paste0("ES/Models/R7/M1/", as.character(i), "/R7_10_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_R7_4.0[[i]][[j]] <- predict(object = my_PCs[[1:7]], type = "cloglog",
                                                 model =  Elqp_models_R7_4.0[[i]][[j]],
                                                 filename = paste0("ES/Models/R7/M4/", as.character(i), "/R7_40_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_R4_0.5[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                                 model =  Elqp_models_R4_0.5[[i]][[j]],
                                                 filename = paste0("ES/Models/R4/M0.5/", as.character(i), "/R4_05_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_R4_1.0[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                                 model =  Elqp_models_R4_1.0[[i]][[j]],
                                                 filename = paste0("ES/Models/R4/M1/", as.character(i), "/R4_10_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    Elqp_predictions_R4_4.0[[i]][[j]] <- predict(object = my_PCs[[1:4]], type = "cloglog",
                                                 model =  Elqp_models_R4_4.0[[i]][[j]],
                                                 filename = paste0("ES/Models/R4/M4/", as.character(i), "/R4_40_", as.character(j), ".tif"),
                                                 format = "GTiff", overwrite = T)
    
    if (j == 10) {
      mean_Elqp_predictions_F9_0.5 <- calc(x = stack(Elqp_predictions_F9_0.5[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/F9/M0.5/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_F9_1.0 <- calc(x = stack(Elqp_predictions_F9_1.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/F9/M1/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_F9_4.0 <- calc(x = stack(Elqp_predictions_F9_4.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/F9/M4/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      
      mean_Elqp_predictions_R7_0.5 <- calc(x = stack(Elqp_predictions_R7_0.5[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R7/M0.5/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_R7_1.0 <- calc(x = stack(Elqp_predictions_R7_1.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R7/M1/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_R7_4.0 <- calc(x = stack(Elqp_predictions_R7_4.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R7/M4/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      
      mean_Elqp_predictions_R4_0.5 <- calc(x = stack(Elqp_predictions_R4_0.5[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R4/M0.5/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_R4_1.0 <- calc(x = stack(Elqp_predictions_R4_1.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R4/M1/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
      mean_Elqp_predictions_R4_4.0 <- calc(x = stack(Elqp_predictions_R4_4.0[[i]]), fun = mean, na.rm = T,
                                           filename = paste0("ES/Models/R4/M4/", as.character(i), "/mean.tif"),
                                           format = "GTiff", overwrite = T)
    }
    
  }
}

###################
####  Evaluate  ###
###################

#######
# AUC
#######

AUCes <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                    "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  AUCtemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    AUCtemp[j, 1] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 2] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 3] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 4] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 5] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 6] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 7] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 8] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    AUCtemp[j, 9] <- SDMTools::auc(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))
    
    if (j == 10) {
      AUCes[i, 2] <- mean(AUCtemp[, 1])
      AUCes[i, 3] <- mean(AUCtemp[, 2])
      AUCes[i, 4] <- mean(AUCtemp[, 3])
      AUCes[i, 5] <- mean(AUCtemp[, 4])
      AUCes[i, 6] <- mean(AUCtemp[, 5])
      AUCes[i, 7] <- mean(AUCtemp[, 6])
      AUCes[i, 8] <- mean(AUCtemp[, 7])
      AUCes[i, 9] <- mean(AUCtemp[, 8])
      AUCes[i, 10] <- mean(AUCtemp[, 9])
    }
  }
}

#######
# maxSSS
#######

maxSSSes <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                       "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  maxSSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    maxSSStemp[j, 1] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 2] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 3] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 4] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 5] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 6] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 7] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 8] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    maxSSStemp[j, 9] <- SDMTools::optim.thresh(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)))[[5]][[1]]
    
    if (j == 10) {
      maxSSSes[i, 2] <- mean(maxSSStemp[, 1])
      maxSSSes[i, 3] <- mean(maxSSStemp[, 2])
      maxSSSes[i, 4] <- mean(maxSSStemp[, 3])
      maxSSSes[i, 5] <- mean(maxSSStemp[, 4])
      maxSSSes[i, 6] <- mean(maxSSStemp[, 5])
      maxSSSes[i, 7] <- mean(maxSSStemp[, 6])
      maxSSSes[i, 8] <- mean(maxSSStemp[, 7])
      maxSSSes[i, 9] <- mean(maxSSStemp[, 8])
      maxSSSes[i, 10] <- mean(maxSSStemp[, 9])
    }
  }
}

#########################
#    compute weights    #
#########################

WEIGHTS_E_C <- raster::stack(my_PCs[[1]])

WEIGHTS_E_O <- list(raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]),
                    raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]),
                    raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]), raster::stack(my_PCs[[1]]))


mean_radius_es <- c()
mean_radius_es[c(1,3,4,5,9,11,12)] <- 400 # operational scale rare species == radius
mean_radius_es[c(2,6,7,8,10,13)] <- 400 # operational scale common species
# i want to change scale of common species to 500 and also change their number of cell in neighbours from 9 to 5

mean_weight_e_C <- list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())
mean_weight_e_O <- list(c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c(),c())

PRESENCES_XY <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())
PRESENCES_CN <- list(list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list(), list())

SCE <- c("F9/M0.5/", "F9/M1/", "F9/M4/",
         "R7/M0.5/", "R7/M1/", "R7/M4/",
         "R4/M0.5/", "R4/M1/", "R4/M4/")

for (i in (1:13)) { #13
  PRESENCES_XY[[i]] <- read.csv(file = paste0("ES/Locations/subsamples/occurrences_", species[i], ".csv"))[,3:4]
  PRESENCES_CN[[i]] <- unique(cellFromXY(object = my_PCs, xy = PRESENCES_XY[[i]]))
  set.seed(2607)
  bg <- sample(x = unique(unlist(bg_sample[[i]])), size = 10000, replace = F)
  WEIGHTS_E_C[[i]] <- distanceFromPoints(object = my_PCs[[1]], xy = PRESENCES_XY[[i]])
  WEIGHTS_E_C[[i]] <- WEIGHTS_E_C[[i]] / (mean_radius_es[i] * 2)
  WEIGHTS_E_C[[i]][WEIGHTS_E_C[[i]] < 1] <-  1
  WEIGHTS_E_C[[i]] <- (1 - (1/WEIGHTS_E_C[[i]]))
  WEIGHTS_E_C <- raster::mask(x = WEIGHTS_E_C, mask = my_PCs[[1]])
  
  for (j in (1:9)) { #9
    tempprediction <- raster(x = paste0("ES/Models/", SCE[j], as.character(i), "/mean.tif"))
    
    tempweightC <- c()
    tempweightO <- c()
    
    WEIGHTS_E_O[[i]][[j]] <- focal(x = tempprediction, w = foc_weights, fun = my_foc_fun)
    WEIGHTS_E_O[[i]][[j]] <- WEIGHTS_E_O[[i]][[j]] - minValue(WEIGHTS_E_O[[i]][[j]])
    WEIGHTS_E_O[[i]][[j]] <- WEIGHTS_E_O[[i]][[j]] * 1/maxValue(WEIGHTS_E_O[[i]][[j]]) # scaling to 0..1
    
    COMMISSIONS_cells <- raster::extract(x = tempprediction, y = xyFromCell(object = tempprediction, cell = bg), cellnumbers = T)
    COMMISSIONS_cells <- COMMISSIONS_cells[COMMISSIONS_cells[,2] > maxSSSes[i,j], 1]
    
    OMISSIONS_cells <- raster::extract(x = tempprediction, y = xyFromCell(object = tempprediction, cell = PRESENCES_CN[[i]]), cellnumbers = T)
    OMISSIONS_cells <- OMISSIONS_cells[OMISSIONS_cells[,2] < maxSSSes[i,j], 1]
    
    if (length(COMMISSIONS_cells) != 0) {
      tempweightC <- raster::extract(x = WEIGHTS_E_C[[i]], y = COMMISSIONS_cells)
      (mean_weight_e_C[[i]][j] <- mean(tempweightC))
    }
    
    else {
      (mean_weight_e_C[[i]][j] <- 1) # weight irrelevant if there are no such errors, so we just set it to 1 to avoid NAs
    }
    
    if (length(OMISSIONS_cells) != 0) {
      tempweightO <- raster::extract(x = WEIGHTS_E_O[[i]][[j]], y = OMISSIONS_cells)
      (mean_weight_e_O[[i]][j] <- mean(tempweightO))
    }
    
    else {
      (mean_weight_e_O[[i]][j] <- 1) # weight irrelevant if there are no such errors, so we just set it to 1 to avoid NAs
    }
  }
  WEIGHTS_E_O[[i]] <- raster::mask(x = WEIGHTS_E_O[[i]], mask = my_PCs[[1]])
}

#######
# TSS
#######

TSSes <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                    "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSes[i, 2] <- mean(TSStemp[, 1])
      TSSes[i, 3] <- mean(TSStemp[, 2])
      TSSes[i, 4] <- mean(TSStemp[, 3])
      TSSes[i, 5] <- mean(TSStemp[, 4])
      TSSes[i, 6] <- mean(TSStemp[, 5])
      TSSes[i, 7] <- mean(TSStemp[, 6])
      TSSes[i, 8] <- mean(TSStemp[, 7])
      TSSes[i, 9] <- mean(TSStemp[, 8])
      TSSes[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS C
#######

TSSesC <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSesC[i, 2] <- mean(TSStemp[, 1])
      TSSesC[i, 3] <- mean(TSStemp[, 2])
      TSSesC[i, 4] <- mean(TSStemp[, 3])
      TSSesC[i, 5] <- mean(TSStemp[, 4])
      TSSesC[i, 6] <- mean(TSStemp[, 5])
      TSSesC[i, 7] <- mean(TSStemp[, 6])
      TSSesC[i, 8] <- mean(TSStemp[, 7])
      TSSesC[i, 9] <- mean(TSStemp[, 8])
      TSSesC[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS O
#######

TSSesO <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSesO[i, 2] <- mean(TSStemp[, 1])
      TSSesO[i, 3] <- mean(TSStemp[, 2])
      TSSesO[i, 4] <- mean(TSStemp[, 3])
      TSSesO[i, 5] <- mean(TSStemp[, 4])
      TSSesO[i, 6] <- mean(TSStemp[, 5])
      TSSesO[i, 7] <- mean(TSStemp[, 6])
      TSSesO[i, 8] <- mean(TSStemp[, 7])
      TSSesO[i, 9] <- mean(TSStemp[, 8])
      TSSesO[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# TSS B
#######

TSSesB <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  TSStemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    pres <- presence_cn[[i]][[j]]
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][1])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][1])
    TSStemp[j, 1] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][2])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][2])
    TSStemp[j, 2] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][3])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][3])
    TSStemp[j, 3] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][4])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][4])
    TSStemp[j, 4] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][5])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][5])
    TSStemp[j, 5] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][6])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][6])
    TSStemp[j, 6] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][7])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][7])
    TSStemp[j, 7] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][8])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][8])
    TSStemp[j, 8] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][9])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][9])
    TSStemp[j, 9] <- TSS(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])
    
    if (j == 10) {
      TSSesB[i, 2] <- mean(TSStemp[, 1])
      TSSesB[i, 3] <- mean(TSStemp[, 2])
      TSSesB[i, 4] <- mean(TSStemp[, 3])
      TSSesB[i, 5] <- mean(TSStemp[, 4])
      TSSesB[i, 6] <- mean(TSStemp[, 5])
      TSSesB[i, 7] <- mean(TSStemp[, 6])
      TSSesB[i, 8] <- mean(TSStemp[, 7])
      TSSesB[i, 9] <- mean(TSStemp[, 8])
      TSSesB[i, 10] <- mean(TSStemp[, 9])
    }
  }
}

#######
# SEDI
#######

SEDIes <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                     "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIes[i, 2] <- mean(SEDItemp[, 1])
      SEDIes[i, 3] <- mean(SEDItemp[, 2])
      SEDIes[i, 4] <- mean(SEDItemp[, 3])
      SEDIes[i, 5] <- mean(SEDItemp[, 4])
      SEDIes[i, 6] <- mean(SEDItemp[, 5])
      SEDIes[i, 7] <- mean(SEDItemp[, 6])
      SEDIes[i, 8] <- mean(SEDItemp[, 7])
      SEDIes[i, 9] <- mean(SEDItemp[, 8])
      SEDIes[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI C
#######

SEDIesC <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                      "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIesC[i, 2] <- mean(SEDItemp[, 1])
      SEDIesC[i, 3] <- mean(SEDItemp[, 2])
      SEDIesC[i, 4] <- mean(SEDItemp[, 3])
      SEDIesC[i, 5] <- mean(SEDItemp[, 4])
      SEDIesC[i, 6] <- mean(SEDItemp[, 5])
      SEDIesC[i, 7] <- mean(SEDItemp[, 6])
      SEDIesC[i, 8] <- mean(SEDItemp[, 7])
      SEDIesC[i, 9] <- mean(SEDItemp[, 8])
      SEDIesC[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI O
#######

SEDIesO <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                      "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIesO[i, 2] <- mean(SEDItemp[, 1])
      SEDIesO[i, 3] <- mean(SEDItemp[, 2])
      SEDIesO[i, 4] <- mean(SEDItemp[, 3])
      SEDIesO[i, 5] <- mean(SEDItemp[, 4])
      SEDIesO[i, 6] <- mean(SEDItemp[, 5])
      SEDIesO[i, 7] <- mean(SEDItemp[, 6])
      SEDIesO[i, 8] <- mean(SEDItemp[, 7])
      SEDIesO[i, 9] <- mean(SEDItemp[, 8])
      SEDIesO[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#######
# SEDI B
#######

SEDIesB <- data.frame("ES" = c(1,2,3,4,5,6,7,8,9,10,11,12,13), "F9_M0.5" = NA,  "F9_M1.0" = NA, "F9_M4.0" = NA,
                      "R7_M0.5" = NA, "R7_M1.0" = NA, "R7_M4.0" = NA, "R4_M0.5" = NA, "R4_M1.0" = NA, "R4_M4.0" = NA)

# obs <- c(rep(x = 1, times = 1000), rep(x = 0, times = 10000))
obs <- c(rep(x = 1, times = 100), rep(x = 0, times = 10000))

for (i in (1:13)) { # for each species
  SEDItemp <- matrix(0, nrow = 10, ncol = 9)
  
  for (j in (1:10)) {
    
    pres <- presence_cn[[i]][[j]]
    
    bg <- bg_sample[[i]][[j]]
    
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 2])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][1])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][1])
    SEDItemp[j, 1] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 3])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][2])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][2])
    SEDItemp[j, 2] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/F9/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 4])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][3])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][3])
    SEDItemp[j, 3] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 5])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][4])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][4])
    SEDItemp[j, 4] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 6])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][5])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][5])
    SEDItemp[j, 5] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R7/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 7])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][6])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][6])
    SEDItemp[j, 6] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M0.5/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 8])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][7])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][7])
    SEDItemp[j, 7] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M1/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 9])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][8])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][8])
    SEDItemp[j, 8] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    M <- SDMTools::confusion.matrix(obs = obs, pred = extract(x = raster(x = paste0("ES/Models/R4/M4/", as.character(i),"/mean.tif")), y = c(pres, bg)), threshold = maxSSSes[i, 10])
    M[1,2] <- round(M[1,2] * mean_weight_e_O[[i]][9])
    M[2,1] <- round(M[2,1] * mean_weight_e_C[[i]][9])
    SEDItemp[j, 9] <- sedi(a = M[2,2], b = M[2,1], c = M[1,2], d = M[1,1])[[2]]
    
    if (j == 10) {
      SEDIesB[i, 2] <- mean(SEDItemp[, 1])
      SEDIesB[i, 3] <- mean(SEDItemp[, 2])
      SEDIesB[i, 4] <- mean(SEDItemp[, 3])
      SEDIesB[i, 5] <- mean(SEDItemp[, 4])
      SEDIesB[i, 6] <- mean(SEDItemp[, 5])
      SEDIesB[i, 7] <- mean(SEDItemp[, 6])
      SEDIesB[i, 8] <- mean(SEDItemp[, 7])
      SEDIesB[i, 9] <- mean(SEDItemp[, 8])
      SEDIesB[i, 10] <- mean(SEDItemp[, 9])
    }
  }
}

#########################
# SUMMARY of ES RESULTS #
#########################

ES_results_NO <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("ES01", "ES02", "ES03", "ES04", "ES05", "ES06", "ES07",
                                                                              "ES08", "ES09", "ES10", "ES11", "ES12", "ES13"),
                                                                            c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                            c("TSS", "SEDI")))

ES_results_C <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("ES01", "ES02", "ES03", "ES04", "ES05", "ES06", "ES07",
                                                                             "ES08", "ES09", "ES10", "ES11", "ES12", "ES13"),
                                                                           c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                           c("TSS", "SEDI")))

ES_results_O <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("ES01", "ES02", "ES03", "ES04", "ES05", "ES06", "ES07",
                                                                             "ES08", "ES09", "ES10", "ES11", "ES12", "ES13"),
                                                                           c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                           c("TSS", "SEDI")))

ES_results_B <- array(data = numeric(), dim = c(13, 9, 2), dimnames = list(c("ES01", "ES02", "ES03", "ES04", "ES05", "ES06", "ES07",
                                                                             "ES08", "ES09", "ES10", "ES11", "ES12", "ES13"),
                                                                           c("F9_M0.5", "F9_M1.0", "F9_M4.0", "R7_M0.5", "R7_M1.0", "R7_M4.0", "R4_M0.5", "R4_M1.0", "R4_M4.0"),
                                                                           c("TSS", "SEDI")))

ES_results_NO[,,1] <- as.matrix(TSSes)[,2:10]
ES_results_NO[,,2] <- as.matrix(SEDIes)[,2:10]

ES_results_C[,,1] <- as.matrix(TSSesC)[,2:10]
ES_results_C[,,2] <- as.matrix(SEDIesC)[,2:10]

ES_results_O[,,1] <- as.matrix(TSSesO)[,2:10]
ES_results_O[,,2] <- as.matrix(SEDIesO)[,2:10]

ES_results_B[,,1] <- as.matrix(TSSesB)[,2:10]
ES_results_B[,,2] <- as.matrix(SEDIesB)[,2:10]
