### Original idea conceived by Fielding and Bell (1997)
### Pseudo-code and R code developed by Rainer Wunderlich with comments from Yu-Pin Lin, Johnathen Anthony
### Manuscript in preparation

require(raster)
rasterOptions(tmpdir = "/home/affu/R_TMP/") # avoid filling /
rasterOptions(tolerance = 0)
rasterOptions(progress = "text")
require(rgdal)
require(rgeos)
require(SDMTools)
# require(spacetime) #to use environmental distances

###############################################################
# define EPSG:6933
###############################################################
epsg6933 <- CRS("+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +ellps=WGS84 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs")


###############################################################
# read a single raster and determine cell indices of NA
###############################################################
my_result <- raster(x = "../Roadkill/Predictions/Mean_prediction_Trimeresurus_stej.tif")
NA_cells <- Which(x = is.na(my_result) == T, cells = T)

###############################################################
# read occurrences, we need x, y, radius, the latter being the
# larger of home range radius and positional uncertainty
# compute median radius
# cell indices of presences 
###############################################################
occurrences <- read.csv(file = "../Roadkill/Occurrences/RoadKill_EPSG_6933_NoNAs.csv", header = T, sep = ",", quote = "\"", dec = ".")
occurrences_Tri_ste <- unique(occurrences[occurrences$species == "Trimeresurus stejnegeri",]) # unique removes exact duplicates
set.seed(2607)
occurrences_Tri_ste$radius <- rnorm(mean = 250, sd =  3, n = NROW(occurrences_Tri_ste))
XY_all <- matrix(cbind(occurrences_Tri_ste$x..EPSG_69, occurrences_Tri_ste$y..EPSG_69), ncol = 2)
Cell_numbers_all <- cellFromXY(object = my_result, xy = XY_all)
occurrences_Tri_ste$cell <- Cell_numbers_all

duplicates <- which(duplicated(Cell_numbers_all) == T)
Cell_numbers_duplicated <- Cell_numbers_all[duplicates]
Cell_numbers_unique <- Cell_numbers_all[-duplicates]

occurrences_Tri_ste_unique <- occurrences_Tri_ste[is.element(occurrences_Tri_ste$cell, Cell_numbers_duplicated) == F, ]
XY_unique <- matrix(cbind(occurrences_Tri_ste_unique$x..EPSG_69, occurrences_Tri_ste_unique$y..EPSG_69), ncol = 2)

median_radius <- round(median(x = occurrences_Tri_ste$radius))
###############################################################
# calculate distance raster, reuse XY to fasten computation
###############################################################
distRaster <- distanceFromPoints(object = my_result, xy = XY_all)
spplot(distRaster)

###############################################################
# mask distRaster with prediction and divide by median_radius
# plot for visual control
###############################################################
distRaster <- mask(x = distRaster, mask = my_result)
distRaster <- distRaster / median_radius
spplot(distRaster)

###############################################################
# scale result to 0.x .. y, with x * y = 1 to maintain symmetry
# around 1
# suggest a default of 0.2 and 5 or 0.1 and 10
# plot for visual control
###############################################################
distRaster[distRaster <= 1] <-  0 # We don't care AT ALL if we 'wrongly' predict into an actual home range or within the range of uncertainty
distRaster[(distRaster > 1) * (distRaster <= 2)] <- 2/3 # We care a bit if we are a bit off
distRaster[(distRaster > 2) * (distRaster <= 3)] <- 1 # We care a bit more if we are a bit more off
distRaster[distRaster > 3] <- 4/3 # We do care yet more if were are even more off
spplot(distRaster)

###############################################################
# calculate maxSSS and confusion matrix and AUC and load sedi.R
###############################################################
NA_Pres_cells <- unique(c(NA_cells, Cell_numbers_unique))
AvailableForBG_cells <- setdiff(1:ncell(my_result), NA_Pres_cells)

obs <- c(rep(1, NROW(x = occurrences_Tri_ste_unique)), rep(0, 5000))
pred <- c(extract(x = my_result, y = XY_unique),
          extract(x = my_result, y = sample(x = AvailableForBG_cells, size = 5000)))
(maxSSS <- optim.thresh(obs = obs, pred = pred)[5][1])
(CM_1 <- confusion.matrix(obs = obs, pred = pred, threshold = as.vector(unlist(maxSSS))))
(AUC <- auc(obs = obs, pred = pred))
source("../Code/R.Code/sedi.R")

###############################################################
# evaluate #1 (TSS, SEDI)
###############################################################
(TSS_1 <- TSS(a = CM_1[1,1], b = CM_1[1,2], c = CM_1[2,1], d = CM_1[2,2]))
(SEDI_1 <- sedi(a = CM_1[1,1], b = CM_1[1,2], c = CM_1[2,1], d = CM_1[2,2]))

###############################################################
# extract scaled weight for all commissions and compute mean
###############################################################
COMMISSIONS_pred <- extract(x = my_result, y = XY_unique)
COMM_index <- which(x = COMMISSIONS_pred < as.vector(unlist(maxSSS)))
COMMISSIONS_pred <- COMMISSIONS_pred[COMM_index]
COMMISSIONS_dist <- extract(x = distRaster, y = XY_unique, cellnumbers = T)
COMMISSIONS_dist <- COMMISSIONS_dist[COMM_index, ]

median_weight <- median(COMMISSIONS_dist[,2])

###############################################################
# multiply 'b' in the confusion matrix with the mean weight
###############################################################
CM_2 <- CM_1
(CM_2[1,2] <- round(CM_2[1,2] * median_weight))
CM_2

###############################################################
# evaluate #2 (TSS, SEDI)
###############################################################
(TSS_2 <- TSS(a = CM_2[1,1], b = CM_2[1,2], c = CM_2[2,1], d = CM_2[2,2]))
(SEDI_2 <- sedi(a = CM_2[1,1], b = CM_2[1,2], c = CM_2[2,1], d = CM_2[2,2]))

###############################################################
# COMPARE!!
###############################################################
