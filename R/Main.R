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
spplot(my_result)
NA_cells <- Which(x = is.na(my_result) == T, cells = T)

###############################################################
# read occurrences, we need x, y, radius, the latter being the
# larger of home range radius and positional uncertainty
# compute mean radius
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

mean_radius <- round(mean(x = occurrences_Tri_ste$radius))



###############################################################
# calculate maxSSS and confusion matrix and AUC and load sedi.R
###############################################################
NA_Pres_cells <- unique(c(NA_cells, Cell_numbers_unique))
AvailableForBG_cells <- setdiff(1:ncell(my_result), NA_Pres_cells)
set.seed(2607)
BG <- sample(x = AvailableForBG_cells, size = 5000)

obs <- c(rep(1, length(x = Cell_numbers_unique)), rep(0, length(BG)))
pred <- c(extract(x = my_result, y = Cell_numbers_unique),
          extract(x = my_result, y = BG))
(maxSSS <- optim.thresh(obs = obs, pred = pred)[5][[1]])
(CM_1 <- confusion.matrix(obs = obs, pred = pred, threshold = as.vector(unlist(maxSSS))))
source("../Code/R.Code/sedi.R")

###############################################################
# evaluate AUC
###############################################################
(AUCscore <- SDMTools::auc(obs = obs, pred = pred))

###############################################################
# evaluate #1 (TSS, SEDI)
###############################################################
(TSS_1 <- TSS(a = CM_1[2,2], b = CM_1[2,1], c = CM_1[1,2], d = CM_1[1,1]))
(SEDI_1 <- sedi(a = CM_1[2,2], b = CM_1[2,1], c = CM_1[1,2], d = CM_1[1,1])[[2]])

###############################################################
# calculate distance raster
###############################################################
distRaster <- distanceFromPoints(object = my_result, xy = XY_all) # cell number duplicates are not relevant in this context
spplot(distRaster)

###############################################################
# mask distRaster with prediction and divide by mean_diameter
# plot for visual control
###############################################################
distRaster <- mask(x = distRaster, mask = my_result)
distRaster <- distRaster / (mean_radius * 2)
spplot(distRaster)

###############################################################
# closely following Fielding & Bell 1997
# but we also consider omission weights based on neighboorhood predictions,
# if we miss entirely, i.e. dont say yes anywhere in the neighborhood, it is VERY bad
###############################################################

Prev <- length(Which(x = my_result > maxSSS[[1]], cells = T)) / length(my_result[is.na(my_result) == F])

distRaster[distRaster < 1] <-  1
spplot(distRaster)
weightRaster_C <- (1 - (1/distRaster))

foc_weights <- matrix(data = 1, nrow = 3, ncol = 3)
my_foc_fun <- function(x, y = maxSSS[[1]], na.rm = T) {
  CELLS_nh <- length(which(x >= 0))
  CELLS_p <- length(which(x > y)) # no need to differentiaate between TRUE FALSE here
  RES <- (CELLS_nh + 1) / (CELLS_p + 1)
  return(RES)
}

weightRaster_O <- focal(x = my_result, w = foc_weights, fun = my_foc_fun)
weightRaster_O <- mask(x = weightRaster_O, mask = my_result)
weightRaster_O <- weightRaster_O - minValue(weightRaster_O)
weightRaster_O <- weightRaster_O * 1/maxValue(weightRaster_O) # scaling to 0..1
weightRaster_O <- weightRaster_O * 1/Prev # taking into account prevalence: if very low omissions are more serious

spplot(distRaster)
spplot(weightRaster_C)
spplot(weightRaster_O)


###############################################################
# extract scaled weight for all commissions/omissions and compute mean
###############################################################
COMMISSIONS_pred <- extract(x = my_result, y = xyFromCell(object = my_result, cell = BG), cellnumbers = T)
COMMISSIONS_cells <- COMMISSIONS_pred[COMMISSIONS_pred[,2] > maxSSS, 1]

COMMISSIONS_weight <- extract(x = weightRaster_C, y = COMMISSIONS_cells)
(mean_weight_C <- mean(COMMISSIONS_weight))

OMISSIONS_pred <- extract(x = my_result, y = xyFromCell(object = my_result, cell = Cell_numbers_unique), cellnumbers = T)
OMISSIONS_cells <- OMISSIONS_pred[OMISSIONS_pred[,2] < maxSSS, 1]

OMISSIONS_weight <- extract(x = weightRaster_O, y = OMISSIONS_cells)
(mean_weight_O <- mean(OMISSIONS_weight))

## ONLY CONSIDER COMMISSION WEIGHTS:

###############################################################
# multiply 'b' in the confusion matrix with the mean weight
###############################################################
CM_2 <- CM_1
(CM_2[2,1] <- round(CM_2[2,1] * mean_weight_C))
CM_2

###############################################################
# evaluate #2 (TSS, SEDI)
###############################################################
(TSS_2 <- TSS(a = CM_2[2,2], b = CM_2[2,1], c = CM_2[1,2], d = CM_2[1,1]))
(SEDI_2 <- sedi(a = CM_2[2,2], b = CM_2[2,1], c = CM_2[1,2], d = CM_2[1,1])[[2]])

## ONLY CONSIDER OMISSION WEIGHTS:

###############################################################
# multiply 'c' in the confusion matrix with the mean weight
###############################################################
CM_3 <- CM_1
(CM_3[1,2] <- round(CM_3[1,2] * mean_weight_O))
CM_3

###############################################################
# evaluate #3 (TSS, SEDI)
###############################################################
(TSS_3 <- TSS(a = CM_3[2,2], b = CM_3[2,1], c = CM_3[1,2], d = CM_3[1,1]))
(SEDI_3 <- sedi(a = CM_3[2,2], b = CM_3[2,1], c = CM_3[1,2], d = CM_3[1,1])[[2]])

## CONSIDER BOTH ERROR TYPE WEIGHTS

###############################################################
# multiply 'b' and 'c' in the confusion matrix with the mean weights
###############################################################
CM_4 <- CM_1
(CM_4[2,1] <- round(CM_4[2,1] * mean_weight_C))
(CM_4[1,2] <- round(CM_4[1,2] * mean_weight_O))
CM_4

###############################################################
# evaluate #4 (TSS, SEDI)
###############################################################
(TSS_4 <- TSS(a = CM_4[2,2], b = CM_4[2,1], c = CM_4[1,2], d = CM_4[1,1]))
(SEDI_4 <- sedi(a = CM_4[2,2], b = CM_4[2,1], c = CM_4[1,2], d = CM_4[1,1])[[2]])
