# Script to perform species' distribution models. The outputs of this file are used in script 01A_organising-data. 

# This includes 
# 1. Combining and standardising covariate raster data
# 2. Functions to fit SDMs 
# 3. Fitting of SDMs (suggest use in parallel )

# Initiated 08/10/2018
# Author: Conor Waldock

# This script will not run fully without obtaining all the necessary raster files and shape polygons and storing in an appropriate folder. 
# Then all file paths to the raster files must be changed within this script. 
# Script will run if loading RasterStackV4.grd which has the necessary rasters post-processing - but will need to obtain MEOW shapefiles. 

# PRE-AMBLE ----
# Load libraries for spatial analysis ----
library(sp)
library(raster)
library(rgdal)
library(sdm)
library(rgeos)
library(RStoolbox)
library(tryCatchLog)
library(sdmvspecies)
library(doParallel)
library(filesstrings)
library(dplyr)
library(maptools)
library(plyr)

# Custom functions ----

# This set of functions  replaces all NAs that are found in the raster with a nearest neighbour. Very slow. 
sample_raster_NA <- function(r, xy){
  apply(X = xy, MARGIN = 1, 
        FUN = function(xy) r@data@values[which.min(replace(distanceFromPoints(r, xy), is.na(r), NA))])
}
FillSiteNAs <- function(StackLayer, # Environmental stack layer
                        RLS_Sites   # Coordiantes as a spatial points database.  
){
  
  # Do first extraction that finds the NA values inside the subset (extract(StackLayer, RLS_Sites))
  # Find the coordinates at the values
  xy <- RLS_sites_SPDF@coords[is.na(extract(StackLayer, RLS_Sites)),]
  
  # Match to the nearest raster cell using a pre-defined function 
  NA_coords <- sample_raster_NA(StackLayer, xy)
  StackLayer[cellFromXY(StackLayer, as.matrix(xy))] <- NA_coords
  return(StackLayer)
}

# ----

# SET UP ENVIRONMENTAL DATA ----
# Load in all rasters and combine into a stack ----
# Get rasters
path  <- '/Volumes/Untitled/Raster_files/BioOricleV2/'
lst   <- list.files(path=path,pattern='asc$',full.names = T)
path2 <- '/Volumes/Untitled/Raster_files/MSECData_Yeager2017'
lst2  <- list.files(path=path2,pattern='nc$',full.names = T)[c(3)]#,4,8,12,13,14)]
lst3  <- list.files(path=path2,pattern='nc$',full.names = T)[c(5,8,12,13,14)]

path3 <- '/Volumes/Untitled/Raster_files/GEBCO-30Sec'
lst4  <- list.files(path=path3,pattern='nc$',full.names = T)

# World ocean mask. 
install.packages('spaMM')
library(spaMM)
data("oceanmask")
plot(oceanmask)

preds_BIOORCLE  <- stack(lst)
preds_MSEChuman <- stack(lst2)
preds_MSEChuman <- preds_MSEChuman$X2015
preds_MSEC      <- stack(lst3)
preds_MSEChuman <- rotate(preds_MSEChuman)
preds_MSEC      <- rotate(preds_MSEC)
extent(preds_MSEC)       <- extent(preds_BIOORCLE)
extent(preds_MSEChuman)  <- extent(preds_BIOORCLE)
crs(preds_BIOORCLE)      <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
crs(preds_MSEC)          <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'
crs(preds_MSEChuman)     <- '+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0'

# Handle depth data
preds_DEPTH     <- raster(lst4)
preds_DEPTH_reprojected <- projectRaster(preds_DEPTH, preds_BIOORCLE[[1]])
#preds_DEPTH <- mask(preds_DEPTH_reprojected, oceanmask)
# preds_DEPTH_reprojected[preds_DEPTH_reprojected > 100] <- NA
plot(preds_DEPTH_reprojected)


# Stack all the environmental data in one place for use in loop. 
preds_MSEC_Crop_reprojected      <- projectRaster(preds_MSEC, preds_BIOORCLE[[1]])
preds_MSEChuman_Crop_reprojected <- projectRaster(preds_MSEChuman, preds_BIOORCLE[[1]])

# Stack the re-scaled rasters together. 
preds_Stack1 <- stack(preds_MSEC_Crop_reprojected, preds_BIOORCLE)
preds_Stack2 <- stack(preds_MSEChuman_Crop_reprojected, preds_Stack1)
preds_Stack3 <- stack(preds_DEPTH_reprojected, preds_Stack2)
plot(preds_Stack3)

# Drop correlated layers
preds_Stack3 <- dropLayer(preds_Stack3, c(6, 9, 11, 16))
names(preds_Stack3) <- c('Depth', 'HumanPop', 'LandArea', 'NPP', 'ReefArea', 'WaveEnergy', 'CurrentVelocity',
                         'O2', 'Iron', 'Nitrate', 'pH', 'Phosphate', 'Salinity', 'Silicate', 'Temperature')

hist(preds_Stack3$Depth)
hist(log10(preds_Stack3$HumanPop+1))
hist(log10(preds_Stack3$LandArea + 1))
hist(log10(preds_Stack3$NPP + 1))
hist(log10(preds_Stack3$ReefArea + 1))
hist(log10(preds_Stack3$WaveEnergy + 1))
hist(log10(preds_Stack3$CurrentVelocity+1))
hist(preds_Stack3$O2)
hist(log10(preds_Stack3$Iron+0.000001))
hist(preds_Stack3$Nitrate)
hist(preds_Stack3$pH)
hist(preds_Stack3$Phosphate)
hist(preds_Stack3$Salinity)
hist(preds_Stack3$Silicate)
hist(preds_Stack3$Temperature)

# Apply transformations to data layers 
?calc
preds_Stack3$HumanPop <- calc(preds_Stack3$HumanPop, function(x) log10(x+1))
preds_Stack3$LandArea <- calc(preds_Stack3$LandArea, function(x) log10(x+1))
preds_Stack3$NPP <- calc(preds_Stack3$NPP, function(x) log10(x+1))
preds_Stack3$ReefArea <- calc(preds_Stack3$ReefArea, function(x) log10(x+1))
preds_Stack3$WaveEnergy <- calc(preds_Stack3$WaveEnergy, function(x) log10(x+1))
preds_Stack3$CurrentVelocity <- calc(preds_Stack3$CurrentVelocity, function(x) log10(x+1))
preds_Stack3$Iron <- calc(preds_Stack3$Iron, function(x) log10(x+0.000001))



# This type of file save means it can be read in by another computer without the file links. 
#writeRaster(preds_Stack3, 
#            file = '/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack_V2.grd', 
#            options="INTERLEAVE=BAND", 
#            overwrite=TRUE)

preds_Stack3 <- stack('/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack_V2.grd')
png(file = 'figures_extra/depthMap.png', res=600, height = 5000, width = 10000)
preds_Stack3$Depth[preds_Stack3$Depth < -200] <- NA
image(preds_Stack3$Depth)
dev.off()

# Read in reef shapefiles to buffer environmental data so that not extrapolating into pelagic ocean ----

# Read in reef spatial polygon dataset
ReefPoly <- readOGR(dsn = '/Volumes/Untitled/Raster_files/WCMC_ReefData/01_Data')
#plot(ReefPoly[1,])
#p1 <- ReefPoly[1,]
#p1_100 <- ReefPoly[1:1000,]

# Create simpler reef shapes as too complex to manage at a global scale. 
#ReefPoly_simple <- gSimplify(ReefPoly, tol=0.1, topologyPreserve=T)

# Create a subset of smaller polygons and see if these get stuck in buffer too. 
Areas_Poly <- do.call(rbind, lapply(ReefPoly@polygons, function(x){x@area}))
SubsetOut <- round(0.9*length(Areas_Poly)):length(Areas_Poly)

# TESTING
#ReefPoly_simple_small <- ReefPoly_simple[order(Areas_Poly)[-SubsetOut],]
#result <- list()
#for(i in 1:length(ReefPoly_simple_small)){result[[i]] <- gBuffer(ReefPoly_simple_small[i,], byid = T, width = 0.5); print(i/length(ReefPoly_simple_small))}
#result <- do.call("rbind", result)

# Function to disaggregate large polygons as causes issues with memory. 
ConvertPolygonToBuffer <- function(Polygon, ID, i){
  Polygon_disaggregated <- sp::disaggregate(Polygon)
  Polygon_disaggregated_simple <- gSimplify(Polygon_disaggregated, tol=0.1, topologyPreserve=T)
  Polygon_disaggregated_simple_Buffer <- gBuffer(Polygon_disaggregated_simple, byid = T, width = 0.5)
  Polygon_disaggregated_simple_Buffer_aggregate <- aggregate(Polygon_disaggregated_simple_Buffer)
  Polygon_disaggregated_simple_Buffer_aggregate <- SpatialPolygonsDataFrame(Polygon_disaggregated_simple_Buffer_aggregate, data = data.frame(ID=ID[i]))
  return(Polygon_disaggregated_simple_Buffer_aggregate)
}

# Run buffer function over all .
result <- list()
SmallPolygonIDs <- order(Areas_Poly)[-SubsetOut]
for(i in 1:length(ReefPoly)){
  result[[i]] <- ConvertPolygonToBuffer(ReefPoly[i,], ID = SmallPolygonIDs, i=i); 
  print(i/length(ReefPoly))}
result <- do.call("rbind", result)

# Run buffer function over each polygon (and feature set). 
result_BIG <- list()
BigPolygonIDs <- order(Areas_Poly)[SubsetOut]
for(i in 1:length(ReefPoly_BIG)){
  result_BIG[[i]] <- ConvertPolygonToBuffer(ReefPoly_BIG[i,], ID = BigPolygonIDs, i=i); print(i/length(ReefPoly_BIG))}
result_BIG_bind <- do.call("rbind", result_BIG)

reef_polygonsAll <- rbind(result_BIG_bind, result)
reef_polygonsAll_V1 <- aggregate(reef_polygonsAll)
reef_polygonsAll_V1 <- as(reef_polygonsAll_V1, 'SpatialPolygonsDataFrame')
writeOGR(reef_polygonsAll_V1, "/Volumes/Untitled/Raster_files/WCMC_ReefData/BufferedPolygons", "Reef_AsBufferedPolygon", driver="ESRI Shapefile")

# png('figures_extra/TestPoly.png', res = 600, width = 4000, height = 2000)
# plot(reef_polygonsAll_V1, cex = 0.1)
# dev.off()

# Test extraction over raster layer. 
#DepthRaster <- mask(preds_Stack3$Depth, reef_polygonsAll_V1)
#plot(DepthRaster)

# Read in world shapefiles to buffer environmental data so that not extrapolating into pelagic ocean ----

WorldPoly <- readOGR(dsn = '/Volumes/Untitled/Raster_files/TM_WORLD_BORDERS_SIMPL-0')
WorldPoly <- aggregate(WorldPoly)
plot(WorldPoly)
WorldPoly_buffer <- buffer(WorldPoly, width = 1)
WorldPoly_buffer2 <- as(WorldPoly_buffer, 'SpatialPolygonsDataFrame')
#plot(WorldPoly_buffer)
#DepthRaster <- mask(preds_Stack3$Depth, WorldPoly_buffer)
#plot(DepthRaster)
#

# COMBINE REEF AND LAND and RLS buffered site MASKS ----
# Read in RLS data and create buffer around sites. 
Site_Location  <- read.csv(file = 'data_upload/RLS_SiteLocations.csv')
RLS_sites_SPDF <- SpatialPointsDataFrame(cbind(Site_Location$SiteLong, Site_Location$SiteLat), Site_Location, proj = crs('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
RLS_sites_buffer <- buffer(RLS_sites_SPDF, width = 0.5)

# Make full reef mask. 
FullReefMask <- aggregate(rbind(reef_polygonsAll_V1, WorldPoly_buffer2))
FullReefMask2 <- raster::aggregate(rbind(RLS_sites_buffer, FullReefMask, makeUniqueIDs = TRUE))
preds_Stack3_MASK <- mask(preds_Stack3, FullReefMask2)
plot(preds_Stack3_MASK)

# Extract RLS reef and land data from masks to check no missing surveys TESTING. 
#ExtractionTest_AllCells <- extract(preds_Stack3$Depth, RLS_sites_SPDF)
#ExtractionTest_RefinedCells <- extract(DepthRaster, RLS_sites_SPDF)
#sum(is.na(ExtractionTest_AllCells))
#sum(is.na(ExtractionTest_RefinedCells))
#plot(DepthRaster)
#points(RLS_sites_SPDF, col = 'black')
#points(RLS_sites_SPDF[is.na(ExtractionTest_RefinedCells),], col = 'red')


# Run the filling of NAs over all values (this is generally ~50 sites so is a tiny fraction). 
preds_Stack3_MASK$Depth           <- FillSiteNAs(preds_Stack3_MASK$Depth, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$HumanPop        <- FillSiteNAs(preds_Stack3_MASK$HumanPop, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$LandArea        <- FillSiteNAs(preds_Stack3_MASK$LandArea, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$NPP             <- FillSiteNAs(preds_Stack3_MASK$NPP, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$ReefArea        <- FillSiteNAs(preds_Stack3_MASK$ReefArea, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$WaveEnergy      <- FillSiteNAs(preds_Stack3_MASK$WaveEnergy, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$CurrentVelocity <- FillSiteNAs(preds_Stack3_MASK$CurrentVelocity, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$O2              <- FillSiteNAs(preds_Stack3_MASK$O2, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Iron            <- FillSiteNAs(preds_Stack3_MASK$Iron, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Nitrate         <- FillSiteNAs(preds_Stack3_MASK$Nitrate, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$pH              <- FillSiteNAs(preds_Stack3_MASK$pH, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Phosphate       <- FillSiteNAs(preds_Stack3_MASK$Phosphate, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Salinity        <- FillSiteNAs(preds_Stack3_MASK$Salinity, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Silicate        <- FillSiteNAs(preds_Stack3_MASK$Silicate, RLS_sites_SPDF); beepr::beep()
preds_Stack3_MASK$Temperature     <- FillSiteNAs(preds_Stack3_MASK$Temperature, RLS_sites_SPDF); beepr::beep()

# Save this masked raster. 
#writeRaster(preds_Stack3_MASK, 
#            file = '/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack_V3.grd', 
#            options="INTERLEAVE=BAND", 
#            overwrite=TRUE)

# Need to re-mask the depth layer to exclude land - was not performed above properly. 
preds_Stack3_MASK <- stack('/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack_V3.grd')


# Depth raster still contains land. 
plot(preds_Stack3_MASK$Depth)

WorldPoly <- readOGR(dsn = '/Volumes/Untitled/Raster_files/TM_WORLD_BORDERS_SIMPL-0')
WorldPoly <- aggregate(WorldPoly)
plot(WorldPoly)
WorldPoly_buffer <- buffer(WorldPoly, width = -11)
WorldPoly_buffer2 <- as(WorldPoly_buffer, 'SpatialPolygonsDataFrame')

preds_Stack3_MASK$Depth <- mask(preds_Stack3_MASK$Depth, WorldPoly, inverse = T)
preds_Stack3_MASK$Depth[preds_Stack3_MASK$Depth > 100] <- NA

# How many RLS sites are actually NAs? 

# Fill the NAs. 
preds_Stack3_MASK$Depth <- FillSiteNAs(preds_Stack3_MASK$Depth, RLS_sites_SPDF); beepr::beep()

# Need to ensure that no RLS sites with these masks produce NAs. 
writeRaster(preds_Stack3_MASK, 
            file = '/Volumes/Untitled/Raster_files/MarineRasterStack/RasterStack_V4.grd', 
            options="INTERLEAVE=BAND", 
            overwrite=TRUE)
# ----

# FITTING SPECIES' DISTRIBUTION MODELS ----
# Read in MEOW shapefile for setting projection limits (use provinces) ----
ogrInfo(dsn = '/Volumes/Untitled/Raster_files/MEOW')
MEOW <- readOGR(dsn = '/Volumes/Untitled/Raster_files/MEOW', layer = 'PROVINCE')
MEOW_Province <- unionSpatialPolygons(MEOW, MEOW$PROVINCE)
plot(MEOW_Province)

# Function to fit SDM models and extract upper and lower niche limits ----
# Function fits SDMs, and saves RDS file on hard-drive. 
ExtractSDMs <- function(OccurrenceData, EnvironmentalData, MEOW_Province, SaveDirectory = wd()){ 
  
  # Convert lats and longs into spatial points data frame. 
  Range <- SpatialPointsDataFrame(coords = cbind(OccurrenceData$SiteLong, OccurrenceData$SiteLat), 
                                  data = data.frame(OccurrenceData %>% dplyr::select(Occurrence)))
  
  # Assign spatial projection 
  crs(Range) <- '+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'
  
  # Get the province mask. 
  MEOW_Mask <- MEOW_Province[unique(over(Range[Range$Occurrence==1,], MEOW_Province))]
  
  # Crop the environmental layers to provinces 
  print('Cropping environmental data to species range')
  EnvironmentalData_Cropped  <- mask(EnvironmentalData, MEOW_Mask)
  
  # Ensure there are no singluar raster layers
  Check_SD <- cellStats(EnvironmentalData_Cropped, "sd")
  if(sum(Check_SD == 0) != 0){
    #If a column is singluar it breaks the pc. Remove all singluar columns
    EnvironmentalData_Cropped <- dropLayer(EnvironmentalData_Cropped, which(Check_SD == 0))
  }else{NULL}
  
  # Perform PCA on selected raster stack variables using RStoolbox::rasterPCA
  print('Fitting PCA to environmental layers')
  EnvironmentalData_Cropped_PCA <- rasterPCA(dropLayer(EnvironmentalData_Cropped, c(1, 15)), nComp = 5, spca = T, maskCheck = TRUE)
  vars <- summary(EnvironmentalData_Cropped_PCA$model)$sdev^2 
  vars <- vars/sum(vars) 
  EnvironmentalData_Cropped_PCA2 <- dropLayer(EnvironmentalData_Cropped_PCA$map, which(vars < 0.1))
  
  # Stack PCA onto temperature and depth. 
  EnvironmentalData_Cropped_FINAL <- stack(EnvironmentalData_Cropped_PCA2, dropLayer(EnvironmentalData_Cropped, which(!1:15 %in% c(1, 15))))
  
  # Convert to a SDM model object 
  SDM_Data <- sdmData(formula = Occurrence~., train = Range, predictors = EnvironmentalData_Cropped_FINAL)
  
  # Create SDM model
  print('Fitting SDM')
  SDM_Model <- sdm(Occurrence~., data=SDM_Data, 
                   methods=c('brt', 'gam', 'svm', 'rf'), 
                   replication='cv',cv.folds=5) # These are the 4 most flexible models. 
  
  # Aggregate models 
  print('Creating model ensemble')
  dir.create(SaveDirectory)
  SDM_Ensemble <- ensemble(SDM_Model, newdata = EnvironmentalData_Cropped_FINAL, 
                           #filename = 'test.img',
                           filename = paste0(
                             SaveDirectory, 
                             OccurrenceData$SpeciesName[1], '.grd'), 
                           overwrite = T, 
                           mean = T, 
                           setting=list(method='weighted',stat='AUC'))
  

  # Create new raster to aggregate from
  SDM_Ensemble <- sdmvspecies::rescale(SDM_Ensemble)
  T_Lower_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.05, na.rm = T), 3)
  T_Lower_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.05, na.rm = T), 3)
  T_Lower_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.05, na.rm = T), 3)
  T_Upper_0.1  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.1], 0.95, na.rm = T), 3)
  T_Upper_0.25 <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.25], 0.95, na.rm = T), 3)
  T_Upper_0.5  <- round(quantile(EnvironmentalData_Cropped$Temperature[SDM_Ensemble > 0.5], 0.95, na.rm = T), 3)
  
  ModelSummary <- getEvaluation(SDM_Model, stat = c('AUC', 'specificity', 'sensitivity', 'TSS'))
  ModelSummaries <- colMeans(ModelSummary)[2:5]
  ModelSummaries <- t(data.frame(ModelSummaries))
  
  # Return summaries of temperature
  return(as_data_frame(cbind(SpeciesName = as.character(OccurrenceData$SpeciesName[1]),
                             T_Lower_0.1, T_Lower_0.25,T_Lower_0.5, 
                             T_Upper_0.1, T_Upper_0.25, T_Upper_0.5, 
                             round(ModelSummaries,3))))
  
  
}
# Run SDMs for all 704 species (regardless of confidence criteria) and save output models and T_Upper and T_Lower ----

# Read in environmental data
preds_Stack3 <- raster('data_upload/RasterStack_V4.grd')

# Read in and create occupancy data
Site_Location  <- read.csv('data_upload/RLS_SiteLocations.csv')

# Read in RLS actualy data
RLS_20 <- readRDS(file = 'data_derived/RLS_20-For-Qgams-2019-09-28.rds')

# Create occupancy data for all species. 
HighQualityOccupancies <- left_join(RLS_20 %>% dplyr::select(SpeciesName, SiteCode, Presence) %>% dplyr::rename(., Occurrence = Presence), Site_Location)

# Run the SDM function in parallel. Each species takes ~ 20 mins (702 species). 
cl <- makeCluster(4)
registerDoParallel(cl)
sdmModels <- 
  foreach(i=1:length(unique(HighQualityOccupancies$SpeciesName)), 
          .packages=c('tidyr', 'dplyr', 'sdm', 'sp', 'raster', 'rgeos', 'sdmvspecies', 'tryCatchLog', 'RStoolbox', 'filesstrings')) %dopar% 
          {
            tryCatch(ExtractSDMs(OccurrenceData = HighQualityOccupancies %>% filter(SpeciesName == unique(.$SpeciesName)[i]), 
                                 EnvironmentalData = preds_Stack3, 
                                 MEOW_Province = MEOW_Province, 
                                 SaveDirectory = 'data_derived/SDM-Models/'
                                   ), 
                     error = function(e) NA)
          }
stopCluster(cl)

# Save the output created by the SDM model function (T_Upper and T_Lower)
saveRDS(do.call(rbind, sdmModels), file = 'data_derived/sdmModelOutputs.rds')
sdmModels <- readRDS('data_derived/sdmModelOutputs.rds')

# END OF SCRIPT ---------------------------------------------------------------------------------------------------------


