#Download data from https://climate-scenarios.canada.ca/?page=cmip6-scenarios#gridded_data

setwd("/Volumes/OneTouch5GB/ClimateModels/")
dropbox_dir<-file.path("/Users/ruthrivkin/Library/CloudStorage/Dropbox/Postdoc_2021-2024/Polar_Bears/WholeGenomes/") #set directory for saving in dropbox 

# loads the packages used in this guide
library(ncdf4)
library(raster)
library(ggplot2)

##HISTORICAL DATA
#Load monthly historical data
temp <- brick("Temp (K)/tas_Global_ensemble_historical_r1i1p1f1_mean.nc") #1980 layers
wind <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_historical_r1i1p1f1_mean.nc") #1980 layers
snow <- brick("Snow depth (m)/snd_Global_ensemble_historical_r1i1p1f1_mean.nc") #1980 layers
ice.conc <- brick("Sea ice conc (%)/siconc_Global_ensemble_historical_r1i1p1f1_mean.nc") #1980 layers
ice.thick <- brick("Sea ice thickness (m)/sithick_Global_ensemble_historical_r1i1p1f1_mean.nc") #1980 layers

# plots the variable
plot(temp)

# Find the time points that make sense to average (2000-2014)
plot(temp[[1979:1980]]) #Last time point is Dec 2014
plot(temp[[1800:1801]]) #1801 corresponds to Jan 2000

meantemp <- mean(temp[[1801:1980]])
plot(meantemp)
#Convert mean temp to celsius
meantemp <- meantemp - 273.15

#Get mean values
meanwind <- mean(wind[[1801:1980]])
meansnow <- mean(snow[[1801:1980]])
meaniceconc <- mean(ice.conc[[1801:1980]])
meanicethick <- mean(ice.thick[[1801:1980]])

#Stack the layers
historical <- brick(meantemp, meanwind, meansnow, meaniceconc, meanicethick)
names(historical) <- c("AirTemp","WindSpeed", "SnowDepth", "IceConc", "IceThick")

#convert longitude from 0-360 to -180-180, already done for the StableCLim dataset
historical180 <- raster::rotate(historical)
plot(historical180)


#save as new raster
writeRaster(historical180, filename=file.path("HistoricalEnvironment.tif"), format="GTiff", overwrite=TRUE)



##Projections

#Temp
temp_ssp1 <- brick("Temp (K)/tas_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
temp_ssp2 <- brick("Temp (K)/tas_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
temp_ssp3 <- brick("Temp (K)/tas_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
temp_ssp4 <- brick("Temp (K)/tas_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
temp_ssp5 <- brick("Temp (K)/tas_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

plot(temp_ssp1[[1019:1020]]) #1020 corresponds to Dec2099
plot(temp_ssp1[[900:901]]) #901 corresponds to Jan 2090

temp_ssp1 <- mean(temp_ssp1[[901:1020]])
temp_ssp2 <- mean(temp_ssp2[[901:1020]])
temp_ssp3 <- mean(temp_ssp3[[901:1020]])
temp_ssp4 <- mean(temp_ssp4[[901:1020]])
temp_ssp5 <- mean(temp_ssp5[[901:1020]])


temp_ssp1 <- temp_ssp1 - 273.15
temp_ssp2 <- temp_ssp2 - 273.15
temp_ssp3 <- temp_ssp3 - 273.15
temp_ssp4 <- temp_ssp4 - 273.15
temp_ssp5 <- temp_ssp5 - 273.15

temp_ssp <- brick(temp_ssp1, temp_ssp2, temp_ssp3, temp_ssp4, )

names(temp_ssp) <- c("temp_ssp1","temp_ssp2", "temp_ssp3", "temp_ssp4", "temp_ssp5")
temp_ssp180 <- raster::rotate(temp_ssp)

writeRaster(temp_ssp180, filename=file.path("Rasters/SSP_AirTemp_2100.tif"), format="GTiff", overwrite=TRUE)

#Temp
temp_ssp1 <- brick("Temp (K)/tas_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
temp_ssp2 <- brick("Temp (K)/tas_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
temp_ssp3 <- brick("Temp (K)/tas_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
temp_ssp4 <- brick("Temp (K)/tas_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
temp_ssp5 <- brick("Temp (K)/tas_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

plot(temp_ssp1[[1019:1020]]) #1020 corresponds to Dec2099
plot(temp_ssp1[[900:901]]) #901 corresponds to Jan 2090

temp_ssp1 <- mean(temp_ssp1[[901:1020]])
temp_ssp2 <- mean(temp_ssp2[[901:1020]])
temp_ssp3 <- mean(temp_ssp3[[901:1020]])
temp_ssp4 <- mean(temp_ssp4[[901:1020]])
temp_ssp5 <- mean(temp_ssp5[[901:1020]])


temp_ssp1 <- temp_ssp1 - 273.15
temp_ssp2 <- temp_ssp2 - 273.15
temp_ssp3 <- temp_ssp3 - 273.15
temp_ssp4 <- temp_ssp4 - 273.15
temp_ssp5 <- temp_ssp5 - 273.15

temp_ssp <- brick(temp_ssp1, temp_ssp2, temp_ssp3, temp_ssp4, temp_ssp5)

names(temp_ssp) <- c("temp_ssp1","temp_ssp2", "temp_ssp3", "temp_ssp4", "temp_ssp5")
temp_ssp180 <- raster::rotate(temp_ssp)

writeRaster(temp_ssp180, filename=file.path("Rasters/SSP_AirTemp_2100.tif"), format="GTiff", overwrite=TRUE)


#Wind
wind_ssp1 <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
wind_ssp2 <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
wind_ssp3 <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
wind_ssp4 <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
wind_ssp5 <- brick("Near-surface wind speed (m:s)/sfcWind_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

wind_ssp1 <- mean(wind_ssp1[[901:1020]])
wind_ssp2 <- mean(wind_ssp2[[901:1020]])
wind_ssp3 <- mean(wind_ssp3[[901:1020]])
wind_ssp4 <- mean(wind_ssp4[[901:1020]])
wind_ssp5 <- mean(wind_ssp5[[901:1020]])

wind_ssp <- brick(wind_ssp1, wind_ssp2, wind_ssp3, wind_ssp4, wind_ssp5)

names(wind_ssp) <- c("wind_ssp1","wind_ssp2", "wind_ssp3", "wind_ssp4", "wind_ssp5")
wind_ssp180 <- raster::rotate(wind_ssp)

writeRaster(wind_ssp180, filename=file.path("Rasters/SSP_WindSpeed_2100.tif"), format="GTiff", overwrite=TRUE)


#Snow Depth
snow_ssp1 <- brick("Snow depth (m)/snd_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
snow_ssp2 <- brick("Snow depth (m)/snd_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
snow_ssp3 <- brick("Snow depth (m)/snd_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
snow_ssp4 <- brick("Snow depth (m)/snd_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
snow_ssp5 <- brick("Snow depth (m)/snd_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

snow_ssp1 <- mean(snow_ssp1[[901:1020]])
snow_ssp2 <- mean(snow_ssp2[[901:1020]])
snow_ssp3 <- mean(snow_ssp3[[901:1020]])
snow_ssp4 <- mean(snow_ssp4[[901:1020]])
snow_ssp5 <- mean(snow_ssp5[[901:1020]])

snow_ssp <- brick(snow_ssp1, snow_ssp2, snow_ssp3, snow_ssp4, snow_ssp5)

names(snow_ssp) <- c("snow_ssp1","snow_ssp2", "snow_ssp3", "snow_ssp4", "snow_ssp5")
snow_ssp180 <- raster::rotate(snow_ssp)

writeRaster(snow_ssp180, filename=file.path("Rasters/SSP_SnowDepth_2100.tif"), format="GTiff", overwrite=TRUE)


#Sea ice concentration
ice.conc_ssp1 <- brick("Sea ice conc (%)/siconc_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
ice.conc_ssp2 <- brick("Sea ice conc (%)/siconc_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
ice.conc_ssp3 <- brick("Sea ice conc (%)/siconc_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
ice.conc_ssp4 <- brick("Sea ice conc (%)/siconc_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
ice.conc_ssp5 <- brick("Sea ice conc (%)/siconc_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

ice.conc_ssp1 <- mean(ice.conc_ssp1[[901:1020]])
ice.conc_ssp2 <- mean(ice.conc_ssp2[[901:1020]])
ice.conc_ssp3 <- mean(ice.conc_ssp3[[901:1020]])
ice.conc_ssp4 <- mean(ice.conc_ssp4[[901:1020]])
ice.conc_ssp5 <- mean(ice.conc_ssp5[[901:1020]])

ice.conc_ssp <- brick(ice.conc_ssp1, ice.conc_ssp2, ice.conc_ssp3, ice.conc_ssp4, ice.conc_ssp5)

names(ice.conc_ssp) <- c("ice.conc_ssp1","ice.conc_ssp2", "ice.conc_ssp3", "ice.conc_ssp4", "ice.conc_ssp5")
ice.conc_ssp180 <- raster::rotate(ice.conc_ssp)

writeRaster(ice.conc_ssp180, filename=file.path("Rasters/SSP_IceConc_2100.tif"), format="GTiff", overwrite=TRUE)


#Sea ice thickness
ice.thick_ssp1 <- brick("Sea ice thickness (m)/sithick_Global_ensemble_ssp119_r1i1p1f1_mean.nc") #1020 layers
ice.thick_ssp2 <- brick("Sea ice thickness (m)/sithick_Global_ensemble_ssp245_r1i1p1f1_mean.nc") 
ice.thick_ssp3 <- brick("Sea ice thickness (m)/sithick_Global_ensemble_ssp370_r1i1p1f1_mean.nc") 
ice.thick_ssp4 <- brick("Sea ice thickness (m)/sithick_Global_ensemble_ssp460_r1i1p1f1_mean.nc") 
ice.thick_ssp5 <- brick("Sea ice thickness (m)/sithick_Global_ensemble_ssp585_r1i1p1f1_mean.nc") #

ice.thick_ssp1 <- mean(ice.thick_ssp1[[901:1020]])
ice.thick_ssp2 <- mean(ice.thick_ssp2[[901:1020]])
ice.thick_ssp3 <- mean(ice.thick_ssp3[[901:1020]])
ice.thick_ssp4 <- mean(ice.thick_ssp4[[901:1020]])
ice.thick_ssp5 <- mean(ice.thick_ssp5[[901:1020]])

ice.thick_ssp <- brick(ice.thick_ssp1, ice.thick_ssp2, ice.thick_ssp3, ice.thick_ssp4, ice.thick_ssp5)

names(ice.thick_ssp) <- c("ice.thick_ssp1","ice.thick_ssp2", "ice.thick_ssp3", "ice.thick_ssp4", "ice.thick_ssp5")
ice.thick_ssp180 <- raster::rotate(ice.thick_ssp)

writeRaster(ice.thick_ssp180, filename=file.path("Rasters/SSP_IceThickness_2100.tif"), format="GTiff", overwrite=TRUE)


#Convert to stack

rastcurrent <- raster::stack(meantemp)
rast2050 <- raster::stack(rcp85.2050)
rast2100 <- raster::stack(rcp85.2100)

# Initial crop- probably not necessary but it works for this pipeline. can make more efficient later
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

#crop rasters
crop.temp <- crop(rastcurrent,crop.extent) 
crop.2050 <- crop(rast2050, crop.extent)
crop.2100 <- crop(rast2100, crop.extent)




#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../Species_Dist_Modelling/Gradient_Forest/Environment/chip12_ind_coord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(crop.temp)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
temp.data = tibble(ID = 1:nrow(pts@coords),
                   Lon = pts$x,
                   Lat = pts$y
)
temp.data


#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
#Create raster list with 50 km buffer around each sample
temp = c(crop.temp, crop.2050, crop.2100)

rasters <- raster::stack(temp)
nlayers(rasters)

#Create raster list with 50 km buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

names(store_data) = names(rasters)
temp.data = bind_cols(temp.data, as_tibble(store_data))
temp.data

# Check each column for NA values
na.check = map_int(temp.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
temp.data.nona = temp.data %>% drop_na


#Export data to a csv file and individual id coordinates
sites$ID <- c(1:1476)
head(sites)
temp_withcoord <- merge(sites, temp.data.nona)
head(temp_withcoord)
str(temp_withcoord)

# Rename columns and convert to celsius, not necessary for StableCLim data
library(dplyr)

temp_mean <- temp_withcoord %>% 
  # mutate(MeanTemp = layer.1 - 273.15, Temp2050 = layer.2 - 273.15, Temp2100 = layer.3 - 273.15) %>% 
  dplyr::rename( MeanTemp = layer.1 , Temp2050 = layer.2, Temp2100 = layer.3)
head(temp_mean)
str(temp_mean)

#write to csv
write_csv(temp_mean, "23.07.20_temp_withallcoord.csv")


#Download other RCP scenarios
#Download StableClim Dataset (https://adelaide.figshare.com/articles/dataset/StableClim/12197976?file=24142217)
tsannual_26 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp26_ts.nc") #251 layers corresponding to 1850-2100
tsannual_45 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp45_ts.nc") #251 layers corresponding to 1850-2100
tsannual_60 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp60_ts.nc") #251 layers corresponding to 1850-2100
tsannual_85 <- brick("/Volumes/OneTouch/PostDoc_2021-2025/Polar Environmental Data/StableClim_V1.0.1/ncdf/annual/StableClim_AnnMean_rcp85_ts.nc") #251 layers corresponding to 1850-2100


# plot time points
plot(tsannual[[240:251]])

#Create new rasterbrick for the same time period as bio-oracle (2000-2014) and future environments (2040-2050 annual mean and 2090-2100 annual means)

#first can we rename layers so they make sense
names(tsannual_26)
names(tsannual_26) <- c(1850:2100)
names(tsannual_26) #okay that worked well

names(tsannual_45) <- c(1850:2100)
names(tsannual_60) <- c(1850:2100)
names(tsannual_85) <- c(1850:2100)


#Get mean value for 2090-2100

rcp26 <- tsannual_26[[241:251]]
plot(rcp26)
rcp26 <- mean(rcp26)

rcp45 <- tsannual_45[[241:251]]
rcp45 <- mean(rcp45)

rcp60 <- mean(tsannual_60[[241:251]])

rcp85 <- mean(tsannual_85[[241:251]])

#convert longitude from 0-360 to -180-180, already done for the StableCLim dataset
#raster_mean180 <- rotate(meantemp)

#save as new raster
writeRaster(rcp26, filename=file.path("Rasters/rcp26_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp45, filename=file.path("Rasters/rcp45_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp60, filename=file.path("Rasters/rcp60_StableClim.tif"), format="GTiff", overwrite=TRUE)
writeRaster(rcp85, filename=file.path("Rasters/rcp85_StableClim.tif"), format="GTiff", overwrite=TRUE)

#Convert to stack

rastcurrent <- raster::stack(meantemp)
rast2050 <- raster::stack(rcp85.2050)
rast2100 <- raster::stack(rcp85.2100)

# Initial crop- probably not necessary but it works for this pipeline. can make more efficient later
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

#crop rasters
crop.rcp26 <- crop(rcp26,crop.extent) 
crop.rcp45 <- crop(rcp45, crop.extent)
crop.rcp60 <- crop(rcp60, crop.extent)
crop.rcp85 <- crop(rcp85, crop.extent)




#4. Extract data from raster for current conditions
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../GFAllEnvironments/23.11.01_tempandice_withcoord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(crop.temp)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
temp.data = tibble(ID = 1:nrow(pts@coords),
                   Lon = pts$x,
                   Lat = pts$y
)
temp.data


#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
#Create raster list with 50 km buffer around each sample
temp = c(crop.rcp26, crop.rcp45, crop.rcp60, crop.rcp85)
names(temp) = c("rcp26", "rcp45", "rcp60", "rcp85")
rasters <- raster::stack(temp)
nlayers(rasters)

#Create raster list with 50 km buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

names(store_data) = names(rasters)
temp.data = bind_cols(temp.data, as_tibble(store_data))
temp.data

# Check each column for NA values
na.check = map_int(temp.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
temp.data.nona = temp.data %>% drop_na


#Export data to a csv file and individual id coordinates
sites$ID <- c(1:411)
head(sites)
temp_withcoord <- merge(sites, temp.data.nona)
head(temp_withcoord)
str(temp_withcoord)


#write to csv
write_csv(temp_withcoord, "../GFAllEnvironments/24.04.15_temp_withallcoord.csv")


#Do it again for  ice thickness (pull other scenarios)

ice.future = c("BO22_RCP26_2100_icethickmean_ss","BO22_RCP45_2100_icethickmean_ss","BO22_RCP60_2100_icethickmean_ss","BO22_RCP85_2100_icethickmean_ss"
)

# Combine present-day and future vectors

# Download rasters  and import into R
ice.raster = load_layers(ice.future)

# Crop rasters to boundary extent
crop.extent <- extent(-145,-45, 50, 90) #define boundary/extent

ice.raster.crop = crop(ice.raster, crop.extent)
names(ice.raster.crop) <- c("ice.rcp26", "ice.rcp45", "ice.rcp60", "ice.rcp85")


#4. Extract data from raster
# Load site locations (needs coordinates for each sampling/target site)
sites <- read.csv("../GFAllEnvironments/24.04.15_temp_withallcoord.csv")

# set x and y for lat and long
x <- sites$Longitude
y <- sites$Latitude

# set up coordinates as a tibble
pts <- cbind(x, y) %>% as_tibble()
pts <- SpatialPoints(pts, proj4string = CRS("+proj=longlat +datum=WGS84"))
pts

#Check points CRS matches raster CRS.
projection(pts) == projection(ice.raster)

#Create a tibble or data.frame to store Bio-ORACLE marine data for each point.
ice.data = tibble(ID = 1:nrow(pts@coords),
                  Lon = pts$x,
                  Lat = pts$y
)
ice.data

rasters <- raster::stack(ice.raster.crop)
nlayers(rasters)

#Create raster list with 20 km buffer around each sample
store_data = list()
for (i in 1:nlayers(rasters)){
  store_data[[i]] = raster::extract(rasters[[i]], pts, buffer=2500
                                    , fun=mean)
}

#Add the extracted data as new columns to ice.data.
# Name variables in the list and then combine data
names(store_data) = names(rasters)
ice.data = bind_cols(ice.data, as_tibble(store_data))
ice.data
# Check each column for NA values
na.check = map_int(ice.data, ~sum(is.na(.)))
summary(na.check > 0)

# Remove NA records
ice.data.nona = ice.data %>% drop_na



#Export data to a csv file and individual id coordinates
sites$ID <- c(1:411)
head(sites)
ice_thickness_withcoord <- merge(sites, ice.data.nona)
head(ice_thickness_withcoord)
write_csv(ice_thickness_withcoord, "../GFAllEnvironments/24.04.15_temp_withallcoord.csv")



#Make some plots
#3. Plot rasters
# Define colour scheme
cols = colorRampPalette(c("#5E85B8","#EDF0C0","#C13127"))

# Define a boundary
boundary = extent(c(xmin = -180, xmax = 20, ymin = 50, ymax = 90))


#load rasters
sites <- read.csv("../GFAllEnvironments/23.11.01_tempandice_withcoord.csv")
coord <- dplyr::select(sites, c("Longitude", "Latitude"))

ice.thick <- c("BO22_icethickmean_ss")
ice.thick <- load_layers(ice.thick, datadir="Rasters/")

ice.cov <- c("BO22_icecovermean_ss")
ice.cov <- load_layers(ice.cov, datadir="Rasters/")

temp <- raster("../Environment/MeanTemp_Current.tif")


# Crop rasters to boundary extent
it.crop = crop(ice.thick, boundary)
ic.crop = crop(ice.cov, boundary)
t.crop = crop(temp, boundary)
t.crop = t.crop -273.15
# Plot ice variables
par(mar = c(5.1, 4.1, 4.1, 5))

pdf("Ice Thickness.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(it.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
     xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()

pdf("Ice Cover.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(ic.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
       xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()

pdf("Temperature.pdf",  width = 6.78, height = 5.3)
plot(NULL, xlim = c(-135, -50), ylim = c(52, 85),
     xlab = "Longitude", ylab = "Latitude", main = "Ice Thickness")
plot(t.crop, col=cols(200), add = TRUE)
points(coord, pch = 19, cex = .75, 
       xlab = "Longitude (°W)", ylab = "Latitude (°N)")
maps::map(add = T, interior = F, fill = T,col = "lightgrey" )
dev.off()
