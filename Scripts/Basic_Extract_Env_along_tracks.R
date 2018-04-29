# Script to add environmental variables to GLS tracking data
# from Red-legged Kittiwakes in the Bering Sea and western North Pacific Ocean
# The oceanographic datasets must be in raster form stored locally.
#
# Masters Research
#
# Abram Fleishman
# Updated 27 Oct 2017
rm(list= ls())

# Load Packages -----------------------------------------------------------

library(tidyverse)
library(lubridate)
library(raster)
library(velox)
library(sp)
library(rgeos)
library(ncdf4)
library(stringr)
library(maptools)
library(marmap)
library(sf)

# Sett point buffer in meters
bufferM<-185000

# Set main dir: Sys.info()[6] is the username for the computer.
# fill in the "" with your user name
if(Sys.info()[["user"]]=="abramfleishman") {
  Migration<-"~/Dropbox/RLKI/RLKI_Migrations"
  BS_ENV<-"~/Desktop/Bering_Environmental_Data"
  LandMaskPath<-"~/Desktop/Bering_Environmental_Data/LandMask"
  OutDataPath<-paste0("~/Dropbox/RLKI/SGRK_processedDATA/tracks_11-17_vars_",bufferM/1000,"km_25.rds")
  DataInPath<-paste0(Migration, "/Data/processedLocs/probGLS/RLKI_GLS_allyearsalldata.rda")
}

# load tracks file
tracks<-readRDS(DataInPath)
unique(tracks$band)

# make date col
tracks$date<-as.Date(tracks$datetime)
tracks$week<-week(tracks$datetime)

# Load the resampled vars (all 0.25 degre resolution)
ssh<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/ssh_0.25.rds"))
uo<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/uo_0.25.rds"))
vo<-readRDS(paste0(BS_ENV,"/global-analysis-forecast-phy-001-024/vo_0.25.rds"))

# Create a date table for the global-analysis-forecast-phy vars
# This table will look up the name of the layer in the raster stack and connect
# it with a real world date

DateTable<-NULL

for(i in 1:length(ssh@layers)) DateTable[i]<-ssh[[i]]@data@names

DateTable<-data.frame(names=DateTable,  DateTime=ymd_hms(gsub("X","",DateTable))) %>%
  mutate(
    year=year(DateTime),month=month(DateTime),day=day(DateTime),
    hour=hour(DateTime),minute=minute(DateTime),second=second(DateTime),
    Date=as.Date(DateTime)
  )

# Loop for global-analysis-forecast-phy vars ---------------------------------


# Create empty columns to populate for all the vars you are interested in
tracks$eke<-NA
tracks$eke_g<-NA

tracks$uo<-NA
tracks$vo<-NA

tracks$ssh<-NA
tracks$ssh_g<-NA

# Get the unique dates in the tracks data
dates<-sort(unique(tracks$date))

# go through each date in the tracks data, extract the vars at each of the
# points with a buffer, aggragate the buffered data to end up with a single
# value for each point.
#
# Currently using the scale of error 185km to buffer and taking the mean.
bufferM<-185000

# install.packages("velox")

for(i in 1:length(dates)){
  print(dates[i])

  # select the date you want
  datex<-as.character(DateTable$names[DateTable$Date==dates[i]])[1]

  # skip the date if it is not in the raster data
  if(length(datex)<1|is.na(datex)) next

  # make a spDataFrame for that date only
  pt<-data.frame(lon360=tracks$lon360[which(tracks$date==dates[i] & is.na(tracks$ssh))],
                 lat=tracks$lat[which(tracks$date==dates[i] & is.na(tracks$ssh))])

  coordinates(pt)<-cbind(tracks$lon360[which(tracks$date==dates[i]&is.na(tracks$ssh))],
                         tracks$lat[which(tracks$date==dates[i]&is.na(tracks$ssh))])

  crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180")

  # Reproject to an equal area projection so that you can buffer
  # This projection was for the Bering Sea
  pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                            +ellps=WGS84 +units=m +no_defs"))

  # Buffer the points by the specified raduis
  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  # Make the buffered points into a Sptail Polygon Data Frame
  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)

  # Reproject to the coordinates of your dataset
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))

  # Createa velox raster dataset for the date of the layer
  # Velox is very fast (20x) for extraction
  vx<- velox(ssh[[datex]])

  # extract out the mean env var at every point from the date of interest
  # using function(x)mean(x,na.rm=T) allows to ignore NA data from clouds or coasts
  # Here you could change for median, sd, or max etc... depending on your interest
  tracks$ssh[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  # This is calculating the slope of the layer see terrain for details.
  # terrain is designed for calculating topographic slope but I have been using
  # it to calculate gradients in SST and SSH as well
  vx<- velox(terrain(ssh[[datex]],opt = "slope",neighbors = 8))
  tracks$ssh_g[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  # eke
  vx<-velox(1/2*(vo[[datex]]^2+uo[[datex]]^2))
  tracks$eke[which(tracks$date==dates[i])]<-vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  vx<-velox(terrain((1/2*(vo[[datex]]^2+uo[[datex]]^2)),opt = "slope",neighbors = 8))
  tracks$eke_g[which(tracks$date==dates[i])]<-vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  # uo
  vx<- velox(uo[[datex]])
  tracks$uo[which(tracks$date==dates[i])]<-vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  # vo
  vx<- velox(vo[[datex]])
  tracks$vo[which(tracks$date==dates[i])]<-vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

}

saveRDS(tracks,OutDataPath)


# Bathymetry     --------------------------------------------------------------
#
# Download NOAA ETOPO1 data
# #resolution 1.85km*27=49.9km2
# res 15 = 15 minutes = .25 degrees
# res 1 = 1 arcminute = 0.0166 degrees
aleu <- getNOAA.bathy(130, -140, 30, 75, resolution = 15,
                      antimeridian = T,keep = T)
# Make it a raster
bathy<-as.raster(aleu)

# plot it
plot(bathy)

# set values >0 to NA
bathy[bathy>=0]<-NA

# Make it a velox
bathyvx<- velox(bathy)

# calculate Aspect, Slope and TPI (topographic Position Index)
aspectvx<- velox(terrain(bathy,opt="aspect",neighbors = 8))
slopevx<- velox(terrain(bathy,opt="slope",neighbors = 8))
TPIvx<- velox(terrain(bathy,opt="TPI",neighbors = 8))


# Add the empty columns to populate
tracks$Depth<-NA
tracks$Aspect<-NA
tracks$Slope<-NA
tracks$BPI<-NA


# loop through the day of the year and extract for the tracking data
# This isa loop because I had memory issues when trying it all at once
for(i in unique(tracks$doy)){

  print(paste("DOY:",i))

  # make a spatial points DF
  tracks_sp <- SpatialPointsDataFrame(coords=cbind(tracks$lon360[tracks$doy==i&is.na(tracks$Depth)],
                                                   tracks$lat[tracks$doy==i&is.na(tracks$Depth)]),
                                      data=data.frame(tracks$band[tracks$doy==i&is.na(tracks$Depth)]),
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_wrap=180"))

  # Reproject in to an equal area projection to extract with the buffer
  pt_sp<-spTransform(tracks_sp,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                                   +ellps=WGS84 +units=m +no_defs"))

  # Buffer
  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  # Spatial polygons
  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)

  # Reproject back to your native env var projects
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))

  # Attach Depth ------------------------------------------------------------
  tracks$Depth[tracks$doy==i] <- bathyvx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$Aspect[tracks$doy==i] <- aspectvx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$Slope[tracks$doy==i] <- slopevx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  tracks$BPI[tracks$doy==i] <- TPIvx$extract(spdf, fun=function(x)mean(x,na.rm=T))
}

saveRDS(tracks,OutDataPath)

# Distance to Shore in m-------------------------------------------------------

# use highest resolution
aleu <- getNOAA.bathy(130, -140, 30, 75, resolution = 1,
                      antimeridian = T,keep = T)
# Make it a raster
bathy<-as.raster(aleu)

# plot it
plot(bathy)

# set values >0 to NA
bathy[bathy>=0]<-NA

# Extract the shoreline
Shore<-rasterToContour(bathy,level=-1)

# reproject
Shore<-spTransform(Shore,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

# calculate distances to ever shore segment
Shoredist<-  gDistance(tracks_sp_laea,Shore,byid = T)

# find minimum distance to shore
tracks$dist2shore<-apply(Shoredist, 2, function(x) min(x, na.rm = TRUE))

saveRDS(tracks,OutDataPath)

# Inside bering sea? ------------------------------------------------------
# use Bering_sea.shp to calculate if a point is in the Bering sea
# http://www.marineregions.org/gazetteer.php?p=details&id=4310

# Make Spatial points DF
tracksTemp<-data.frame(tracks)
coordinates(tracksTemp)<-cbind(tracksTemp$lon,tracksTemp$lat)
proj4string(tracksTemp)<-CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")

# read in shape file
ber<-read_sf(paste0(BS_ENV,"/iho"),layer="iho",crs="+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0") %>% as("Spatial")
proj4string(ber)
plot(ber)

# use over function to extract overlap points
InBer<-over(tracksTemp,ber)

# add a var to the tracks DF
tracks$bs<-as.character(InBer$name)



# load packages for downloading environmental vars ------------------------
devtools::install_github('ropensci/rerddap')
devtools::install_github('ropensci/plotdap')
devtools::install_github("rmendels/rerddapXtracto")
# library(xtractomatic)

library(rerddapXtracto)
# version()
library(rerddap)

# esrlNcepRe NCEP Reanalysis 1 --------------------------------------------

SST<-info( "erdMBsstd8day")
xtractomatic::searchData(searchList = "varname:sst")

# ypos<-c(16,26)
# xpos<-c(-162,-150)+360
# tpos<-c("2017-04-24","2017-05-28")

# Aqua on MODIS - CHL - monthly composites "2010-06-01","2017-04-15"
erdMBsstd8day<-rxtracto_3D(dataInfo = SST,parameter = "sst",
                           xcoord = c(-162,-150)+360,
                           ycoord = c(16,26),
                           zcoord = c(0.0,0.0),
                           zName="altitude",
                           tName = "time",
                           tcoord = c("2017-04-24T00:00:00Z","2017-05-28T00:00:00Z"),
                           urlbase = "http://coastwatch.pfeg.noaa.gov/erddap/",verbose = T,)
str(erdMBsstd8day)

library(raster)
sst<-stack(x = 'sst.nc')
names(sst)<-paste0("X",seq(as.Date("2017-04-24"),as.Date("2017-05-28"),1))

tracks<-read_csv("~/Desktop/Olivias_data/")
tracks<-data.frame(date=seq(as.Date("2017-04-24"),as.Date("2017-05-28"),1),lat=18,lon=-155+360)
dates<-unique(tracks$date)
# Loop to add sst and sst_g
tracks$sst<-NA
tracks$sst_g<-NA
for(i in 1:length(dates)){
  print(dates[i])
  datex<-paste0("X",year(dates[i]),".",
                str_pad(string = month(dates[i]),width = 2,pad = "0",side = "left"),".",
                str_pad(string = day(dates[i]),width = 2,pad = "0",side = "left"))

  pt<-data.frame(lon=tracks$lon[which(tracks$date==dates[i]&is.na(tracks$sst))],
                 lat=tracks$lat[which(tracks$date==dates[i]&is.na(tracks$sst))])

  coordinates(pt)<-cbind(tracks$lon[which(tracks$date==dates[i]&is.na(tracks$sst))],
                         tracks$lat[which(tracks$date==dates[i]&is.na(tracks$sst))])
  crs(pt)<-CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")


  pt_sp<-spTransform(pt,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

  spol <- gBuffer(pt_sp, width=bufferM, byid=TRUE)

  spdf <- SpatialPolygonsDataFrame(spol, data.frame(id=1:length(spol)), FALSE)
  spdf <- spTransform(spdf,CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +lon_wrap=180"))

  # sst
  vx<- velox(sst[[datex]])
  tracks$sst[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))

  vx<- velox(terrain(sst[[datex]],opt = "slope",neighbors = 8))
  tracks$sst_g[which(tracks$date==dates[i])] <- vx$extract(spdf, fun=function(x)mean(x,na.rm=T))
}
