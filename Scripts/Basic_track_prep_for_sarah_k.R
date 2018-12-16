# Script to calculate:
# - distance from each point along a track to a fixed point (colony)
# - distance from each point  to a fixed feature like a shoreline,
#   bathymetric feature (200m isobath), or any other line (ship track)
# - distance between points on a track
# - time between points on a track


# For Sarah K. @ UCSC
# Abram Fleishman
# Updated 15 Dec 2018


# Install trakR packages from github --------------------------------------


# install.packages("devtools") # for installing packages from github
# devtools::install_github("abfleishman/trakR") # install my package

library(trakR)

# Load tracks -------------------------------------------------------------


# This should be one file for all the tracks or create a data.frame with all the tracks
# tracks<-read.csv("path/to/tracks.csv")
tracks<-readRDS("~/Dropbox (Personal)/RLKI/SGRK_processedDATA/GPS_RLKI_rawGPS_compiled_24Jun18.rda")
head(tracks)

# Add interpoint time interval
tracks$InterpointTime<-InterpointTime(tracks = tracks,
                                      ID = "CaptureID",
                                      DateTime = "DateTime")

# Add interpoint distance interval uses argosfilter::distanceTrack
# install.packages("argosfilter")
tracks$InterpointDist<-InterpointDist(tracks = tracks,
                                      ID = "CaptureID",
                                      lat = "Latitude",
                                      lon = "Longitude")

# Calculate distance to colony
tracks$Dist2Col<-Dist2Colony(tracks = tracks,
                             dataLat = "Latitude",
                             dataLon = "Longitude",
                             ColonyLat = 56.6,
                             ColonyLong = -169.5)


# Distance to Shore in m-------------------------------------------------------
#
# This can be used to fine the distance to any SpatialLinesDataFrame
#
# install.packages("marmap") #for bathymetry
library(marmap)
library(raster)
# Download NOAA ETOPO1 data
# #resolution 1.85km*27=49.9km2
# res 15 = 15 minutes = .25 degrees
# res 1 = 1 arcminute = 0.0166 degrees
# use highest resolution
aleu <- marmap::getNOAA.bathy(130, -140, 30, 75, resolution = 15,
                      antimeridian = T, keep = F)
# Make it a raster
bathy<-marmap::as.raster(aleu)

# plot it
raster::plot(bathy)

# Extract the shoreline (could do 200m isobath using level=-200)
Shore<-rasterToContour(bathy,level=-1)

# reproject into a equal area projection (might should do equal distance...)
Shore_laea<-spTransform(Shore,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
                             +ellps=WGS84 +units=m +no_defs"))

# make a spatial points DF from your tracks
tracks_sp <- SpatialPointsDataFrame(coords=cbind(tracks$Longitude,
                                                 tracks$Latitude),
                                    data=data.frame(tracks$CaptureID),
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_wrap=180"))

# Reproject in to an equal area projection
pt_sp<-spTransform(tracks_sp,CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500000 +y_0=1000000
+ellps=WGS84 +units=m +no_defs"))

# calculate distances to ever shore segment (this is slow with lots of points or a complex line file)
Shoredist<-  rgeos::gDistance(pt_sp,Shore_laea,byid = T)

# find minimum distance to shore
tracks$dist2shore<-apply(Shoredist, 2, function(x) min(x, na.rm = TRUE))
