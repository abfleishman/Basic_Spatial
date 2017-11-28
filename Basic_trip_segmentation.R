# Trip segmentation for animal tracking data

# Abram Fleishman
# 28 Nov 2017
# Prepared for Olivia T. @ SJSU

library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)

# load all the functions
for(i in list.files("Functions/",full.names = T)) source(i)

# Load tracks -------------------------------------------------------------
# This should be one file for all the tracks or create a data.frame with all the tracks
tracks<-read.csv("path/to/tracks.csv")
tracks<-readRDS("~/Dropbox/SGRK_processedDATA/GPSrlki_raw.rda")
head(tracks)

# Add interpoint time interval
tracks$InterpointTime<-InterpointTime(tracks = tracks,ID = "CaptureID",DateTime = "DateTime")

# Add interpoint distance interval
install.packages("argosfilter")
tracks$InterpointDist<-InterpointDist(tracks = tracks,ID = "CaptureID",lat = "Latitude",lon = "Longitude")

# Calculate distance to colony
tracks$Dist2Col<-Dist2Colony(tracks = tracks,ColonyLat = 56.6,ColonyLong = -169.5)

# segment into tripos
tracks_w_trips<-MakeTrip(tracks = tracks,ID = "CaptureID",DistCutOff = 100,Dist2Colony = "Dist2Col")
head(tracks_w_trips)

# Plot a bird to check
ggplot(tracks_w_trips[tracks_w_trips$CaptureID=="B27",],aes(Latitude,Longitude,col=factor(TripNum)))+
  geom_path(size=.2)+
  geom_point(data=tracks_w_trips[tracks_w_trips$CaptureID=="B27"&tracks_w_trips$ColonyMovement%in%c("Out","In"),])
