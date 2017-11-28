# Extract the "most probable track" from the probGLS output and plot them together
#
# This script extracts the "most probable track" from the probGLS output for
# each bird and compiles them into a single "flat" data file.  One row
# corrisponds to a single estimated location.

library(readr) #for read_csv that inturprets dates for you
library(lubridate)
library(dplyr)
library(tidyr)
library(stringr)

library(mapdata)
library(ggplot2)


# List all the twl_0.5_csv
# Files<-list.files(file.path(Migration,"Data/processedTwl"),pattern = "twl_0.5.csv",full.names = T,recursive = T)
Files<-list.files("~/Dropbox/RLKI/RLKI_Migrations/Data/processedTwl",pattern = "twl_0.5.csv",full.names = T,recursive = T)
Files<-Files[str_detect(Files,"15-16")]

#  loop to extract most probable track from pr data
outData<-NULL

for(i in 1:length(Files)){
  #if(i==17) next
  twlPath<-Files[i]

  # Extract season which is used to tell what kind of tag
  Season<-basename(dirname(twlPath))

  pr<-readRDS(file = file.path("~/Desktop/RKLI_FinalGLSTracks_2010-16/",Season,gsub("_twl_0.5.csv","_0.5_noSST_pGLS.rda",basename(twlPath))))

  # Seperate metadata saved in the file name
  Meta<-separate(data.frame(File=gsub("_twl_0.5.csv","_0.5_pGLS",basename(twlPath))),File,c("Tag","downDate","downTime","BirdID","Band","twlThresh","Fun"),sep="_")

  Meta$Band<-gsub("driftadj","",Meta$Band)

  print(Meta)
  # Make a data frame with the most probabl tracks
  mpt<-data.frame(pr$`most probable track`, coordinates(pr$`most probable track`))
  mpt$lon360<-ifelse(mpt$lon<0,mpt$lon+360,mpt$lon)
  mpt$BirdID<-Meta$BirdID
  mpt$Band<-Meta$Band
  mpt$Tag<-Meta$Tag
  outData<-bind_rows(outData,mpt)
}

mpt<-outData

# Save compiled data as a csv
write_csv(mpt,"~/Dropbox/RLKI/RLKI_Migrations/Data/processedLocs/probGLS/RLKI_16-17_pGLS_part10000_iter200Adj5min.csv")

# Save compiled data as a compressed R file (much smaller, good for emailing)
saveRDS(mpt,"~/Dropbox/RLKI/RLKI_Migrations/Data/processedLocs/probGLS/RLKI_16-17_pGLS_part10000_iter200Adj5min.rds")

# Make some point and line plots ------------------------------------------

# function to deal with antimeridian
wrap360 = function(lon) {
  lon360<-ifelse(lon<0,lon+360,lon)
  return(lon360)
}

w2hr<-map_data('world')
w2hr_sub<-w2hr[w2hr$region%in%c("USA","Russia","China","Japan","Canada"),]
mpt$month<-month(mpt$dtime)


ggplot(mpt, aes(x=lon360, y=lat, group=as.factor(BirdID))) +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_path(aes(colour=as.factor(month))) +
  #scale_color_viridis(discrete = T) +
  coord_fixed(ratio=1.7,xlim = c(142,210),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+facet_wrap(~BirdID)

ggsave(filename = "~/Dropbox/RLKI/RLKI_Migrations/Maps/KDplots/Tracks_by_bird.jpg",scale = 2,height=8,width=10)


ggplot(mpt, aes(x=lon360, y=lat, group=as.factor(BirdID))) +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_path(aes(colour=as.factor(BirdID))) +
  #scale_color_viridis(discrete = T) +
  coord_fixed(ratio=1.7,xlim = c(142,210),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+facet_wrap(~month)

ggsave(filename = "~/Dropbox/RLKI/RLKI_Migrations/Maps/KDplots/Tracks_by_month.jpg",scale = 2,height=8,width=10)

ggplot(mpt, aes(x=lon360, y=lat, group=as.factor(BirdID))) +
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="grey60",size=0.1)+
  geom_path(aes(colour=as.factor(month))) +
  #scale_color_viridis(discrete = T) +
  coord_fixed(ratio=1.7,xlim = c(142,210),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw()+facet_grid(month~BirdID)

ggsave(filename = "~/Dropbox/RLKI/RLKI_Migrations/Maps/KDplots/Tracks_by_month_bird.jpg",scale = 2,height=8,width=10)

# Kernal Density Visualzation  ------------------------------------------------------------------
library(adehabitatHR)
library(sp)
library(adehabitatMA)

# For all bird kernel set id=1
mpt$id<-1

# Create month column
mpt$month<-month(mpt$dtime)

# Select the migration and wintering months
mpt1<-mpt[mpt$month%in% c(1:3,9:12),]

# mpt1<-mpt[mpt$month==2,]
unique(mpt1$BirdID)
unique(mpt1$month)

# create a SpatialPointsDataFrame
mptyrs.spdf <- SpatialPointsDataFrame(coords=as.data.frame(cbind(mpt1$lon,mpt1$lat)),
                                      data=mpt1, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Reproject data into laea (Lambert Azmuthal Equal Area) to preserve distances
# and areas (units are meters)
mptyrs.spdf.t <- spTransform(mptyrs.spdf,
                             CRS("+proj=laea +lat_0=50 +lon_0=-180 +x_0=2500 +y_0=1000
                                 +ellps=WGS84 +units=km +no_defs"))

# add the laea coordinates to the dataframe
head(coordinates(mptyrs.spdf.t))

mptyrs.spdf.t$x_laea<-round(coordinates(mptyrs.spdf.t)[,1],0)
mptyrs.spdf.t$y_laea<-round(coordinates(mptyrs.spdf.t)[,2],0)

# Regular grid with no land mask ---------------------------
# define the coordinates (all the detections fall inside the polygon)
# 2. Compute a grid around the fixes.
buffer_x <- as.integer((max(mptyrs.spdf.t$x_laea) - min(mptyrs.spdf.t$x_laea)) * 0.5/100) * 10
buffer_y <- as.integer((max(mptyrs.spdf.t$y_laea) - min(mptyrs.spdf.t$y_laea)) * 0.5/100) * 10
buffer <- max(buffer_x, buffer_y)
xy_sp <- SpatialPoints(data.frame(x = c((as.integer((max(mptyrs.spdf.t$x_laea) + 100)/100) * 100 + buffer),
                                        (as.integer((min(mptyrs.spdf.t$x_laea) - 100)/100) * 100 - buffer)),
                                  y = c((as.integer((max(mptyrs.spdf.t$y_laea) + 100)/100) * 100 + buffer),
                                        (as.integer((min(mptyrs.spdf.t$y_laea) - 100)/100) * 100 - buffer))))
customGrid <- ascgen(xy_sp, cellsize = 50)

# compute kernel density for all birds combined
udall <- kernelUD(mptyrs.spdf.t, h = 88,grid=customGrid)

# Extract the various utilization distribution polygons
ud25 <- getverticeshr(udall, percent=25, standardize=TRUE)
ud50 <- getverticeshr(udall, percent=50, standardize=TRUE)
ud75 <- getverticeshr(udall, percent=75, standardize=TRUE)
ud95 <- getverticeshr(udall, percent=95, standardize=TRUE)

# Transform them back into lon/lat
ud25 <- spTransform(ud25, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
ud50 <- spTransform(ud50, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
ud75 <- spTransform(ud75, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))
ud95 <- spTransform(ud95, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Plot the kernel
years<-"2015-16"
ggplot()+
  geom_polygon(data=fortify(ud95),aes(x=wrap360(long),y=lat,group=group,fill="95% UD "),show.legend=T)+
  geom_polygon(data=fortify(ud75),aes(x=wrap360(long),y=lat,group=group,fill="75% UD "),show.legend=T)+
  geom_polygon(data=fortify(ud50),aes(x=wrap360(long),y=lat,group=group,fill="50% UD "),show.legend=T)+
  geom_polygon(data=fortify(ud25),aes(x=wrap360(long),y=lat,group=group,fill="25% UD "),show.legend=T)+
  geom_polygon(data=w2hr_sub,aes(wrap360(long),lat,group=group),fill="black",color="black",size=0.1)+
  scale_fill_manual(values = rev(c("grey35","grey50","grey75","grey95")))+
  annotate("text",x=200,y=40,label=years,color="black",size=8)+
  coord_map("lambert",45,60,xlim = c(142,210),ylim=c(40,65))+
  xlab("Longitude (0-360)")+
  ylab("Latitude")+
  theme_bw(base_size = 18)+
  theme(axis.text  = element_text(color="black",hjust = 0.5),
        strip.background = element_rect(fill="white",colour="white"),legend.title = element_blank(),legend.position = "top")

ggsave(filename = "~/Dropbox/RLKI/RLKI_Migrations/Maps/KDplots/Kernel_UD_1516.jpg",scale = 2,height=8,width=10)

