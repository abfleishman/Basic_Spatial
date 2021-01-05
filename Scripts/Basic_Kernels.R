# Basic UD and mapping
# Abram Fleishman & Morgan Gilmour
# 4 Jan 2021
# Prepared for Scott S. @ SJSU

library(adehabitatHR)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)

# Load the scaleARS() function 
# Modified from Lascalles et al. 2016 
# Citation: Lascelles BG, Taylor PR, Miller MGR, Dias MP, Oppel S, Torres L, Hedd A, Le
# Corre M, Phillips RA, Shaffer SA, et al. 2016. Applying global criteria to
# tracking data to define important areas for marine conservation. Divers
# Distribut. 22(4):422–431. doi:10.1111/ddi.12411.

scaleARS <- function(DataGroup=f, Latitude="lat",Longitude="Lon_360", ID = "CaptureID",TrackTime="datetime",
                     Scales = c(seq(1, 25, 1), seq(30, 50, 5), 75, seq(100, 250, 50)), Peak = "Flexible")
{
  
  
  require(geosphere)
  require(sp)
  require(rgdal)
  require(rgeos)
  require(adehabitatLT)
  
  if(!Latitude %in% names(DataGroup)) stop("Latitude field does not exist")
  if(!Longitude %in% names(DataGroup)) stop("Longitude field does not exist")
  if(!ID %in% names(DataGroup)) stop("ID field does not exist")
  if(!TrackTime %in% names(DataGroup)) stop("TrackTime field does not exist")
  
  if(class(DataGroup)!= "SpatialPointsDataFrame")     ## convert to SpatialPointsDataFrame and project
  {
    mid_point<-data.frame(centroid(cbind(DataGroup[Longitude], DataGroup[Latitude])))
    DataGroup.Wgs <- SpatialPoints(data.frame(DataGroup[Longitude], DataGroup[Latitude]), 
                                   proj4string=CRS("+proj=longlat +datum=WGS84 +lon_wrap=180"))
    
    DgProj <- CRS(paste("+proj=laea +lon_0=", mid_point$lon, " +lat_0=", mid_point$lat, sep=""))
    DataGroup.Projected <- spTransform(DataGroup.Wgs, CRS=DgProj)
    DataGroup <- SpatialPointsDataFrame(DataGroup.Projected, data = DataGroup)
  }else{DgProj<-DataGroup@proj4string}
  
  DataGroup$X <- DataGroup@coords[,1]
  DataGroup$Y <- DataGroup@coords[,2]
  
  
  DataGrouplt <- as.ltraj(data.frame(DataGroup$X, DataGroup$Y), date=as.POSIXct(DataGroup[[TrackTime]], origin="1970/01/01", tz="GMT"), id=DataGroup[[ID]], typeII = TRUE)
  
  Scales <- Scales * 1000
  
  fpt.out <- fpt(DataGrouplt, radii = Scales, units = "seconds")
  fpt.scales <- varlogfpt(fpt.out, graph = FALSE)
  Temp <- as.double(fpt.scales[1,])
  plot(Scales, Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)))
  
  ars.scales <- NULL
  UIDs <- unique(DataGroup[[ID]])
  for(i in 1:length(UIDs))
  {
    if(length(Scales) == length(which(is.na(fpt.scales[i,])))) {print(paste("Warning: ID", UIDs[i], "is smaller than smallest scale and will be ignored")); next}
    Temp <- as.double(fpt.scales[i,])
    #lines(Scales,Temp)
    plot(Scales, Temp, type="l")
    
    q <- which(!is.na(Temp))
    p <- 2
    while(!is.na(Temp[q[p]]) & Temp[q[p]] < Temp[q[p-1]] & q[p] != length(Temp)) {p <- p + 1}
    while(!is.na(Temp[q[p]]) & Temp[q[p]] > Temp[q[p-1]]) {p <- p + 1}
    
    rfpt <- Scales[q[p-1]]
    if(suppressWarnings(min(which(is.na(Temp))) == p)) {print(paste("ID", UIDs[i], "has no peak")); next}
    FirstPeak <- Scales[q[p-1]]
    MaxPeak <- Scales[which(Temp == max(Temp[q[p-1]:length(Temp)], na.rm=T))]
    if(Peak == "Flexible")
    {
      if(FirstPeak < MaxPeak[1])
      {
        MaxPeak <- MaxPeak[MaxPeak >= FirstPeak]
        ifelse(MaxPeak[1] < FirstPeak + (max(Scales)/3), ars.sc <- MaxPeak[1], ars.sc <- FirstPeak)
      }  else  {ars.sc <- FirstPeak}
    }
    if(Peak == "Max") {ars.sc <- MaxPeak}
    if(Peak == "First")  {ars.sc <- FirstPeak}
    if(Peak == "User")
    {
      print("Select Peak on Graph")
      N <- identify(Scales, Temp, n=1)
      ars.sc <- Scales[N]
    }
    abline(v=ars.sc, col="red", lty=2)
    ars.scales <- c(ars.scales, ars.sc)
    #print(ars.sc)
    #readline("proceed?")
  }
  
  AprScale <- mean(ars.scales)
  AprScale <- round(AprScale/1000,3)
  plot((Scales/1000), Temp, type="l", ylim=c(0, max(fpt.scales, na.rm=T)), xlab="Scales (km)", ylab="")
  for(i in 1:length(UIDs))
  {
    Temp <- as.double(fpt.scales[i,])
    lines((Scales/1000),Temp)
  }
  abline(v=ars.scales/1000, col="red", lty=2)
  abline(v=AprScale, col="darkred", lty=1, lwd=3)
  #print(ars.scales)
  #print(AprScale)
  text(max(Scales/1000)/2, 1, paste(AprScale, "km"), col="darkred", cex=3)
  return(AprScale)
}


# Load tracks -------------------------------------------------------------
# This should be one file for all the tracks or create a data.frame with all the tracks
tracks<-read.csv("path/to/tracks.csv")
#
# Make a kernel for an individual track ---------------------------------

# Make a spatialPointsDataFrame from the tracking data
tracks.ind<-tracks %>% filter(tripID=="IDnum") #Filter tracks to just one bird or trip
tracks.spdf.ind <- SpatialPointsDataFrame(coords=data.frame(tracks.ind$Longitude,tracks.ind$Latitude),
                                      data=data.frame(id=tracks.ind$tripID),
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# Kernel analyses

# Determine the smoothing parameter, h
# Several options to determine the scale at which animals search.
# 1) h can either be a default value from the function h="href";
# 2) h can be determined via a least squares cross validation h="LSCV";
# 3) h can be a numeric value determined by the user. Sometimes published studies
# report the h-value that they used, and that can be used here.
# 4) h can be determined via the function scaleARS() from Lascelles et al. (2016)
scale_animal<-scaleARS(DataGroup = tracks.ind,Latitude = "Latitude",Longitude = "Longitude",
                   ID = "tripID",TrackTime = "datetime")

# 
# This makes the UD
ud_ind <- kernelUD(tracks.spdf.ind, h = "href")
# These get out the polygons for the 25-95% utilization distributions
ud25 <- getverticeshr(ud_ind, percent=25, standardize=TRUE)
ud50 <- getverticeshr(ud_ind, percent=50, standardize=TRUE)
ud75 <- getverticeshr(ud_ind, percent=75, standardize=TRUE)
ud95 <- getverticeshr(ud_ind, percent=95, standardize=TRUE)

#
# Make a kernel for the whole colony ----------------------------------

# Make a spatialPointsDataFrame from the tracking data
tracks.spdf <- SpatialPointsDataFrame(coords=data.frame(tracks$Longitude,
                                                           tracks$Latitude), # coords= your location data (lat/long).
                                      data=data.frame(id=rep(1,length=length(tracks$Latitude))), # To make UD for whole colony, tell the function the IDs are all the same.
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # proj4string# is the projection of the data.

# Make the kernel for the whole colony
# See notes about different methods for smoothing and "href" above.
ud_island <- kernelUD(tracks.spdf, h = "href")

# These get out the polygons for the 25-95% utilization distributions
ud25_island <- getverticeshr(ud_island, percent=25, standardize=TRUE)
ud50_island <- getverticeshr(ud_island, percent=50, standardize=TRUE)
ud75_island <- getverticeshr(ud_island, percent=75, standardize=TRUE)
ud95_island <- getverticeshr(ud_island, percent=95, standardize=TRUE)

#
# Plotting time! ----------------------------------------------------------

w<-map_data("world") %>% rename(world_group=group) # load map for ggplot
#
# plot for individual trips ---------------------------------------------------
uds<-bind_rows(
  mutate(fortify(ud25), ud="25"), #fortify() makes polygons of the UD objects so that they can be plotted with ggplot
  mutate(fortify(ud50), ud="50"),
  mutate(fortify(ud75), ud="75"),
  mutate(fortify(ud95), ud="95")
) %>% arrange(group,ud)

p1<-ggplot(uds)+ 
  coord_fixed(xlim = c(-116,-113),ylim=c(30,32))+ # extent of your map (should be bigger than your furthest points in all directions)
  geom_polygon(data=w,aes(long,lat,group=world_group),fill="black",color="black")+ # plots the islands and world
  geom_polygon(data=uds%>% filter(ud=="95"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds
  geom_polygon(data=uds%>% filter(ud=="75"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds
  geom_polygon(data=uds%>% filter(ud=='50'),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds
  geom_polygon(data=uds%>% filter(ud=="25"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds
  facet_wrap(~id)+ # this makes sub-plots for each bird, if applicable (i.e. trips for the same bird)
  scale_fill_grey()+
  geom_point(aes(y=31,x=-114),pch=17,color="black",size=0.8)+ # Plots your colony location
  annotate("text",x=-114,y=-31,label="Isla San Jorge",color="black",size=6)+ #labels your colony
  theme_bw()
  #
p1

ggsave(p1, "my_map.jpg",width = 6, height = 6,units = "in")

#
# plot for whole colony ---------------------------------------------------

p1<-ggplot()+
  coord_fixed(xlim = c(-116,-113),ylim=c(30,32))+ # extent of your map (should be bigger than your furthest points in all directions)
  geom_polygon(data=w,aes(long,lat,group=world_group),fill="black",color="black")+ # plots the islands and world
  geom_polygon(data=fortify(ud95_island),aes(x=long,y=lat,group=group),fill="grey35",show.legend=T)+ #plots the uds
  geom_polygon(data=fortify(ud75_island),aes(x=long,y=lat,group=group),fill="grey50",show.legend=T)+
  geom_polygon(data=fortify(ud50_island),aes(x=long,y=lat,group=group),fill="grey75",show.legend=T)+
  geom_polygon(data=fortify(ud25_island),aes(x=long,y=lat,group=group),fill="grey95",show.legend=T)+
  geom_point(aes(y=31,x=-114),pch=17,color="black",size=0.8)+ # Plots your colony  location
  annotate("text",x=-114,y=-31,label="Isla San Jorge",color="black",size=6)+ #labels your colony
  theme_bw()#
p1

ggsave(p1, "my_map.jpg",width = 6, height = 6,units = "in")
