# Basic ud and mapping
# Abram Fleishman
# 13 Oct 2017
# Prepared for Olivia T. @ SJSU

library(adehabitatHR)
library(ggplot2)
library(maps)
library(mapdata)
library(dplyr)

# Load tracks -------------------------------------------------------------
# This should be one file for all the tracks or create a data.frame with all the tracks
tracks<-read.csv("path/to/tracks.csv")

# Make the kernel for an individual track ---------------------------------

# Make a spatialPointsDataFrame from the tracking data
tracks.spdf.ind <- SpatialPointsDataFrame(coords=as.data.frame(cbind(tracks$Longitude,tracks$Latitude)),
                                      data=data.frame(id=tracks$id),
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

# This makes the UD
ud_ind <- kernelUD(tracks.spdf.ind, h = "href")
# These get out the polygons for the 25-95% utilization distributions
ud25 <- getverticeshr(ud_ind, percent=25, standardize=TRUE)
ud50 <- getverticeshr(ud_ind, percent=50, standardize=TRUE)
ud75 <- getverticeshr(ud_ind, percent=75, standardize=TRUE)
ud95 <- getverticeshr(ud_ind, percent=95, standardize=TRUE)


# Make the kernel for an the whole colony ----------------------------------

# Make a spatialPointsDataFrame from the tracking data
tracks.spdf <- SpatialPointsDataFrame(coords=as.data.frame(tracks$Longitude,
                                                           tracks$Latitude), # coords= you position data
                                      data=data.frame(id=rep(1,length=length(tracks$Latitude))), # data is a single column of bird ids saying which points for to which birds
                                      proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")) # proj4string# is the projection of the data

# Make the kernel for an the whole colony
# # there are lots of different methods for smoothing and "href" is one.  You will have to figure out which to use
ud_island <- kernelUD(tracks.spdf, h = "href")

# These get out the polygons for the 25-95% utilization distributions
ud25_island <- getverticeshr(ud_island, percent=25, standardize=TRUE)
ud50_island <- getverticeshr(ud_island, percent=50, standardize=TRUE)
ud75_island <- getverticeshr(ud_island, percent=75, standardize=TRUE)
ud95_island <- getverticeshr(ud_island, percent=95, standardize=TRUE)


# Plotting time! ----------------------------------------------------------

w<-map_data("worldHires") # load map for ggplot

# plot for individual trips ---------------------------------------------------
uds<-bind_rows(
  mutate(fortify(ud25), ud="25"),
  mutate(fortify(ud50), ud="50"),
  mutate(fortify(ud75), ud="75"),
  mutate(fortify(ud95), ud="95")
) %>% arrange(group,ud)

p1<-ggplot(data=uds)+ facet_wrap(~group)+ # this make sub-plots for each bird
  coord_fixed(xlim = c(-116,-113),ylim=c(30,32))+ # extent of your map (should be bigger than your furthest points in all dirs)
  geom_polygon(data=w,aes(long,lat,group=group),fill="black",color="black")+ # plots the islands and world
  geom_polygon(data=uds %>% filter(ud=="95"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds\
  geom_polygon(data=uds%>% filter(ud=="75"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds\
  geom_polygon(data=uds%>% filter(ud=='50'),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds\
  geom_polygon(data=uds%>% filter(ud=="25"),aes(x=long,y=lat,group=group,fill=ud),show.legend=T)+ #plots the uds\
  scale_fill_grey()+
  geom_point(aes(y=31,x=-114),pch=17,color="black",size=.8)+ # Plots you colony  location
  annotate("text",x=-114,y=-31,label="San Jorge",color="black",size=6)+ #labels your colony
  theme_bw()
  #
p1

ggsave(p1, "my_map.jpg",width = 6, height = 6,units = "in")


# plot for whole colony ---------------------------------------------------

p1<-ggplot()+
  coord_fixed(xlim = c(-116,-113),ylim=c(30,32))+ # extent of your map (should be bigger than your furthest points in all dirs)
  geom_polygon(data=w,aes(long,lat,group=group),fill="black",color="black")+ # plots the islands and world
  geom_polygon(data=fortify(ud95_island),aes(x=long,y=lat,group=group),fill="grey35",show.legend=T)+ #plots the uds
  geom_polygon(data=fortify(ud75_island),aes(x=long,y=lat,group=group),fill="grey50",show.legend=T)+
  geom_polygon(data=fortify(ud50_island),aes(x=long,y=lat,group=group),fill="grey75",show.legend=T)+
  geom_polygon(data=fortify(ud25_island),aes(x=long,y=lat,group=group),fill="grey95",show.legend=T)+
  geom_point(aes(y=31,x=-114),pch=17,color="black",size=.8)+ # Plots you colony  location
  annotate("text",x=-114,y=-31,label="San Jorge",color="black",size=6)+ #labels your colony
  theme_bw()#
p1

ggsave(p1, "my_map.jpg",width = 6, height = 6,units = "in")
