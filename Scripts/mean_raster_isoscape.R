# Make a raster with mean values per grid cell for some type of data.
#
# In this example I am making a raster with mean d13C in each grid cell for
# my red-legged kittiwakes in the western north pacific. You must start with a
# dataframe with three columns, spatail coordinants, and some values you want the
# mean of.
#
# Abram Fleishman
# 1 Apr 2018
# San Jose State University


library(ggplot2)
library(gridExtra)
library(raster)
library(dplyr)
# start with a df
d13C_data<-read.csv("")

# Make a spatial Points Dataframe
tracks_sp <- SpatialPointsDataFrame(coords=cbind(d13C_data$medlon360,d13C_data$medlat),
                                    data=as.data.frame(d13C_data),
                                    proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +lon_wrap=180"))

# Set the grid cell resolution for your out raster (we are in degrees in this example)
res<-1.25

# Make a raster with the same extent as as your SPDF above with a buffer in the
#  units that the spdf is in
e<-raster::extent(raster::buffer(tracks_sp,width=res))

# Make a raster from the extent with a resolution in the units you have
aa<-raster::raster(e,resolution=res)

# make turn the points into a raster with the mean of all the points in each cell
d13Cb<-raster::rasterize(x =tracks_sp[!is.na(tracks_sp[["d13C"]]),],
                         y = aa,
                         field="d13C",
                         fun=function(x,...) mean(x,na.rm=T))

# Make a raster with the count of birds in each cell
d13Cb_n<-raster::rasterize(x =tracks_sp[!is.na(tracks_sp[["d13C"]]),],
                           y = aa,
                           field="d13C",
                           fun=function(x,...) length(x))

# change to a dataframe for ggplot (this might explode if you have a ton of cells)
d13Ca1<-raster::as.data.frame(d13Cb,xy=T)
d13Ca1n<-raster::as.data.frame(d13Cb_n,xy=T)


grid.arrange(
ggplot()+
  coord_map(projection = "azequalarea",
            orientation = c(45,170,1),
            xlim=c(140,180),
            ylim=c(40,60))+
  geom_polygon(data=map_data("world2"),
               aes(x=long,y=lat,group=group),
               col="grey")+
  geom_tile(data=d13Ca1 %>% filter(),
            aes(x=x,y=y,fill=layer))+
  geom_text(data=d13Ca1n,aes(x=x,y=y,label=layer),size=2)+
  scale_x_continuous(breaks=c(140,150,160,170,180,190,200),
                     labels=c(140,150,160,170,180,-170,-160))+
  scale_colour_gradientn(colours = rev(rainbow(n=8,end = .7)))+
  scale_fill_gradientn(na.value = NA,
                       colours = rev(rainbow(n=8,end = .7)))+
  guides(fill = guide_colorbar(title = expression(paste(delta^{13}, "C")),
                               barwidth = 18,
                               barheight = 1,
                               direction = 'horizontal'))+
  labs(x="Longitude",
       y="Latitude")+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color="grey80",size = .5)),

# plot the points for comparison
ggplot(a ,aes(y=medlat,x=medlon360,col=d13C))+
  geom_polygon(data=map_data("world2"),aes(x=long,y=lat,group=group),col="grey")+
  geom_point(size=5)+
  scale_color_gradientn(colours = rev(rainbow(10,end=0.8)))+
  guides(color=guide_colorbar(barwidth = 15))+
  coord_map(projection = "azequalarea",orientation = c(45,170,1),xlim=c(140,202),ylim=c(40,67))+
  scale_x_continuous(breaks=c(140,150,160,170,180,190,200),labels=c(140,150,160,170,180,-170,-160))+
  theme_classic(base_size = 18)+
  theme(legend.position = "top",
        legend.direction = "horizontal",
        axis.ticks = element_blank(),
        panel.grid.major = element_line(color="grey80",size = .5))+
  labs(color=expression(paste(delta^{13}, "C")), x = "Longitude", y="Latitude"),
nrow=1)
