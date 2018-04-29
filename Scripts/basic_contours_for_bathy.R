# Script to add contours to bathymatry map

library(raster)
library(marmap)
library(ggplot2)

aleu <- getNOAA.bathy(130, -140, 30, 75, resolution = 15,
                      antimeridian = T,keep = T) #resolution 1.85km*27=49.9km2
# Make it a raster
bathy<-as.raster(aleu)

# Create a xyz table for ggplot
bath<-fortify(aleu)

#
# Extract the contours you want
Shore<-rasterToContour(bathy,level=c(0,-50, -100, -200, -1500))

# use fortify from the ggplot2 package to transform the SpatialLinesDataFrame
# to a long format table
Shore<-fortify(Shore)

# Plot with ggplot
ggplot()+
  geom_raster(data=bath,aes(x,y,fill=z))+
  geom_path(data=Shore,aes(long,lat,group=group))+
  scale_fill_etopo()
