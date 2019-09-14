# Generate random points while using land mask
# 
# For Luke Halpin
# Date writen: 14 Sep 2019

library(raster)
library(dismo)

# download a landmask (raster with values that can identify the water)
# here is one: ftp://ftp.cdc.noaa.gov/Datasets/noaa.oisst.v2.highres/lsmask.oisst.v2.nc
# 
# In this case 1 == water; 0 == land
mask<-raster("C:/Users/ConservationMetrics/Downloads/lsmask.oisst.v2.nc")

# set the land as NA.  The dismo package needs a layer with NAs to identify
# places not to make a point
mask[mask==0]<-NA

# limit to the extent that you need (I made this up)
mask_crop<-crop(mask,extent(130,210,-55,-10))

# then use randomPoints from the dismo package to generate points
bg <- dismo::randomPoints(mask = mask_crop, n = 1000 )
plot(bg)
