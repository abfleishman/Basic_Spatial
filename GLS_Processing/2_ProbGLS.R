# This script is to prepare and  run data through the probGLS package to
# to resolve winter tracking data for Integeo GLS tags
#
# Script was modified from code downloaded from the probGLS github website by:
# Abram B. Fleishman and Rachael Orben
# 15 June 2016
#
# This script is intended for use with the Integeo C65 GLS tags set to collect
# data on setting 7 (5 minute clipped light, sum wet over five minutes checking
# every 6 seconds).
#
# NOTE::: This script must be modified if data were not collected with setting 7.


# Load (Install) Packages --------------------------------------------------------

# install.packages("devtools")
library(devtools)
# install_github("eldarrak/FLightR")
# install_github("SLisovski/TwGeos")
# install_github("benjamin-merkel/probGLS")

library(FLightR) # for the twGeos2TAGS()
library(TwGeos) # for readLig, readLux, readAct
library(probGLS)
library(readr) #for read_csv that inturprets dates for you
library(lubridate) # easy date time manipulation
library(dplyr)


# Load data ---------------------------------------------------------------

# Set the path to the twilight file
twlPath<-"/Users/abramfleishman/Dropbox/RLKI/RLKI_Migrations/Data/processedTwl/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj_twl_0.5.csv"
# Set activity path
actPath<-"/Users/abramfleishman/Dropbox/RLKI/RLKI_Migrations/Data/rawgls/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj.deg"
# Set raw light path
lightPath<-"/Users/abramfleishman/Dropbox/RLKI/RLKI_Migrations/Data/rawgls/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj.lux"


# Read in all the data files
twl<-read_csv(twlPath)
act<-read.table(actPath,skip = 19,header = T, sep="\t")
lux<-readMTlux(lightPath)

# this code is set up for loggers that did not collect temperature data
# Currently we do not have any teperature data from the 15-16 tags so I have set
# sen and tem to NULL which allows the probGLS function to work without changing it
sen<-NULL
tem<-NULL

# Add column named "wetdry" for probGLS
act<-act %>%
  mutate(dtime=as.POSIXct(DD.MM.YYYY.HH.MM.SS,format = "%m/%d/%Y %H:%M:%S",tz="GMT")) %>%
  select(dtime,wetdry=wets0.50) %>%
  filter(!is.na(dtime)) # get rid of any NAs in the activity data

# Confirm and reset timezones to GMT
tz(twl$Twilight)<-"GMT" # problems if tz=UTC so set tz=GMT to avoid any issues.
tz(act$dtime)<-"GMT" # problems if tz=UTC so set tz=GMT to avoid any issues.

# Transform in to trn format (I think this is a geolight format)
twltemp<-twl
twltemp$Twilight<-twltemp$Twilight3 #the twilight3 column has the edited twilights
twl3<-export2GeoLight(twltemp)




# OPTIONAL: Outlier removal -----------------------------------------------

# This plots twilight "outliers" k interquartile ranges from the median
# the parameter "loess.quartile" in the prob_algorithm function is the k parameter here
# this allows you to see which twilights would be excluded
# We did not use this
# a<-loessFilter(twl=twl3,k=2.5)

# see ?prob_algorithm for details

# Set some parameters for the probGLS function ----------------------------

# Longitude, latitude of the capture site
CapLoc<-c( -169.66, 56.6)

# Set wet/dry temporal sampling resolution (in seconds)
wetDryRes<-6

# twilight error distribution estimation
tw            <- twilight_error_estimation()

# Run the probGLS function  -----------------------------------------------

pr <- prob_algorithm(
  particle.number             = 10000,
  # We used 10000 (probably overkill)
  iteration.number            = 200,
  # We used 200 (probably overkill)
  trn                         = twl3,
  # This is the twilight data
  sensor                      = NULL,
  # this is temperature data from the GLS looger that has been corrected to approximate SST values
  act                         = act,
  # The activity data
  loess.quartile              = NULL,
  # We did not use this. Not used if set to NULL
  tagging.location            = CapLoc,
  # change to colony lon and lat
  tagging.date                = as.Date(min(twl3$tFirst), tz = "GMT"),
  # change to tagging date
  retrieval.date              = as.Date(max(twl3$tFirst), tz = "GMT"),
  # change to retrieval date
  wetdry.resolution           = wetDryRes,
  # change to sampling resolution of conductivity switch in sec
  sunrise.sd                  = tw,
  # from twilight_error_estimation()
  sunset.sd                   = tw,
  # from twilight_error_estimation()
  range.solar                 =  c(-7, -1),
  # This is the angle of the sun below the horizon at which twilight occurs. the defaults are c(-7,-1) and we used the default successfully for the Integeos
  speed.wet                   = c(.5 , 0.25, 1.7),
  # optimal speed, speed standard deviation and max speed allowed if logger is wet in m/s
  #4km/hr = max = 1.11m/s: Shamoun-Baranes, J., Bouten, W., Camphuysen, C.J. & BAAIJ, E. (2011) Riding the tide: intriguing observations of gulls resting at sea during breeding. Ibis, 153, 411–415.
  #6km/hr=max: Conners, M.G., Hazen, E.L., Costa, D.P. & Shaffer, S.A. (2015) Shadowed by scale: subtle behavioral niche partitioning in two sympatric, tropical breeding albatross species. Movement Ecology, 1–20.
  speed.dry                   = c(10.6, 5.3, 22.22),
  # optimal speed, speed standard deviation and max speed allowed if logger is dry in m/s
  # prediction for average ground speed of blkis
  # Elliott, K.H., Chivers, L.S. & Bessey, L. (2014) Windscapes shape seabird instantaneous energy costs but adult behavior buffers impact on offspring. Movement Ecology.
  # maximum speed allowed = 80km/hr = 22.22
  # Paredes, R., Harding, A.M.A., Irons, D.B., Roby, D.D., Suryan, R.M., Orben, R., Renner, H.M., Young, R. & Kitaysky, A.S. (2012) Proximity to multiple foraging habitats enhances seabirds’ resilience to local food shortages. Marine Ecology Progress Series, 471, 253–269.
  # currently no good rational for SD of speed.  Jsut set at half optimum speed

  sst.sd                      = 0.5,
  # SST standard deviation in degree C
  max.sst.diff                = 3,
  # max difference in SST allowed in degree C
  days.around.spring.equinox  = c(10, 10),
  # days before the Spring equinox and days after the Spring equinox. (20 March)
  days.around.fall.equinox    = c(10, 10),
  ice.conc.cutoff             = NULL,
  # ranging from 0 to 1 max pct of seaice in which the animal is believed to be
  boundary.box                = c(130, -140, 30, 75),
  #constrain the analysis to a box
  land.mask                   = T,
  # true for only use ocean

  med.sea                     = F,
  black.sea                   = F,
  baltic.sea                  = F,
  caspian.sea                 = F,

  east.west.comp              = T,
  #if TRUE apply biotrack east west movement compensation (Biotrack manual v11 page 31pp.)
  backward                    = F,
  # run algorithm from end to start
  NOAA.OI.location            = '/Users/rachaelorben/Dropbox/LandMask'
) #dir with all .nc files from NOAA (sea ice, land, mean daily sst, sd daily sst)


# Save (or load Previously saved)  output for furture use -----------------

# Save compressed rds file with the results of the prob_algorithm()
saveRDS(pr,file = "~/Dropbox/RLKI/RLKI_Migrations/Data/processedLocs/probGLS/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj_0.5_noSST_pGLS_AdjTwl.rds")

# this line reads in a previously saved probfile.
# pr<-readRDS("~/Desktop/RKLI_FinalGLSTracks_2010-16/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj_0.5_noSST_pGLS.rda")

# Plot the results and save a few graphs ----------------------------------
library(maps)

# Load world2 maps.  these maps have 0-360 longitude
data("world2MapEnv")
plot_timeline(pr = pr,degElevation = NULL,center.longitude = 180)

# transform logitude to 0-360 for plotting
lon360<-ifelse(coordinates(pr$`most probable track`)[,1]<0,
               coordinates(pr$`most probable track`)[,1]+360,
               coordinates(pr$`most probable track`)[,1])

lat<-coordinates(pr$`most probable track`)[,2]

# Save plot with most probable track
jpeg(filename =  "~/Dropbox/RLKI/RLKI_Migrations/Maps/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj_0.5_noSST_pGLS_AdjTwl.jpg",
     width = 10,height = 8,units = "in",res = 200)

  map("world2",ylim=c(35,75),xlim=c(140,220),fill=T,xlab="Lon",ylab="Lat")
  map.axes(xlab="Lon",ylab="Lat")
  lines(lon360,lat,col="cyan")
  points(lon360,lat,col="red",cex=.2)

dev.off()

# This is the code from the plot_map function which uses the world map instead of
# world2.  I have changed this to use 0-360 lon
opar <- par(mfrow = c(1, 1), mar = c(4, 4, 0, 0))
pr_df<-data.frame(pr[[1]])
pr_df2<-data.frame(pr[[2]])
pr_df$lon360<-ifelse(pr_df$lon<0,pr_df$lon+360,pr_df$lon)
pr_df2$lon360<-ifelse(pr_df2$lon<0,pr_df2$lon+360,pr_df2$lon)

coordinates(pr_df)<-cbind(pr_df$lon360,pr_df$lat)
coordinates(pr_df2)<-cbind(pr_df2$lon360,pr_df2$lat)

#  Save plot with 0-360 lon; most probable track; all possible points
jpeg(filename = "~/Dropbox/RLKI/RLKI_Migrations/Maps/RLKI2015-16GLS/R799_18Jun16_232241_B108_00774driftadj_all_points_noSST_pGLS_AdjTwl.jpg",
     width = 10,height = 8,units = "in",res = 200)

  map("world2", col = 4, lwd = 0.5,xlim =pr_df@bbox[1,],ylim = pr_df@bbox[2,])
  for (s in 1:length(unique(pr_df2$step))) {
    plot(pr_df[pr_df$step == unique(pr_df$step)[s],
               ], col = colorRampPalette(c("grey90", "grey50"))(nrow(pr_df2))[s],
         add = T, pch = 19, cex = 0.3)
  }
  map("world2", add = T, col = 4, lwd = 0.5)
  map.axes()
  map.scale(pr_df@bbox[1,1]*1.03,pr_df@bbox[2,1]*1.1,relwidth = .15,ratio = F,metric = T)
  mm2 <- pr_df2
  lines(mm2$lon360, mm2$lat, col = 1)
  points(mm2$lon360, mm2$lat, cex = .5, pch = 21, bg = colorRampPalette(c("yellow",
                                                                          "darkred"))(nrow(pr_df2)))
  mm3 <- mm2[is.na(mm2$median.sun.elev), ]
  points(mm3$lon360, mm3$lat, cex = 0.7, pch = 3)
dev.off()
par(opar)


