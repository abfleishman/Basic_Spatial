# Geolocation analysis with Open Source Tools
# 2016 North American Ornithological Congress, Washington D.C.
# Downloaded from:
#
# Sponsored by:
# Migrate Technology LLC.-- www.migratetech.co.uk
# The Cooper Ornithological Society
# The National Science Foundation

# Edited by Abram Fleishman 22 Dec 16 for use with
# GLS tag data from RLKI on St.George ISland, AK


# Read in data from a .lux file from a Migrate Technology tag  ------------

#First reduce the data down to just datestamps and light levels
#use the readMTlux function in TwGeos to read in the data

#to install TwGeos
#library(devtools)
#install_github("SLisovski/TwGeos")

library(TwGeos)
library(lubridate)
library(tidyr)

# clear your environment
rm(list=ls())

# So that I can just go through the list of files and keep track of where I am,
# I define a TagNum and change that as IF I were using a loop.  But I am not
# since there are manual steps

#ABF: RLKI2015-16GLS  threshold: 0.5

# list all the lux files (we used the driftadj file)
# I decided to process all of them so we can evaluate whether or not to use the driftadj files

# path to the directory with all the raw data
Migration<-"~/Dropbox/RLKI/RLKI_Migrations/Data/rawgls/RLKI2015-16GLS"

# List the files
FilesLux<-list.files(Migration,pattern="driftadj.lux",full.names = T)
FilesLux

# We are going to process one file at a time manually to deterime twilights.
# To do this we will go through each file (TagNum)
TagNum=1
Band<-"173300774"
# Or you can extract the band from the file name but depends of how you named them
# Tag<-unlist(strsplit(basename(FilesLux[TagNum]),split = "_"))
#
# Tag
#
# Tag[5]<-strsplit(Tag[5],"drift")[[1]] # Warning OK

#read the data into a dataframe called d.lux
d.lux<-readMTlux(FilesLux[TagNum])


# Read in Handling data  and subset for deployment/retrieval ----------------
# Handling data must havethe following columns:
#   Band: 5 digit suffix
#   Deploy_AK_datetime: format= mm/dd/yyyy HH:MM
#   Retrieval_AK_datetime: format= mm/dd/yyyy HH:MM
#
# Read in handling data
hand<-read.csv(file.path("~/Dropbox/RLKI/RLKI_Migrations/Data/rawgls/RLKI_deploymentmatrix_2010-2016.csv"),
               stringsAsFactors = F )

# Subset hadling data for the correct tag
hand<-hand[hand$Band==Band,]
hand # make sure it is only a single row and that the deploy and retrive times are not NA
# hand$Deploy_AK_datetime<-"08-04-2015 8:00" # if needed manually set deploy and retrive times

# subset light data to on bird data only correcting for the tag data in GMT
d.lux<-d.lux[d.lux$Date>with_tz(mdy_hm(hand$Deploy_AK_datetime,tz="America/Anchorage"),tzone = "GMT")&
               d.lux$Date<with_tz(mdy_hm(hand$Retrieval_AK_datetime,tz="America/Anchorage"),tzone = "GMT"),]

head(hand)
head(d.lux)         #view the first few lines
max(d.lux$Light)    #note the maximum value
min(d.lux$Light)    #note the minimum value


# Plot the light data -----------------------------------------------------

#Lets view the light data.
#You can use the plot function to look at small pieces of the dataset.
dev.off()
plot(d.lux$Date[7000:8000], d.lux$Light[7000:8000], type = "l")

#In lux files light values go very high, so we can log transform data before selecting twilights.
# this is handy because it increases the contrast between light and dark values
# and allows more precision with determining twilights
d.lux$Light<-log(d.lux$Light)

#Now try the simple plot again
plot(d.lux$Date[7000:8000], d.lux$Light[7000:8000], type = "l")

#For a more complete view use the lightimage() function in the TwGeos package.
#In this graph each vertical line is a day (24 hours) of data.
#Low light levels are shown in dark shades of gray and high levels are light shades.
#The offset value of 12 adjusts the y-axis to put night (dark shades) in the middle.
lightImage(d.lux, offset = 0, zlim = c(0, max(d.lux$Light)), dt = 300)
#Note the dark areas near the beginning and end of the dataset.
#These are probably associated with nesting and the pre/post deployment period. <- Incubation!

#------------------------------------------``
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\``
#------------------------------------------``

#Options for editing twilights.
#   preprocessLight() from the TwGeos package
#   findTwilights() from the TwGeos package
#   twilightCalc() from the GeoLight package
#   TAGS - a web-based editor

#Establish a threshold for determining twilights (what light value separates day and night?)
#The best choice is the lowest value that is consistently above any noise in the nighttime light levels
#For this log transformed LUX tag a good choice appears to be 1.5 which is about the same as log(4.5)
#threshold = log(4.5)

threshold = 0.5 # prefered as it catches little light spikes This is what We used
#------------------------------------------``
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
#------------------------------------------``

# Process light data ------------------------------------------------------


## preprocessLight() is an interactive function for editing light data and deriving twilights
## Unfortunately, it does not work with an RStudio web interface (i.e. a virtual machine)
## Note: if you are working on a Mac set gr.Device to X11 and install Quartz first (https://www.xquartz.org)
## See help file for details on the interactive process.

M<-max(d.lux$Light)
help(preprocessLight)

## for pc
# twl <- TwGeos::preprocessLight(d.lux, threshold = threshold, offset = 0, lmax = M, gr.Device = "default")

## for mac
twl <- TwGeos::preprocessLight(d.lux, threshold = threshold, offset = 0, lmax = M, gr.Device = "x11")

# step 1: trim to deployment, regular click on beginning, 'a' to accept
#
# step 2: find twis, regular click on beginning to get blue and orange dots, zoom out with'-' to see more
# days at once, skim days, click on those to come back to, 'a' to accept and advance
#
# step 3: edit twis that are marked with green dot
#
# step 4: skim light curves, 'q' to accept

head(twl,20)
# save the Twl file in a new folder with _twl and threshold used at the end of file name
OutFolder<-"~/Dropbox/RLKI/RLKI_Migrations/Data/rawgls/RLKI2015-16GLS/Data/processedTwl/RLKI2015-16GLS"
path1<-paste0(gsub("\\.lux",'',file.path(OutFolder,basename(FilesLux[TagNum]))),"_twl_",threshold,".csv")
write.csv(twl, file = path1, quote = FALSE, row.names = FALSE)
