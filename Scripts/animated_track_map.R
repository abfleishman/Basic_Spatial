
# Script animate gps tracks from gps loggers
# 
# Made for Caitie Kroeger
#
# author: Abram Fleishman
# Updated 22 Feb 2019
# 
library(ggplot2)
library(dplyr)
library(readr)
# install.packages("gganimate")
# install.packages("gifski")
# install.packages("transformr")
library(gifski)
library(transformr)
library(gganimate)
library(lubridate)

# Read in data and set Band as a factor (we use reader to parse datetime)
tracks<-read_csv("C:/Users/ConservationMetrics/Downloads/Kroeger_Tracks_exclusions_2jan18.csv") %>% 
  mutate(Band=factor(Band)) %>% 
  group_by(Band) %>% 
  # I had to make a bird step so that all the birds would happen at the same
  # time.  I think you could do some other things instead like faceting by a
  # tracking session (all the 2011 NOV/Dec birds together for instance.  )
  mutate( bird_step=as.numeric(timestamp-min(timestamp))) 

head(tracks) %>% as.data.frame

sort(unique(as.Date(tracks$timestamp)))
ggplot(tracks ,aes(x=lon, y=lat, color=Band))+
  # add the points
  # geom_point(aes(group = bird_step),show.legend = F)+
  # add the path
  geom_path(show.legend = F)+ 
  # point for Campbell Island
  geom_point(x=169.4207,y=-52.79334,size=5,pch=24,fill="yellow",show.legend = F)+ 
  # label for Campbell Island
  geom_text(x=169.4207,y=-52.79334,label="campbell\nIsland",color="black",hjust=-.1,show.legend = F)+
  # Here comes the gganimate specific bits
  labs(title = 'Date: {frame_along}', x = 'Longitude', y = 'Latitude') +
  # slowly reveal the points and path
  transition_reveal(bird_step) +
  # change this to other options for soother transitions?  I am not sure how to use this 
  ease_aes('linear')+
  theme_classic()+
  coord_fixed()
# save.  I am not sure how to set the size... might be based on your "viewer size at the time of creation
anim_save( filename = "Kroeger_Tracks.gif",path ="C:/Users/ConservationMetrics/Downloads",)
