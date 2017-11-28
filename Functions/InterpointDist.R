#' Calculate the Distance between points on a track for each bird
#' Vector of distances between points
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>
#' @param data data.frame of data to be queried for  latand long values
#' @param ID quoted name of column in data that is and ID field
#' @param lat quoted name of column in data that has latitude values
#' @param lon quoted name of column in data that has longitude values

#' @return A vector of distances in meters between adjacent points in an animal track
#' @examples
#' InterpointDist(data,ID="File",lat="Latitude",lon="Longitude")
#' @export
#'
#############################################################################################
#### Calculate the Distance between points on a track for each bird                      ####
#### Vector of distances between points
#############################################################################################
InterpointDist<-function(data,ID="File",lat="Latitude",lon="Longitude"){
  # Step 1: Two packages for mapping.####
  library(argosfilter)

  # Step 4a: Calculate distance from each point to the San Jorge center point.
  dataOut<-NULL
  Birds<-unique(data[[ID]])
  for(i in 1:length(Birds)){ #this is a for loop
    Data<-data[data[[ID]]==Birds[i],]
    # This is a function to calculate distance between two points from the rgeos package
    InterpointDist<-c(NA,round(distanceTrack(lat = Data[[lat]],lon = Data[[lon]])*1000,digits=1))
    dataOut<-c(dataOut,InterpointDist)
  }
  return(dataOut)
}
