#' BearingFromPoint Calculate bearing between two points in a vector
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>

#' @param dataIn a dataframe
#' @param ID a quoted string indicating the column name for a unique ID key
#' @param lat a quoted string indicating the column name for longitude values
#' @param lon a quoted string indicating the column name for latitude values
#' @return a vector of bearings from adjacent points
#' @examples
#' tracks$Bering<-BearingFromPoint( dataIn=tracks,ID="File", lat="Latitude",lon="Longitude")
#' @export

BearingFromPoint<-function(dataIn=tracks,ID = "File" , lat = "Latitude", lon = "Longitude"){

  # remove bad positions
  dataIn<-as.data.frame(dataIn)

  dataIn<- dataIn[ dataIn[lon] > -180 & dataIn[lon] < 180 & dataIn[lat] > -90 & dataIn[lat] < 90,]
  dataOut<-NULL

  Birds<-unique(dataIn[[ID]])

  for(i in 1:length(Birds)){
    Data<-dataIn[dataIn[[ID]]==Birds[i],]

    BearingfromPoint<-c(round(
      geosphere::bearing(
        cbind(Data[[lon]],Data[[lat]]),
        cbind(lead(Data[[lon]]), lead(Data[[lat]]))),
      digits=1))
    dataOut<-c(dataOut,BearingfromPoint)
  }
  return(dataOut)
}
