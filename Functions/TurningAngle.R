#' Calculate the turning angles between points on a track for each bird
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>
#' @param data data.frame of data to be queried for  latand long values
#' @param ID quoted name of column in data that is and ID field
#' @param Bearing quoted name of column in data that has the bearing between points in an animal track

#' @return A vector of turning angles (degrees) between adjacent points in an animal track
#' @examples
#' TurningAngle(data=data,ID="File", Bearing="PointBearing")
#' @export
TurningAngle<-function(data=data,ID="File", Bearing="PointBearing"){
  dataOut<-NULL
  Birds<-unique(data[[ID]])
  for(i in 1:length(Birds)){
    Data<-data[data[[ID]]==Birds[i],]

    TurnAngle<-lag(Data[[Bearing]])-Data[[Bearing]]

    TurnAngle[TurnAngle>180&!is.na(TurnAngle)]<-TurnAngle[TurnAngle>180&!is.na(TurnAngle)]-360
    TurnAngle[TurnAngle< -180&!is.na(TurnAngle)]<-TurnAngle[TurnAngle< -180&!is.na(TurnAngle)]+360

    dataOut<-c(dataOut,TurnAngle)
  }
  return(dataOut)
}
