#' Calculate the instantanious speed at each point along on a track for each bird
#' Vector of times between points
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>
#' @param Dist vector of distances between points
#' @param Time vector of times between points

#' @return A vector of speeds (km/h) between adjacent points in an animal track
#' @examples
#' Speed(Dist=data$Dist2Point,Time=data$InterPointTime)
#' @export
#############################################################################################
# Calculate instantaious speed with a distance and a duration (set up for input meters and seconds and out put km and hours)
#############################################################################################
Speed<-function(Dist=data$Dist2Point,Time=data$InterPointTime){
  round((as.numeric(Dist)/1000)/(as.numeric(Time)/(60*60)),3)
}
