#' Dist2Colony Calculate the distance of a vector of points (lat,long) to a a fixed point
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>

#' @param data a dataframe
#' @param ColonyLong a longitude value of length = 1
#' @param ColonyLat a latitude value of length = 1
#' @return a vector of bearings from adjacent points
#' @examples
#' Dist2Colony( data=data, ColonyLat=56.22234,ColonyLong=-164.45922)
#'
#' @export
Dist2Colony<-function(data,ColonyLat,ColonyLong){
  # Step 1: Two packages for mapping.####

  # Step 4a: Calculate distance from each point to the San Jorge center point.
  Point2Colony<-vector(mode = "numeric",length = nrow(data))
  for(i in 1:length(data$Latitude)){ #this is a for loop
    # This is a function to calculate distance between two points from the argosfilter package
    Point2Colony[i]<-distance(lat1 = ColonyLat,lon1 = ColonyLong ,lat2 = data$Latitude[i],lon2 = data$Longitude[i])
    Point2Colony
  }
  return(Point2Colony)
}
