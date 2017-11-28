#' AddDist2Colony Calculate
#'
#' @author Abram B. Fleishman <abram.fleishman AT sjsu.edu>

#' @param data a data frame with with latitude and longitude and colony site names
#' @param CaptureSitesData data with colony site names and colony lat long
#' @param dataLon quoted name of column in data that has longitude values
#' @param dataLat quoted name of column in data that has latitude values
#' @param SiteName a quoted string indicating what the column houses your SiteNames.  Must be the same across both your data and the CaptureSitesData
#' @param dataset a commented string indicating what dataset you are using (currently only works with a single variable)
#' @return a ggplot map frame for the data to be ploted with
#' @examples
#' AddDist2Colony(data=data,CaptureSitesData=CapSitesSel,SiteName="SiteShort")
#' @export

AddDist2Colony<-function(data=data,CaptureSitesData=CapSitesSel,SiteName="SiteShort"){
  dataOut<-vector(mode = "numeric",length = nrow(data))
  Sites<-as.character(unique(data[[SiteName]]))
  
  for(j in 1:length(Sites)){

    CapSub<-CapSitesSel[CapSitesSel[SiteName]==Sites[j],]
    dataSub<-data[data[SiteName]==Sites[j],]
    distanceVector<-Dist2Colony(data = dataSub,ColonyLat = CapSub$Lat,ColonyLong = CapSub$Lon)
    dataOut[data[SiteName]==Sites[j]]<- distanceVector
  }
  return(dataOut)
}
