#Morphological strategies of corals in the anthropocene: functions

rm(list=ls())

logit <- function(p) { log(p/(1-p))}
invLogit <- function(p) {exp(p)/(1+exp(p)) }

#Bearing from lat and long

rad2Degree <- function(Rad) (Rad*180)/pi
degree2Rad <- function(Degree) Degree*(pi/180)

bearingLatLong <- function(Lat1, Long1, Lat2, Long2) {
  Lat1 <- degree2Rad(Lat1)
  Long1 <- degree2Rad(Long1)
  Lat2 <- degree2Rad(Lat2)
  Long2 <- degree2Rad(Long2)
  Y <- cos(Lat2) * sin(Long2-Long1)
  X <- cos(Lat1) * sin(Lat2) -
    sin(Lat1) * cos(Lat2) * cos(Long2-Long1)
  Bearing <- atan2(Y, X)
  Bearing <- rad2Degree(Bearing)
  Bearing <- (Bearing + 360) %% 360
  return(Bearing)
}

#rescale

rescale <- function(Dat, Min = NA, Max = NA) {
  if(is.na(Min)) Min <- min(Dat)
  if(is.na(Max)) Max <- max(Dat)
  return((Dat - Min)/(Max - Min))
}


#multi join list of dataframes by key

multiJoin <- function(DFList) Reduce(function(...) merge(..., by = c("site","campaign")),DFList)


#square root transforms -inf to +inf data

squareRootifier <- function(Data) {
  sapply(1:length(Data), function(x) {
    if(Data[x] >=0) {
      sqrt(Data[x])
    } else {
      -(sqrt(abs(Data[x])))
    }
  })
}


