library(RandomFields)
library(ggplot2)
library(plyr)
library(dplyr)
library(spdep)
library(gstat)
library(spacetime)
library(sqldf)
library(CompRandFld)
library(fields)
library(mgcv)
library(MASS)
library(snht)

#Input:
#data: the spatio-temporal data to model.  Must have at least three columns:
#  StationID: The identifier for the location, to tie with sp.
#  Time: The time of the observation
#  Value: The observed value
#Currently, all other columns are ignored.
#Note: If NULL, the Queens' object "ground" is used.
#sp: Define a SpatialPoints object for modeling.  If NULL, defaults to stations.  If
# not NULL, it should be a SpatialPoints object with a column called "StationID"
#time: vector of times for filtering data.  If NULL, no filtering is done
#mod: A variogram model to fit to the data, as created by vgmST.  If NULL, a separable
#  exponential model is used.
#...: Additional arguments to pass to variogramST
#Output:
#Spatio-temporal variogram
fitModel = function(data=NULL, sp=NULL, time=NULL, mod=NULL, ...){
  #data
  if(!is.null(data)){
    stopifnot(is(data, "data.frame"))
    stopifnot(all(c("StationID", "Time", "Value") %in% colnames(data)))
  }
  #sp
  if(!is.null(sp)){
    stopifnot("StationID" %in% names(sp))
    stopifnot(is(sp,"SpatialPoints"))
  }
  #time
  #mod
  if(!is.null(mod)){
    stopifnot(is(mod,"StVariogramModel"))
  }
  
  #If sp is NULL, load the appropriate data and define sp as station locations
  if(is.null(sp)){
    if(!exists("station"))
      load("Data/station_new_coords.RData")
    sp = station
    sp = sp[!is.na(sp$Longitude),]
    sp = sp[,c("StationID", "East2", "North2")]
    colnames(sp) = c("StationID", "Easting", "Northing")
    coordinates(sp) = c("Easting", "Northing")
  }
  
  if(is.null(time))
    time = unique(data$Time)

  if(is.null(data)){
    if(!exists("ground"))
      load("Data/ground_with_distance.RData")
    data = filter(ground, Time %in% time)
  }
  
  data = filter(data, StationID %in% sp$StationID)
  if(nrow(unique( data[,c("StationID", "Time")]))!=nrow(data))
    stop("data has multiple observations for one StationID-Time pair!")
  d = expand.grid(StationID=sp$StationID, Time=time)
  d = merge(d, data, by=c("StationID", "Time"), all.x=TRUE)
  d = arrange(d, Time, StationID)
  stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
  vst = variogramST( Value ~ 1, data=stdf, ... )
  
  #Fit the variogram model
  if(is.null(mod))
    mod = vgmST("separable"
           ,space=vgm(psill=.90,"Exp", range=600, nugget=0.4),
           ,time =vgm(psill=.90,"Exp", range=40, nugget=0.4),
           ,sill=0.9)

  vstModel = fit.StVariogram(vst, model=mod)
  return(vstModel)
}

#Input:
#model: a spatio-temporal variogram, as returned by fitModel
#...?
#Output:
#Simulated data
simulate(vst){
  time = unique(ground$Time)
  grid = as.matrix(station[,c("East2", "North2")])[!is.na(station$East2),]
  data = RFsim(coordx=grid, corrmodel="exponential"
      ,grid=FALSE, param=list(nugget=.693, mean=0, sill=1, scale=1e10) )
  data$data
}

#Input:
#d: output from simulate()
#randPct: % of random errors
#sysPct: % of systematic errors
#...?
#Output:
#d with contaminated observations.  Additionally, output the original data
#so analysis of model performance will be easier.
contaminate(d, randPct, sysPct){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via SNHT applied to each individual
#station.
univariateSNHT(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via pairwise SNHT.
pairwiseSNHT(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via the isolation metric defined in
#Filzmoser's paper: Identification of local multivariate outliers (2014).
#Note: this method may also detect random errors.
isolation(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a homogenization
#algorithm.
#...?
#Output:
#Identification of random errors (probably a vector of scores)
outlier(d){
  
}