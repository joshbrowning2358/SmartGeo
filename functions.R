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

#Output:
#Generates a SpatialPoints object with the coordinates of the Queens' stations
stationSp = function(){
  if(!exists("station")){
    if(Sys.info()[4]=="JOSH_LAPTOP")
      load("~/Documents/Professional Files/Mines/SmartGeo/Queens/Data/station_new_coords.RData")
    if(Sys.info()[4]=="jb")
      load("/media/storage/Professional Files/Mines/SmartGeo/Queens/Data/station_new_coords.RData")
  }
  sp = station
  sp = sp[!is.na(sp$Longitude),]
  sp = sp[,c("StationID", "East2", "North2")]
  colnames(sp) = c("StationID", "Easting", "Northing")
  coordinates(sp) = c("Easting", "Northing")
  return(sp)
}

#Example:
#Output: returns a list meant to look like an StVariogramModel object.  Useful on Linux,
#where we can't fit variograms because of missing packages.
vstMod = function(){
  list(space=data.frame(model=c("Nug", "Exp")
                             ,psill=c(0.2197915,0.7802085)
                             ,range=c(0,609.344) )
            ,time=data.frame( model=c("Nug", "Exp")
                             ,psill=c(0.6188773,0.3811227)
                             ,range=c(0,93.15087) )
            ,sill=0.01647729
            ,stModel="separable")
}

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
fitModel = function(data=ground, sp=stationSp(), time=unique(ground$Time)[1:100], mod=NULL, ...){
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
  
  if(is.null(data)){
    if(!exists("ground"))
      load("Data/ground_with_distance.RData")
    data = ground
  }
  
  if(is.null(time))
    time = unique(data$Time)[1:100]
  
  data = filter(data, StationID %in% sp$StationID, Time %in% time)
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

#fitModel doesn't work in Linux (package issues).  So, manually copy over variogram.
#Input:
#model: a spatio-temporal variogram, as returned by fitModel
#...?
#Output:
#Simulated data
simulate(vst=vstMod(), time=ground$Time[1:100], sp=stationSp()){
  grid = as.matrix(data.frame(sp)[,c("Easting", "Northing")])
  data = RFsim(coordx=grid, corrmodel="exp_exp", grid=FALSE
          ,param=list(
            nugget_s= vst$space$psill[vst$space$model=="Nug"]
           ,scale_s = vst$space$range[vst$space$model=="Exp"]
           ,nugget_t= vst$space$psill[vst$space$model=="Nug"]
           ,scale_t = vst$time$range[vst$time$model=="Exp"]
           ,sill    = vst$sill) )
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