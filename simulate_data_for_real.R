library(nnet)
library(reshape)
library(mvtnorm)
if(Sys.info()[4]=="jb"){
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
  source("/media/storage/Github/SmartGeo/functions.R")
}
if(Sys.info()[4]=="JOSH_LAPTOP"){
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
  source("~/Github/SmartGeo/functions.R")
}
if(grepl("ch120", Sys.info()[4])){
  setwd("~/SmartGeo/Queens")
  source("~/Github/SmartGeo/functions.R")
}

#' @param spatialV The object containing a list of 3 elements (for each period).  Each
#' element is a list of length two, where the first element is the empirical variogram
#' and the second is the list of parameters to pass to modSinExp.
#' @param nnetMod The neural network model fit to the data, from which the errors were
#' computed.
#' @stationPairs A data.frame with 4 columns: s1, s2, dist, and theta.
#' @covariate A data.frame with all the data required to predict for the nnetMod
#' @station The stations corresponding to each observation in covariate
#' @time The time values corresponding to each observation in covariate
#' @randErrorPct The percent of observations that should have random errors added to their
#' value
#' @sysErrorThresh Systematic errors are created by introducing a bias in a sensor after
#' that sensor has been offline for a period of time (at least three days).  This is done
#' because it seems to be how systematic errors generally occur in the data.  To determine
#' systematic error, a statistic similar to the SNHT can be used:
#' (mean period 1 - mean period 2)/s.e.(estimator)
#' This statistic should be N(0,1).  The argument to this function specifies what value
#' this statistic should take to be considered an error (i.e. 5 => all statistics larger
#' in absolute value than 5 are assumed errors) and the number and size of systematic
#' errors are then determined by this threshold.  Specifying 1e5 (or larger) for this
#' argument forces no systematic errors to occur.
simulate = function(spatialV, nnetMod, stationPairs, covariate, station, time
        ,randErrorPct=0
        ,sysErrorThresh=1e5){
  stopifnot(length(spatialV)==3)
  stopifnot(all(lapply(spatialV, length)==2))
  stopifnot(is(nnetMod,"nnet"))
  stopifnot(all(colnames(stationPairs)==c("s1", "s2", "dist", "theta")))
  #Since we cast with s1 ~ s2 and want a symmetric matrix, we need to ensure s1 is of
  #the same type as s2
  stopifnot(is(stationPairs$s1)==is(stationPairs$s2))
  stopifnot(randErrorPct<=1)
  stopifnot(randErrorPct>=0)
  stopifnot(sysErrorThresh>=0)
  
  #Estimate the expected values from nnet
  noError = predict(nnetMod, newdata=covariate)
  noError = data.frame(station=station, time=time, model=noError[,1])
  
  #Create the three covariance matrices
  stationPairs$Cov1 = spatialV[[1]][[2]][1] - modBesExp(spatialV[[1]][[2]], stationPairs$dist, 0)
  stationPairs$Cov2 = spatialV[[2]][[2]][1] - modBesExp(spatialV[[2]][[2]], stationPairs$dist, 0) 
  stationPairs$Cov3 = spatialV[[3]][[2]][1] - modBesExp(spatialV[[3]][[2]], stationPairs$dist, 0)

  Sigma1 = cast(stationPairs, s1 ~ s2, value="Cov1" )
  rownames(Sigma1) = Sigma1[,1]
  Sigma1$s1 = NULL
  #image(as.matrix(Sigma1))
  #eigen(as.matrix(Sigma1))$values

  Sigma2 = cast(stationPairs, s1 ~ s2, value="Cov2" )
  rownames(Sigma2) = Sigma2[,1]
  Sigma2$s1 = NULL

  Sigma3 = cast(stationPairs, s1 ~ s2, value="Cov3" )
  rownames(Sigma3) = Sigma3[,1]
  Sigma3$s1 = NULL

  Sigma = list(Sigma1, Sigma2, Sigma3)

  #Simulate errors for each time step
  t = sort(unique(noError$time))
  n = c(800, 4500-800, 6105-4500)
  sim = lapply(1:3, function(i){
    rand = rmvnorm(n[i], sigma=as.matrix(Sigma[[i]]))
    sim = data.frame( StationID=rep(colnames(Sigma[[i]]),times=n[i])
                     ,Error=as.numeric(t(rand)))
  })
  sim = do.call("rbind", sim)
  sim$Time = rep(t, times=371)
  sim$StationID = as.character( sim$StationID )
  noError$station = as.character( noError$station )
  sim = merge(sim, noError, by.x=c("StationID", "Time"), by.y=c("station", "time") )
  sim$Value = sim$model + sim$Error
  
  #Contaminate, if desired
  if(randErrorPct>0){
    load("~/Professional Files/Mines/SmartGeo/Queens/Results/global_outliers.RData")
    errorRows = sample(1:nrow(sim), size=nrow(sim)*randErrorPct)
    sim[errorRows,"Value"] = sim[errorRows,"Value"] +
      sample(largeErrors, size=length(errorRows), replace=T)
    sim$randErrorFl = (1:nrow(sim)) %in% errorRows
  }
  if(sysErrorThresh<1e5){
    load("~/Professional Files/Mines/SmartGeo/Queens/Results/changepoints.RData")
    sysErrorPct = sum(abs(shift$stat)>sysErrorThresh) / nrow(shift)
    sim = ddply(sim, "StationID", function(df){
      diffT = as.numeric( difftime(df$Time[-1], df$Time[-nrow(df)], units="days") )
      errorFl = rbinom(sum(diffT>3), size=1, p=sysErrorPct)
      errorRows = (1:(nrow(df)-1))[diffT>3][errorFl==1]
      df$sysError = 0
      for(i in errorRows){
        err = sample(shift$diff[abs(shift$stat)>sysErrorThresh], size=1)
        df$Value[(i+1):nrow(df)] = df$Value[(i+1):nrow(df)] + err
        df$sysError[i] = err
      }
      return(df)
    })
    return(sim)
  }
  
  return(sim)
}