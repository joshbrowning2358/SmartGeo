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
simulate = function(spatialV, nnetMod, stationPairs, covariate, station, time){
  stopifnot(length(spatialV)==3)
  stopifnot(all(lapply(spatialV, length)==2))
  stopifnot(is(nnetMod,"nnet"))
  stopifnot(all(colnames(stationPairs)==c("s1", "s2", "dist", "theta")))
  #Since we cast with s1 ~ s2 and want a symmetric matrix, we need to ensure s1 is of
  #the same type as s2
  stopifnot(is(stationPairs$s1)==is(stationPairs$s2))
  
  #Estimate the expected values from nnet
  noError = predict(nnetMod, newdata=covariate)
  noError = data.frame(station=station, time=time, noError=noError[,1])
  
  #Create the three covariance matrices
  stationPairs$Cov1 = spatialV[[1]][[2]][1] - modSinExp(spatialV[[1]][[2]], stationPairs$dist, 0)
  stationPairs$Cov2 = spatialV[[2]][[2]][1] - modSinExp(spatialV[[2]][[2]], stationPairs$dist, 0) 
  stationPairs$Cov3 = spatialV[[3]][[2]][1] - modSinExp(spatialV[[3]][[2]], stationPairs$dist, 0)

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

  #Simulate errors for each time step
  t = sort(unique(noError$time))
  for(i in 1:800){
    rnorm
  }
}

load("Data/ground_for_nnet.RData")
load("Data/station_pairs.RData")
load("Results/nnet_model_all_data_new.RData")
load("Results/sin_sp_variograms.RData")
station = station[!is.na(station$East2),]

plot.empVario(fits[[1]], adj=T, rmAni=T, boundaries=0:65*10, model="sinexp_exp")
plot.empVario(fits[[2]], adj=T, rmAni=T, boundaries=0:65*10, model="sinexp_exp")
plot.empVario(fits[[3]], adj=T, rmAni=T, boundaries=0:65*10, model="exponential")

nnetMod$coefnames
head(ground,2)
covariate = ground[,c("nnetA", "nnetBC", "nnetD", "nnetYL", "nnetEast2", "nnetNorth2", "nnetTime")]
colnames(covariate) = c("A", "BC", "D", "YL", "East2", "North2", "Time")

# sim = RFsim(coordx=as.matrix(station[1:371,c("East2","North2")]), coordt=0:80*24
#       ,corrmodel="exp_exp", grid=FALSE
#       ,param=list(nugget=0.3, mean=0, scale_s=988.5665, scale_t=4.721552, sill=0.03278) )