library(ggplot2)
library(plyr)
library(reshape)
library(dplyr)
library(sp)
library(gridExtra)
#library(spdep)
#library(gstat)
#library(spacetime)
#library(fields)
#library(RandomFields)
library(sqldf)
#library(CompRandFld)
library(mgcv)
library(MASS)
library(snht)

#Output:
#Generates a SpatialPoints object with the coordinates of the Queens' stations
loadSp = function(){
  load = F
  if(!exists("station")) load=TRUE
  if(exists("station"))
    if(any(!c("East2", "North2") %in% colnames(station))) load=TRUE
  if(load){
    if(Sys.info()[4]=="JOSH_LAPTOP")
      load("~/Professional Files/Mines/SmartGeo/Queens/Data/station_new_coords.RData")
    if(Sys.info()[4]=="jb")
      load("/media/storage/Professional Files/Mines/SmartGeo/Queens/Data/station_new_coords.RData")
    if(grepl("^ch", Sys.info()[4]))
      load("~/SmartGeo/Queens/Data/station_new_coords.RData")
  }
  sp = station
  sp = sp[!is.na(sp$Longitude),]
  sp = sp[,c("StationID", "East2", "North2")]
  colnames(sp) = c("StationID", "Easting", "Northing")
  coordinates(sp) = c("Easting", "Northing")
  return(sp)
}

#Input
#timeCnt: Integer.  Number of time observations to include in dataset.  If larger than
#  nrow(ground), then all rows are used.
#Output:
#Generates a data.frame with the ground data
loadGround = function(timeCnt=100){
  stopifnot(timeCnt>0)
  
  if(!exists("ground")){
    if(Sys.info()[4]=="JOSH_LAPTOP")
      load("~/Professional Files/Mines/SmartGeo/Queens/Data/ground_with_distance.RData")
    if(Sys.info()[4]=="jb")
      load("/media/storage/Professional Files/Mines/SmartGeo/Queens/Data/ground_with_distance.RData")
    if(grepl("^ch", Sys.info()[4]))
      load("~/SmartGeo/Queens/Data/ground_with_distance.RData")
  }
  times = unique(ground$Time)[1:timeCnt]
  data = filter(ground, Time %in% times)
  return(data)
}

# #Example:
# #Output: returns a list meant to look like an StVariogramModel object.  Useful on Linux,
# #where we can't fit variograms because of missing packages.
# vstMod = function(){
#   list(space=data.frame(model=c("Nug", "Exp")
#                              ,psill=c(0.2197915,0.7802085)
#                              ,range=c(0,609.344) )
#             ,time=data.frame( model=c("Nug", "Exp")
#                              ,psill=c(0.6188773,0.3811227)
#                              ,range=c(0,93.15087) )
#             ,sill=0.01647729
#             ,stModel="separable")
# }

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
empVario = function(data=loadGround(100), sp=loadSp(), time=unique(data$Time), mod=NULL, ...){
  #data
  stopifnot(is(data, "data.frame"))
  stopifnot(all(c("StationID", "Time", "Value") %in% colnames(data)))
  #sp
  stopifnot("StationID" %in% names(sp))
  stopifnot(is(sp,"SpatialPoints"))
  #time
  stopifnot(all(time %in% data$Time))

  data = filter(data, StationID %in% sp$StationID, Time %in% time)
  if(nrow(unique( data[,c("StationID", "Time")]))!=nrow(data))
    stop("data has multiple observations for one StationID-Time pair!")
  
  #gStat fitting:
  d = expand.grid(StationID=sp$StationID, Time=time)
  d = merge(d, data, by=c("StationID", "Time"), all.x=TRUE)
  d = arrange(d, Time, StationID)
  stdf = STFDF(sp=sp, time=time, data=data[,"deltaValue",drop=F] )
  vst = variogramST( deltaValue ~ 1, data=stdf, ... )
  
#   #CompRandFld fitting:
#   data = cast(data, Time ~ StationID, value="Value")
#   rownames(data) = data$Time
#   data$Time = NULL
#   vst = EVariogram(data=as.matrix(data)
#                 ,coordx=as.matrix(data.frame(sp)[,-1])
#                 ,coordt=as.numeric(time))
#   vst = FitComposite(data=as.matrix(data)
#                 ,coordx=as.matrix(data.frame(sp)[,-1])
#                 ,coordt=as.numeric(time)
#                 ,corrmodel="exp_exp")
#   return(vst)
  
   #Fit the variogram model
   if(is.null(mod))
     mod = vgmST("separable"
            ,space=vgm(psill=.50,"Sph", range=max(vst$dist,na.rm=T)/2, nugget=0.4),
            ,time =vgm(psill=.50,"Sph", range=max(vst$timelag,na.rm=T)/2, nugget=0.4),
            ,sill=max(vst$gamma, na.rm=T))
 
   vstModel = fit.StVariogram(vst, model=mod)
   return(list(vst=vst, vstModel=vstModel))
}

empVario(data=loadGround(800),tlags=0:30*9)

#fit: An object returned from empVario()
#CAUTION: Assumes a spherical variogram model!!!
plot.empVario = function(fit){
  p_s = ggplot(data.frame(fit$vst), aes(x=dist, y=gamma)) +
    geom_line(aes(color=timelag, group=timelag))
  p_t = ggplot(data.frame(fit$vst), aes(x=timelag, y=gamma)) +
    geom_line(aes(color=spacelag, group=spacelag))
  theta_s = fit$vstModel$space
  theta_s = c(nugget=theta_s[1,2]*fit$vstModel$sill, fit$vstModel$sill, range=theta_s[2,3])
  theta_t = fit$vstModel$time
  theta_t = c(nugget=theta_t[1,2]*fit$vstModel$sill, fit$vstModel$sill, range=theta_t[2,3])
  s_grid = seq(0,max(fit$vst$dist,na.rm=T),length.out=100)
  t_grid = seq(0,max(fit$vst$timelag,na.rm=T),length.out=100)
  sPlot = data.frame(x=s_grid, y=sphericalN(theta_s, s_grid))
  tPlot = data.frame(x=t_grid, y=sphericalN(theta_t, t_grid))
  p_s = p_s + geom_line(data=sPlot, aes(x=x, y=y ), color="red" )
  p_t = p_t + geom_line(data=tPlot, aes(x=x, y=y ), color="red" )
  print(arrangeGrob(p_s, p_t))
}

#Variogram models to fit via MLE.  Restricted MLE is preferred to reduce the bias
#in the estimate of Sigma.  But, since we have a ton of data, this is probably not
#a big issue.  The main thing we should check it the assumption of normality.
exponential=function(theta,h){
  stopifnot(length(theta)==2)
  theta[1]*(1-exp(-h/theta[2]))
}

exponentialN=function(theta,h){
  stopifnot(length(theta)==3)
  theta[1] + theta[2]*(1-exp(-h/theta[3]))
}

gaussian=function(theta,h){
  stopifnot(length(theta)==2)
  theta[1]*(1-exp(-h^2/theta[2]^2))
}

gaussianN=function(theta,h){
  stopifnot(length(theta)==3)
  theta[1] + theta[2]*(1-exp(-h^2/theta[3]^2))
}

spherical = function(theta,h){
  stopifnot(length(theta)==2)
  ifelse(h<theta[2], theta[1]*(3*h/(2*theta[2])-(h/theta[2])^3/2), theta[1])
}

sphericalN = function(theta,h){
  stopifnot(length(theta)==3)
  ifelse(h<theta[3], theta[1] + theta[2]*(3*h/(2*theta[3])-(h/theta[3])^3/2)
        ,theta[1]+theta[2])
}

sphN_sphN = function(theta_s, theta_t, h, u){
  stopifnot(length(theta_s)==3)
  stopifnot(length(theta_t)==3)
  sphericalN(theta_s,h)*sphericalN(theta_t,u)
}

#Weighted Residual Sums of Squares.  Used in conjunction with optim to fit model variograms.
#Inputs
#theta: the parameter vector
#modelFunc: A function taking two arguments: the parameter vector theta and the distance.
#emp: An empirical variogram with three columns: np, dist, and gamma.
WRSS=function(theta, modelFunc, emp){
  #theta
  #modelFunc
  stopifnot(is(modelFunc,"function"))
  #emp
  stopifnot(all(c("np","dist","gamma") %in% colnames(emp)))
  
	gam.theta=modelFunc(theta, emp$dist)
	sum((emp$np/(gam.theta^2))*((emp$gamma-gam.theta)^2))
}

#Spatio-temporal version of WRSS
#modelFunc: now, takes four arguments: theta_s, theta_t, h (spatial distance), and u (temporal distance)
#emp: Now need columns np, dist, gamma, and timelag
WRSS_st=function(theta_s, theta_t, emp, modelFunc=sphN_sphN){
  #theta
  #modelFunc
  stopifnot(is(modelFunc,"function"))
  #emp
  stopifnot(all(c("np","dist","gamma", "timelag") %in% colnames(emp)))
  
  gam.theta=modelFunc(theta_s, theta_t, h=emp$dist, u=emp$timelag)
	sum((emp$np/(gam.theta^2))*((emp$gamma-gam.theta)^2), na.rm=T)
}

#Input:
#emp: empirical spatio-temporal variogram, as returned by empVario
#mod: The semivariogram model.  Should take two arguments: theta (vector) and distance.
#initial, lower, upper: parameters passed to optim.
#Output:
#Simulated data
#Notes:
#If a space-time process has constant variance over space and time, then
#Cov(Z(s1,t1),Z(s2,t2)) = Var(Z) - gamma(s1-s2, t1-t2)
#We'll assume (at least in the first iteration) a separable spherical variogram, as it has
#compact support.
fitModel = function(emp, mod=spherical, initial=c(1,1), lower=c(0,0), upper=c(Inf,Inf)){
  stopifnot(is(mod,"function"))
  stopifnot(c("np", "dist", "gamma") %in% colnames(emp))
  
  optim( initial, WRSS, emp=emp, method="L-BFGS-B"
    ,modelFunc=mod, lower=lower, upper=upper)$par
}

#initial_s: parameters are sill and range, so pick reasonable parameters.
#initial_t: parameters are sill and range, so pick reasonable parameters.
fitModelST = function(emp, mod=sphN_sphN
    ,initial_s=c(max(emp$gamma, na.rm=T)/3, max(emp$gamma, na.rm=T), max(emp$dist, na.rm=T)/3)
    ,lower_s=c(0,0,0), upper_s=c(max(emp$gamma, na.rm=T),Inf,Inf)
    ,initial_t=c(max(emp$gamma, na.rm=T)/3, max(emp$gamma, na.rm=T), max(emp$timelag, na.rm=T)/3)
    ,lower_t=c(0,0,0), upper_t=c(max(emp$gamma, na.rm=T),Inf,Inf)){
  stopifnot(is(mod,"function"))
  stopifnot(c("np", "dist", "gamma", "timelag") %in% colnames(emp))
  
  n_s = length(initial_s)
  n_t = length(initial_t)
  WRSS_wrap = function(theta){
    WRSS_st(theta[1:n_s], theta[(n_s+1):(n_s+n_t)], modelFunc=mod, emp=emp)
  }
  opt = optim( c(initial_s,initial_t), WRSS_wrap, method="L-BFGS-B"
    ,lower=c(lower_s, lower_t), upper=c(upper_s,upper_t))$par
  return(list(theta_s=opt[1:n_s], theta_t=opt[(n_s+1):(n_s+n_t)]))
}

#Input:
#d: output from simulate()
#randPct: % of random errors
#sysPct: % of systematic errors
#...?
#Output:
#d with contaminated observations.  Additionally, output the original data
#so analysis of model performance will be easier.
contaminate = function(d, randPct, sysPct){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via SNHT applied to each individual
#station.
univariateSNHT = function(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via pairwise SNHT.
pairwiseSNHT = function(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via the isolation metric defined in
#Filzmoser's paper: Identification of local multivariate outliers (2014).
#Note: this method may also detect random errors.
isolation = function(d){
  
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a homogenization
#algorithm.
#...?
#Output:
#Identification of random errors (probably a vector of scores)
outlier = function(d){
  
}