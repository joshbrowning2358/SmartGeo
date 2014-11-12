library(ggplot2)
library(gridExtra)
library(plyr)
library(reshape)
library(dplyr)
library(sp)
if(!grepl("ch120", Sys.info()[4])){
  library(spacetime)
  library(gstat)
}
if(Sys.info()[1]=="Linux"){
  library(CompRandFld)
}
#library(spdep)
#library(fields)
#library(RandomFields)
library(sqldf)
library(MASS)
library(snht)
library(spam)
library(mvoutlier)

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
#d: Data.frame with columns s, t, and value for location, time, and observed data, respectively.
#s: Matrix.  First column is ID, subsequent columns are coordinates (typically 2 col=(x,y))
#maxs: maximum distance to use between locations as a proportion of maximum difference.
#maxt: maximum time lag to model, as a proportion of maximum difference.
#s_brks: vector of inceasing order partitioning distances.  Output values will be
#  grouped in buckets of (0,s_brks[1]), (s_brks[1],s_brks[2]), etc.  NULL means all
#  pairs will be returned.
#t_brks: vector of increasing order partitioning time, similar to s_brks.
#angle_brks: vector of angles to partition the domain.  Similar to s_brks, and useful
#  for checking anisotropy.  Units in degrees.
#Output:
#data.frame with all the squared differences between observations, or if any of s_brks,
#  t_brks, or angle_brks are defined then aggregated data.
variogramPts = function(d, s, maxs=0.5, maxt=0.25
      ,s_brks=NULL, t_brks=NULL, angle_brks=NULL){
  #d
  stopifnot(c("s","t","value") %in% colnames(d))
  stopifnot(is(d$t,"numeric"))
  d = d[!is.na(d$value),]
  #s
  stopifnot(colnames(s)[1]=="s")
  if(!is(s[,1],"factor") & !is(s[,1],"numeric"))
    stop("s's first column must be a factor or numeric!")
  #maxs
  stopifnot( maxs>0 & maxs<=1 )
  #maxt
  stopifnot( maxt>0 & maxt<=1 )
  #s_brks
  if(!is.null(s_brks)){
    stopifnot(is(s_brks,"numeric"))
    stopifnot(all(diff(s_brks)>0))
    stopifnot(s_brks[1]>=0)
    if(s_brks[1]!=0)
      s_brks = c(0, s_brks)
  }
  #t_brks
  if(!is.null(t_brks)){
    stopifnot(is(t_brks,"numeric"))
    stopifnot(all(diff(t_brks)>0))
    stopifnot(t_brks[1]>=0)
    if(t_brks[1]!=0)
      t_brks = c(0, t_brks)
  }
  #angle_brks
  if(!is.null(angle_brks)){
    stopifnot(is(angle_brks,"numeric"))
    stopifnot(all(diff(angle_brks)>0))
    stopifnot(angle_brks[1]>=0)
    stopifnot(all(angle_brks<=180))
    if(angle_brks[1]!=0)
      angle_brks = c(0, angle_brks)
    if(angle_brks[length(angle_brks)]!=180)
      angle_brks = c(angle_brks, 180)
  }

  #Manipulate distance data.frame
  dist = as.matrix(t(combn(s[,1],2)))
  dist = data.frame(s1=dist[,1], s2=dist[,2])
  dist = merge(dist, s, by.x="s1", by.y="s")
  dist = merge(dist, s, by.x="s2", by.y="s", suffixes=c("1","2"))
  dist$dist = sqrt( (dist$Easting1-dist$Easting2)^2 + (dist$Northing1-dist$Northing2)^2)
  dist$theta = atan( (dist$Northing1-dist$Northing2) / (dist$Easting1-dist$Easting2))*180/pi
  dist$theta[dist$theta<0] = dist$theta[dist$theta<0] + 180
  dist = dist[,c("s1", "s2", "dist", "theta")]
  #Add pairs like (s_i, s_i) since they won't show up with dist
  dist = rbind(dist, data.frame(s1=s[,1], s2=s[,1], dist=0, theta=0))
  
  maxs = max(dist$dist)*maxs
  maxt = (max(d$t)-min(d$t))*maxt
  
  #Merge tables using sqldf
#   d2 = d
#   temp = sqldf(paste("
#     select
#        d.s        as s1
#       ,d.t        as t1
#       ,d.value    as val1
#       ,d2.s       as s2
#       ,d2.t       as t2
#       ,d2.value   as val2
#       ,dist.dist
#     from
#       d join dist --No criteria implies cartesian product 
#         join d2 on d2.s = dist.s2
#     where
#           d.s         = dist.s1
#       and dist.dist  <= ",maxs,"
#       and d2.t       <= d.t + ",maxt,"
#       and d2.t       >= d.t
#   "))
  
  #Merge tables
  out = merge(d[,c("s", "t", "value")], dist, by.x="s", by.y="s1")
  #Filter out pairs beyond maxs
  out = out[out$dist<=maxs,]
  out = merge(out, d[,c("s","t","value")], by.x=c("s2"), by.y=c("s"))
  #Filter out pairs beyond maxt
  out = out[(out$t.x>=out$t.y) & out$t.x-out$t.y<=maxt,]
  out = data.frame( EZ2 = (out$value.y-out$value.x)^2
                   ,delta_t = out$t.x - out$t.y
                   ,delta_s = out$dist
                   ,theta = out$theta )
  #Remove pairs at the same point and time
  out = out[out$delta_t>0 | out$delta_s>0,]

  if(!is.null(s_brks))
    out$delta_s = sapply(out$delta_s, function(x) s_brks[sum(s_brks<x)+1] )
  if(!is.null(t_brks))
    out$delta_t = sapply(out$delta_t, function(x) t_brks[sum(t_brks<x)+1] )
  if(!is.null(angle_brks))
    out$theta = sapply(out$theta, function(x) angle_brks[sum(angle_brks<x)+1] )
  if(!is.null(s_brks) | !is.null(t_brks) | !is.null(angle_brks))
    out = ddply(out, c("delta_t", "delta_s", "theta"), function(df) mean(df$EZ2, na.rm=T) )

  return(out)
}

#Input:
#data: the spatio-temporal data to model.  Must have at least three columns:
#  StationID: The identifier for the location, to tie with sp.
#  Time: The time of the observation
#  <varName>: The character string supplied in the arguments of this function
#Currently, all other columns are ignored.
#Note: If NULL, the Queens' object "ground" is used.
#sp: Define a SpatialPoints object for modeling.  If NULL, defaults to stations.  If
# not NULL, it should be a SpatialPoints object with a column called "StationID"
#time: vector of times for filtering data.  If NULL, no filtering is done
#varName: Name of the dependent variable being modeled.  Must be a column of data.
#mod: A variogram model to fit to the data, as created by vgmST.  If NULL, a separable
#  exponential model is used.  Ignored if !fitModel
#Output:
#Spatio-temporal variogram
empVario = function(data=loadGround(100), sp=loadSp(), time=sort(unique(data$Time))
    ,varName="Value", mod=NULL, ...){
  #data
  stopifnot(is(data, "data.frame"))
  stopifnot(all(c("StationID", "Time", varName) %in% colnames(data)))
  #sp
  stopifnot("StationID" %in% names(sp))
  stopifnot(is(sp,"SpatialPoints"))
  #time
  stopifnot(all(time %in% data$Time))
  #varName
  stopifnot(is(varName, "character"))
  
  data = filter(data, StationID %in% sp$StationID, Time %in% time)
  if(nrow(unique( data[,c("StationID", "Time")]))!=nrow(data))
    stop("data has multiple observations for one StationID-Time pair!")
#  #Remove time points with no data.  NOTE: This removes obs, but will error out later
#   filt = ddply(ground, "Time", function(df){
#     sum(!is.na(df[,varName]))
#   })
#   data = filter(data, Time %in% filt$Time[filt$V1>0])
  
  #gStat fitting:
  d = expand.grid(StationID=sp$StationID, Time=time)
  d = merge(d, data, by=c("StationID", "Time"), all.x=TRUE)
  d = arrange(d, Time, StationID)
  stdf = STFDF(sp=sp, time=time, data=d[,varName,drop=F] )
  vst = variogramST( formula=as.formula( paste0(varName,"~ 1") ), data=stdf, ... )
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

   vstModel = fit.StVariogram(vst, model=mod, lower=rep(0,5), method="L-BFGS-B")
   return(list(vst=vst, vstModel=vstModel))
}

#fit: An object returned from empVario()
#boundaries: Spatial boundaries passed to variogramST.  his argument allows this
#  function to correctly define spacelag.
#model: logical.  Should the model be plotted as well.  CAUTION: Assumes a spherical
#  variogram model!!!
#adj: The variogram fitted by gstat::variogram assumes alpha is defined as clockwise from
#  N, but my code assumes alpha is counterclockwise from E.  The gstat::variogram model
#  was used for tlag=0 only, so if adj=T then the alpha's for tlag=0 are adjusted to my
#  definition.
#scale: argument to be passed to facet_wrap.  If "fixed", then all variogram plots will
#  have the same scale.  Alternatives are "free", "free_x", and "free_y".
#rmAni: Should the grouping by angle be removed for the plot (and hence only one plot created)?
plot.empVario = function(fit, boundaries=NULL, model=F, adj=FALSE, scale="fixed", rmAni=F){  
  toPlot = data.frame(fit$vst)
  theta_s = fit$vstModel$space
  theta_s = c(nugget=theta_s[1,2], 1, range=theta_s[2,3])
  theta_t = fit$vstModel$time
  theta_t = c(nugget=theta_t[1,2], 1, range=theta_t[2,3])
  toPlot$fit = fit$vstModel$sill*sphericalN(theta_s, toPlot$dist)*
                                 sphericalN(theta_t, toPlot$timelag)
  
  if("dir.hor" %in% colnames(toPlot)){
    toPlot = toPlot[!is.na(toPlot$dir.hor),]
    if(adj){
      toPlot$dir.hor[toPlot$timelag==0] = 90-toPlot$dir.hor[toPlot$timelag==0]
      toPlot$dir.hor[toPlot$dir.hor<0] = toPlot$dir.hor[toPlot$dir.hor<0]+180
    }
  }
  
  if(!is.null(boundaries)){
    #Adjust boundaries so first element is 0, next element is very small, and then original
    boundaries = boundaries[boundaries>0]
    boundaries = c(0, boundaries[1]*.00001, boundaries)
    toPlot$spacelag = boundaries[findInterval( toPlot$dist, boundaries )]
  }
  
  if(rmAni & "dir.hor" %in% colnames(toPlot)){
    toPlot = ddply(toPlot, c("timelag", "spacelag"), function(df){
      data.frame(np=sum(df$np, na.rm=T)
                ,dist = sum(df$np*df$dist, na.rm=T)/sum(df$np, na.rm=T)
                ,gamma = sum(df$gamma*df$np, na.rm=T)/sum(df$np, na.rm=T)
                ,fit = sum(df$fit*df$np, na.rm=T)/sum(df$np, na.rm=T)
    ) } )
  }

  p_s = ggplot(toPlot, aes(x=dist, y=gamma, color=timelag, group=timelag)) +
    geom_line(aes(linetype="data"))
  p_t = ggplot(toPlot, aes(x=timelag, y=gamma, color=spacelag, group=spacelag)) +
    geom_line(aes(linetype="data"))
  if(model){
    p_s = p_s + geom_line(aes(y=fit, linetype="model")) +
          scale_linetype_manual(breaks=c("data", "model"), values=c(2,1))
    p_t = p_t + geom_line(aes(y=fit, linetype="model")) +
          scale_linetype_manual(breaks=c("data", "model"), values=c(2,1))
  }
  if("dir.hor" %in% colnames(toPlot) & !rmAni){
    p_s = p_s + facet_wrap( ~ dir.hor, scale=scale )
    p_t = p_t + facet_wrap( ~ dir.hor, scale=scale )
  }
  print(arrangeGrob(p_s, p_t))
}

#Variogram models to fit via MLE.  Restricted MLE is preferred to reduce the bias
#in the estimate of Sigma.  But, since we have a ton of data, this is probably not
#a big issue.  The main thing we should check it the assumption of normality.
# #Commented as it may cause problems with gam()
# exponential=function(theta,h){
#   stopifnot(length(theta)==2)
#   theta[1]*(1-exp(-h/theta[2]))
# }

exponentialN=function(theta,h){
  stopifnot(length(theta)==3)
  theta[1] + theta[2]*(1-exp(-h/theta[3]))
}

# #Commented as it may cause problems with gam()
# gaussian=function(theta,h){
#   stopifnot(length(theta)==2)
#   theta[1]*(1-exp(-h^2/theta[2]^2))
# }

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
  ifelse(h<theta[3], theta[1] + (theta[2]-theta[1])*(3*h/(2*theta[3])-(h/theta[3])^3/2)
        ,theta[2])
}

#theta_s: Nugget psill (i.e. % of sill) and range
#theta_t: Nugget psill (i.e. % of sill) and range
#Model:
#sill*(nugget_s+(1-nugget_s)*spherical(s))*(nugget_t+(1-nugget_t)*spherical(t))
sphN_sphN = function(theta_s, theta_t, sill, h, u){
  stopifnot(length(theta_s)==2)
  stopifnot(length(theta_t)==2)
  theta_s = c(theta_s[1], 1, theta_s[2]) #Force sill to be 1 for each component
  theta_t = c(theta_t[1], 1, theta_t[2]) #Force sill to be 1 for each component
  sill*sphericalN(theta_s,h)*sphericalN(theta_t,u)
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
WRSS_st=function(theta_s, theta_t, sill, emp, modelFunc=sphN_sphN){
  #theta
  #modelFunc
  stopifnot(is(modelFunc,"function"))
  #emp
  stopifnot(all(c("np","dist","gamma", "timelag") %in% colnames(emp)))
  
  gam.theta=modelFunc(theta_s=theta_s, theta_t=theta_t, sill=sill, h=emp$dist, u=emp$timelag)
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
    ,initial_s=c(1/3, max(emp$dist, na.rm=T)/3), lower_s=c(0,0), upper_s=c(0.999,Inf)
    ,initial_t=c(1/3, max(emp$timelag, na.rm=T)/3), lower_t=c(0,0), upper_t=c(0.999,Inf)
    ,initial_sill = max(emp$gamma, na.rm=T), lower_sill=0, upper_sill=Inf){
  stopifnot(is(mod,"function"))
  stopifnot(c("np", "dist", "gamma", "timelag") %in% colnames(emp))
  
  n_s = length(initial_s)
  n_t = length(initial_t)
  WRSS_wrap = function(theta){
    WRSS_st(theta[1:n_s], theta[(n_s+1):(n_s+n_t)], theta[n_s+n_t+1], modelFunc=mod, emp=emp)
  }
  opt = optim( c(initial_s, initial_t, initial_sill), WRSS_wrap, method="L-BFGS-B"
    ,lower=c(lower_s, lower_t, lower_sill), upper=c(upper_s, upper_t, upper_sill))$par
  return(list(theta_s=opt[1:n_s], theta_t=opt[(n_s+1):(n_s+n_t)], sill=opt[n_s+n_t+1]))
}

#Takes an object from fitModelST (assuming sphN*sphN) and computes the covariance.
gamma2cov = function(vstModel){
# gamma(h, u) = Var( Z(s0, t0) - Z(s0+h, t0+u) )
#  = Var(Z(s0, t0)) + Var(Z(s0+h, t0+u)) - 2*Cov(Z(s0,t0), Z(s0+h, t0+u))
#  = 2 (sill/2) - 2*Cov(Z(s0,t0), Z(s0+h, t0+u))
# 2*Cov(h,u) = sill - gamma(h,u)
# Cov(h,u) = sill/2 - gamma(h,u)/2  
  theta_s = vstModel$space
  theta_s = c(nugget=theta_s[1,2], 1, range=theta_s[2,3])
  theta_t = vstModel$time
  theta_t = c(nugget=theta_t[1,2], 1, range=theta_t[2,3])
  
  function(h, u){
    vstModel$sill*( 1 - sphericalN(theta_s, h)*sphericalN(theta_t, u))/2
  }
}

# load("Results/new_sp_variograms.RData")
# load("Data/station_pairs.RData")
# fit = fits[[1]]
# model = fit[[2]]
# #Adjust the space and time ranges down for faster simulation during testing
# model$space$range = c(0, 10)
# model$time$range = c(0, 10)

#Input:
#model: output from fitModelST, i.e. a list of three elements:
#  space: data.frame with columns model, psill, and range and rows corresponding to Nug and Sph
#  time: same as space, but with data for temporal variogram
#  sill: global sill for the covariance model
#station_pairs: Should be a pairwise matrix with columns s1 (station 1), s2 (station 2)
#  dist (distance between s1 and s2) and theta (angle between s1 and s2).  All pairs should
#  be included (i.e. if n stations, should be n^2 entries)
#t: vector of times to simulate on
#Output:
#data.frame with three columns:
#s_id: ID for the spatial location, corresponds to the row of sp
#t_id: ID for the temporal value, corresponds to the element of t
#value: realization of the random variable at the corresponding point/time
simulate = function(model, station_pairs, t=0:10*24, method=2){  
  #model
  stopifnot(is(model,"StVariogramModel"))
  #station_pairs
  stopifnot(all(colnames(station_pairs)==c("s1","s2","dist","theta")))
  #t
  stopifnot(is(t,"numeric"))
  
  sRange = model$space$range[model$space$model=="Sph"]
  tRange = model$time$range[model$time$model=="Sph"]
  tPairs = merge( data.frame(t), data.frame(t), by=NULL )
  tDiffs = abs(tPairs[,1]-tPairs[,2])
  zeros = sum(tDiffs>tRange)*sum(station_pairs$dist>sRange)
  cat("Non-zero entries:\t", nrow(station_pairs)*length(t)^2-zeros,"\n")
  cat("Zero entries:\t\t\t\t\t", zeros,"\n")
  if(interactive()){
    readline("Would you like to continue? Press ENTER (or escape out) ")
  }

  covMod = gamma2cov(model)

  if(method==1){
    #Relabel the station IDs to be 1, 2, 3, ...
    loc = data.frame(name=union(station_pairs$s1, station_pairs$s2))
    loc$id = 1:nrow(loc)
    station_pairs = merge( station_pairs, loc, by.x="s1", by.y="name")
    station_pairs$s1 = station_pairs$id
    station_pairs$id = NULL
    station_pairs = merge( station_pairs, loc, by.x="s2", by.y="name")
    station_pairs$s2 = station_pairs$id
    station_pairs$id = NULL
    
    covPairs = ddply(station_pairs, "s1", function(df){
      underSRange = df[df$dist<=sRange,]
      overSRange = df[df$dist>sRange,]
      out = NULL
      s1 = df$s1[1] #All are the same
      for(t1 in 2:length(t)){
        toBind = data.frame(s1=s1, t1=t[t1], underSRange[,c("s2","dist")])
        toBind = merge(toBind, data.frame(t2=t[1:(t1-1)]))
        toBind$cov = covMod(toBind$dist, toBind$t1-toBind$t2)
        out = rbind(out, toBind)
      }
      for(t1 in 2:length(t)){
        toBind = data.frame(s1=s1, t1=t[t1], overSRange[,c("s2","dist")])
        t2 = data.frame(t2=t[1:(t1-1)])
        #Can remove obs. that are outside of the tRange
        t2 = t2[t2$t2<=tRange,,drop=F]
        toBind = merge(toBind, t2)
        toBind$cov = covMod(toBind$dist, toBind$t1-toBind$t2)
        out = rbind(out, toBind)
      }
      #Add on all pairs that occur at the same time
      for(t1 in 1:length(t)){
        toBind = data.frame(s1=s1, t1=t[t1], s2=df$s2, dist=df$dist, t2=t[t1])
        toBind$cov = covMod(toBind$dist, toBind$t1-toBind$t2)
        out = rbind(out, toBind)
      }
      
      return(out)
    } )
    
    index = data.frame(station=union(station_pairs$s1, station_pairs$s2))
    index = merge(index, data.frame(t=t) )
    index$ID = 1:nrow(index)
    
    covPairs = merge(covPairs, index, by.x=c("s1", "t1"), by.y=c("station", "t") )
    colnames(covPairs)[colnames(covPairs)=="ID"] = "rowID"
    covPairs = merge(covPairs, index, by.x=c("s2", "t2"), by.y=c("station", "t") )
    colnames(covPairs)[colnames(covPairs)=="ID"] = "colID"
    #Add elements below the diagonal, but filter out diagonal entries
    covPairs = covPairs[,c("rowID", "colID", "cov")]
    covPairs2 = covPairs[covPairs$rowID!=covPairs$colID,c("colID", "rowID", "cov")]
    colnames(covPairs2) = colnames(covPairs)
    covPairs = rbind(covPairs, covPairs2)
  } else {
    id = data.frame(s=union(station_pairs[,1], station_pairs[,2]))
    id = merge(id, data.frame(t=t))
    id$ID = 1:nrow(id)
    
    covPairs1 = station_pairs[station_pairs$dist<=sRange,]
    covPairs2 = station_pairs[station_pairs$dist>sRange,]
    timePairs = data.frame(t1=t)
    timePairs = merge(timePairs, data.frame(t2=t))
    covPairs1 = merge(covPairs1, timePairs, by=NULL)
    covPairs2 = merge(covPairs2, timePairs[abs(timePairs$t1-timePairs$t2)<=tRange,], by=NULL)
    covPairs = rbind(covPairs1, covPairs2)
    covPairs$cov = covMod(covPairs$dist, abs(covPairs$t1-covPairs$t2))

    #Add rowID, colID
    covPairs = merge(covPairs, id, by.x=c("s1", "t1"), by.y=c("s", "t"))
    colnames(covPairs)[colnames(covPairs)=="ID"] = "rowID"
    covPairs = merge(covPairs, id, by.x=c("s2", "t2"), by.y=c("s", "t"))
    colnames(covPairs)[colnames(covPairs)=="ID"] = "colID"
  }
  
  Sigma = as.spam(list(i=covPairs$rowID, j=covPairs$colID, covPairs$cov))
  out = rmvnorm.spam(1, Sigma=Sigma)
  return(out)
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
# pairwiseSNHT = function(d){
#   
# }
#NO NEED!  Implemented in your snht package.

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a random error
#algorithm.
#...?
#Output:
#Homogenized dataset, where homogenization is done via the isolation metric defined in
#Filzmoser's paper: Identification of local multivariate outliers (2014).
#Note: this method may also detect random errors.
isolation = function(d){
  library(mvoutlier)
  dat = data.frame(dat)
  dat$X = X
  dat$Y = Y
  ggplot(dat, aes(x=X1, y=X2, color=X) ) + geom_point()
  ggplot(dat, aes(x=X1, y=X2, color=Y) ) + geom_point()
  temp = locoutNeighbor(as.matrix(dat[,1:2]),X,Y,propneighb=0.1,chisqqu=0.975,
    variant="knn",usemax=1,npoints=100,indices=c(1,11,24,36))
}

#Input:
#d: contaminated dataset from contaminate(), possibly after passed to a homogenization
#algorithm.
#...?
#Output:
#Identification of random errors (probably a vector of scores)
outlier = function(d){
  
}









####################################################################################
# Edited functions from gstat (to allow anisotropic space-time variograms)
####################################################################################

variogramST <- function (formula, locations, data, ..., tlags = 0:15, cutoff, 
    width = cutoff/15, boundaries = seq(0, cutoff, width), progress = interactive(), 
    pseudo = TRUE, assumeRegular = FALSE, na.omit = FALSE) 
{
    if (missing(data)) 
        data = locations
    if (missing(cutoff)) {
        ll = !is.na(is.projected(data@sp)) && !is.projected(data@sp)
        cutoff <- spDists(t(data@sp@bbox), longlat = ll)[1, 
            2]/3
    }
    if (is(data, "STIDF")) 
        return(variogramST.STIDF(formula, data, tlags, cutoff, 
            width, boundaries, progress, ...))
    stopifnot(is(data, "STFDF") || is(data, "STSDF"))
    it = index(data@time)
    if (assumeRegular || zoo:::is.regular(zoo:::zoo(matrix(1:length(it)), 
        order.by = it), strict = TRUE)) {
        twidth = diff(it)[1]
        tlags = tlags[tlags <= min(max(tlags), length(unique(it)) - 
            1)]
    }
    else {
        warning("strictly irregular time steps were assumed to be regular")
        twidth = mean(diff(it))
    }
    ret = vector("list", length(tlags))
    obj = NULL
    t = twidth * tlags
    if (progress) 
        pb = txtProgressBar(style = 3, max = length(tlags))
    for (dt in seq(along = tlags)) {
        ret[[dt]] = StVgmLag(formula, data, tlags[dt], pseudo = pseudo, 
            boundaries = boundaries, ...)
#         ret[[dt]] = StVgmLag(formula, data, tlags[dt], pseudo = pseudo, 
#             boundaries = boundaries, alpha=c(0,45,90,135))
        ret[[dt]]$id = paste("lag", dt - 1, sep = "")
        if (progress) 
            setTxtProgressBar(pb, dt)
    }
    if (progress) 
        close(pb)
    v = do.call(rbind, ret)
    v$timelag = rep(t, sapply(ret, nrow))
    if (is(t, "yearmon")) 
        class(v$timelag) = "yearmon"
    b = attr(ret[[2]], "boundaries")
    b = c(0, b[2]/1000000, b[-1])
    b = b[-2]
    v$spacelag = c(0, b[-length(b)] + diff(b)/2)
    class(v) = c("StVariogram", "data.frame")
    if (na.omit) 
        v <- na.omit(v)
    attr(v$timelag, "units") <- attr(twidth, "units")
    if (isTRUE(!is.projected(data))) 
        attr(v$spacelag, "units") = "km"
    return(v)
}

#edit(gstat:::StVgmLag)
StVgmLag <- function (formula, data, dt, pseudo, boundaries, ...) 
{
    dotLst <- list(...)
    .ValidObs = function(formula, data) !is.na(data[[as.character(as.list(formula)[[2]])]])
    d = dim(data)
    ret = vector("list", d[2] - dt)
    if (dt == 0) {
        for (i in 1:d[2]) {
            d0 = data[, i]
            valid = .ValidObs(formula, d0)
            if (sum(valid) <= 1) 
                ret[[i]] <- NULL
            else {
                d0 = d0[valid, ]
                ret[[i]] = variogram(formula, d0, boundaries=boundaries, ...)
            }
        }
    }
    else {
        for (i in 1:(d[2] - dt)) {
            d1 = data[, i]
            valid1 = .ValidObs(formula, d1)
            d2 = data[, i + dt]
            valid2 = .ValidObs(formula, d2)
            if (sum(valid1) == 0 || sum(valid2) == 0) 
                ret[[i]] <- NULL
            else {
                d1 = d1[valid1, ]
                d2 = d2[valid2, ]
#                 obj = gstat(NULL, paste("D", i, sep = ""), formula, 
#                   d1, set = list(zero_dist = 3), beta = 0)
#                 obj = gstat(obj, paste("D", i + dt, sep = ""), 
#                   formula, d2, beta = 0)
#                 ret[[i]] = variogram(obj, cross = "ONLY", pseudo = pseudo,
#                   ...)
                alpha = 0
                if("alpha" %in% names(dotLst)) alpha=dotLst$alpha
                  ret[[i]] = crossVariogram(d1=d1, d2=d2, boundaries=boundaries, alpha=alpha)
            }
        }
    }
    VgmAverage(ret, boundaries, alphaFl="alpha" %in% names(dotLst), ...)
}

#I edited this so it keeps alpha separated
VgmAverage <- function (ret, boundaries, alphaFl, ...)
{
    dotLst <- list(...)
    alpha = 0
    if("alpha" %in% names(dotLst))
      alpha = dotLst$alpha
    ret = ret[!sapply(ret, is.null)]
    #Need access to the new VgmFillNA function (NOT the default in gstat)
    ret = lapply(ret, VgmFillNA, boundaries = c(0, 0.000001 * 
        boundaries[2], boundaries[-1]), alpha=alpha)
    np = apply(do.call(cbind, lapply(ret, function(x) x$np)), 
        1, sum, na.rm = TRUE)
    gamma = apply(do.call(cbind, lapply(ret, function(x) x$gamma * 
        x$np)), 1, sum, na.rm = TRUE)/np
    dist = apply(do.call(cbind, lapply(ret, function(x) x$dist * 
        x$np)), 1, sum, na.rm = TRUE)/np
    #Return different v based on if alpha was supplied
    if(!alphaFl){
      v = data.frame(np = np, dist = dist, gamma = gamma)
    } else {
      #Want NA's removed from dir.hor.  Using mean, max, min will all work, but
      #mean should give weird values if something weird happens (ie. rows misalign).
      #So, this is a helpful safety check.
      dir.hor = apply(do.call(cbind, lapply(ret, function(x) x$dir.hor)),
          1, mean, na.rm = TRUE)
      v = data.frame(np = np, dist = dist, gamma = gamma, dir.hor = dir.hor)
    }
    class(v) = class(ret[[1]])
    attr(v, "boundaries") = attr(ret[[1]], "boundaries")
    v[is.na(v)] = NA
    v
}

#Original function had errors if alpha was supplied to variogramST (original call)
VgmFillNA <- function (x, boundaries, alpha)
{
  if(!"dir.hor" %in% colnames(x))
    x$dir.hor = 0
  if( !all(alpha %in% x$dir.hor) )
    x = rbind.fill(x, data.frame(dir.hor=alpha[!alpha %in% x$dir.hor]))
  x = ddply(x, "dir.hor", function(df){
    n = length(boundaries) - 1
    ix = rep(NA, n)
    ix[findInterval(df$dist, boundaries)] = 1:nrow(df)
    df[ix, ]
  })
  return(x)
}

#Helper function used in StVgmLag, as the variogram() call seems problematic
crossVariogram = function(d1, d2, boundaries, alpha=c(0)){
  d1 = data.frame(d1)
  colnames(d1) = c("z", "x", "y")
  d2 = data.frame(d2)
  colnames(d2) = c("z", "x", "y")

  #If no rows in one data.frame, don't do anything
  if(min(nrow(d1),nrow(d2))==0){
    warning("No rows in one of the desired cross-variograms.  Do you have data at all time pts?")
    return(data.frame(np=NA, dist=NA, gamma=NA, dir.hor=NA))
  }
  
  boundaries = boundaries[boundaries>0]
  #Create a bucket to split 0's from the smallest interval
  boundaries = c(0, boundaries[1]*.000001, boundaries)
  pairs = merge(data.frame(d1), data.frame(d2), by=NULL)
  pairs$dist = sqrt( (pairs$x.x-pairs$x.y)^2 + (pairs$y.x-pairs$y.y)^2 )
  pairs$dist_bucket = findInterval(pairs$dist, boundaries)
  
  pairs$theta = atan((pairs$y.x-pairs$y.y) / (pairs$x.x-pairs$x.y) )*180/pi
  pairs$theta[is.na(pairs$theta)] = 0 #theta between adjacent points
  pairs$theta[pairs$theta<0] = pairs$theta[pairs$theta<0]+180  

  midPts = c(alpha[length(alpha)]-180,alpha,alpha[1]+180)
  midPts = (midPts[2:length(midPts)-1]+midPts[2:length(midPts)])/2
  pairs$theta_bucket = findInterval(pairs$theta, midPts)
  
  #Add points paired with themselves to each angle
  if(sum(pairs$dist==0)>0){
    selfPairs = pairs[pairs$dist==0,]
    pairs = pairs[pairs$dist>0,]
    for(angle in 1:length(alpha)){
      selfPairs$theta_bucket = angle
      pairs = rbind(pairs, selfPairs)
    }
  }

  #Group edges together:
  pairs$theta_bucket[pairs$theta_bucket>length(alpha)] = 1
  pairs$theta_bucket[pairs$theta_bucket==0] = length(alpha)
  
  vg = group_by(pairs, dist_bucket, theta_bucket )
  out = summarize(vg
          ,np=n()
          ,dist=mean(dist)
          ,gamma=mean((z.x-z.y)^2)/2)
  out$dir.hor = alpha[out$theta_bucket]
  #Add in any missing theta, dist combinations so all returned objects have same groups
  missing = expand.grid(dist_bucket=1:length(boundaries), theta_bucket=1:length(alpha)) 
  missing = missing[!paste(missing$dist_bucket, missing$theta_bucket) %in%
                     paste(out$dist_bucket, out$theta_bucket),]
  out = rbind(out, data.frame(missing, np=0
                             ,dist=boundaries[missing$dist_bucket]
                             ,gamma = NA
                             ,dir.hor=alpha[missing$theta_bucket]) )
  out$theta_bucket = NULL
  out$dist_bucket = NULL
  
  out = out[order(out$dir.hor, out$dist),]
  return(out)
}

locoutNeighbor <- function (dat, X, Y, propneighb = 0.1, variant = c("dist", "knn"), 
    usemax = 1/3, npoints = 50, chisqqu = 0.975, indices = NULL, 
    xlab = NULL, ylab = NULL, colall = gray(0.7), colsel = 1, 
    ...) 
{
    if (is.null(ylab)) {
        ylab <- paste("Degree of isolation from ", (1 - propneighb) * 
            100, "% of the neighbors", sep = "")
    }
    n <- nrow(dat)
    p <- ncol(dat)
    covr <- robustbase:::covMcd(dat)
    cinv <- solve(covr$cov)
    MDglobal <- sqrt(mahalanobis(dat, covr$center, cinv, inverted = TRUE))
    qchi <- sqrt(qchisq(chisqqu, p))
    MDglobalTF <- (MDglobal < qchi)
    if (!is.null(indices)) {
        indices.reg = indices[MDglobalTF[indices]]
        indices.out = indices[!MDglobalTF[indices]]
    }
    idx <- matrix(1:n, n, n)
    sel <- as.vector(idx[lower.tri(idx)])
    hlp <- as.matrix(dat[rep(1:(n - 1), seq((n - 1), 1)), ] - 
        dat[sel, ])
    MDij <- sqrt(rowSums((hlp %*% cinv) * hlp))
    MDpair <- matrix(0, n, n)
    MDpair[lower.tri(MDpair)] <- MDij
    MDpair <- t(MDpair)
    MDpair[lower.tri(MDpair)] <- MDij
    Xmat <- matrix(rep(t(X), length(X)), ncol = dim(t(X))[2], 
        byrow = TRUE)
    Ymat <- matrix(rep(t(Y), length(Y)), ncol = dim(t(Y))[2], 
        byrow = TRUE)
    EuclD <- sqrt((c(rep(X, length(X))) - Xmat)^2 + (c(rep(Y, 
        length(Y))) - Ymat)^2)
    sortD <- apply(EuclD, 1, sort, index.return = TRUE)
    neighbmatrix <- matrix(unlist(unlist(sortD, recursive = FALSE)[seq(from = 2, 
        to = 2 * n, by = 2)]), ncol = n, nrow = n, byrow = TRUE)
    if (variant == "dist") {
        maxD <- max(EuclD) * usemax
        neigbound <- matrix(NA, nrow = n, ncol = npoints)
        vec <- seq(from = 0, to = maxD, length = npoints)
        for (i in 1:n) {
            for (j in 1:npoints) {
                MDneig <- sort(MDpair[i, EuclD[i, ] <= vec[j]])
                if (length(MDneig) > 1) {
                  MDneig <- MDneig[-1]
                }
                MDbound <- MDneig[ceiling(length(MDneig) * propneighb)]
                neigbound[i, j] <- pchisq(MDbound^2, p, MDglobal[i]^2)
            }
        }
    }
    else {
        npoints <- min(n, npoints)
        maxn <- trunc(n * usemax)
        neigbound <- matrix(NA, nrow = n, ncol = npoints)
        vec <- ceiling(seq(from = 1, to = maxn, length = npoints))
        for (i in 1:n) {
            for (j in 1:npoints) {
                MDneig <- sort(MDpair[i, neighbmatrix[i, 1:vec[j]]])
                if (length(MDneig) > 1) {
                  MDneig <- MDneig[-1]
                }
                MDbound <- MDneig[ceiling(length(MDneig) * propneighb)]
                neigbound[i, j] <- pchisq(MDbound^2, p, MDglobal[i]^2)
            }
        }
    }
    out = data.frame( MDglobal, neigbound )
    colnames(out) = c("MDglobal", paste0("MDlocal",1:ncol(neigbound)))
    return(out)
}