setwd("~/Professional Files/Mines/SmartGeo/Queens/")
load("Data/station_pairs.RData")
load("Results/new_sp_variograms.RData")

#First, use a model with a very short range (hence sparse) and a subset of stations:
model = fits[[1]][[2]]
model$space$range = c(0,10)
model$time$range = c(0,10)
#Remove nugget:
model$space$psill = c(0,1)
model$time$psill = c(0,1)
filter = as.numeric(as.character(station_pairs$s1)) <=210 &
  as.numeric(as.character(station_pairs$s2)) <=210
station_pairs = station_pairs[filter,]
covMod = gamma2cov(model)

simulate(model, station_pairs, t=0:1*24, method=2)

makeCovMat = function(station_pairs, t, covMod){
  covMat = data.frame( station_pairs[,1:3] )
  covMat = merge( covMat, data.frame(t1=t) )
  covMat = merge( covMat, data.frame(t2=t) )
  
  covMat$cov = covMod( covMat$dist, abs(covMat$t1-covMat$t2) )
  covMat = cast( covMat, s1 + t1 ~ s2 + t2, value="cov" )
  rownames(covMat) = paste0(covMat[,1], "_", covMat[,2])
  covMat = covMat[,-1:-2]
  covMat = as.matrix( covMat )
  return(covMat)
}

#t = 0:1*24
t = 0:0*24
load("Data/station_pairs.RData")
vals = unique(station_pairs$s1)[1:10]
station_pairs = station_pairs[station_pairs$s1 %in% vals & station_pairs$s2 %in% vals,]

Sigma.s = makeCovMat(station_pairs, t, covMod)
all(t(Sigma.s)==Sigma.s)
eigen(Sigma.s)$values


t = 0:9*24
load("Data/station_pairs.RData")
vals = unique(station_pairs$s1)[1:1]
station_pairs = station_pairs[station_pairs$s1 %in% vals & station_pairs$s2 %in% vals,]

Sigma.t = makeCovMat(station_pairs, t, covMod)
all(t(Sigma.t)==Sigma.t)
eigen(Sigma.t)$values


t = 0:2*24
load("Data/station_pairs.RData")
vals = unique(station_pairs$s1)[1:3]
station_pairs = station_pairs[station_pairs$s1 %in% vals & station_pairs$s2 %in% vals,]

Sigma.st = makeCovMat(station_pairs, t, covMod)
all(t(Sigma.st)==Sigma.st)
eigen(Sigma.st)$values
Sigma.st[1:4,1:4]


########################################################################################
# Very Simple Example
########################################################################################

#Assume range is 1 in space and time, and separable spherical model:
gamma = function(h){
  ifelse(h>1, 1, 3*h/2-h^3/2)
}
qplot( 0:12/10, gamma(0:12/10) )
s = 0:10/10
t = 0:10/10

Sigma.h = data.frame( s1=s )
Sigma.h = merge( Sigma.h, data.frame(s2=s) )
Sigma.h$dist = abs(Sigma.h$s1-Sigma.h$s2)
Sigma.h$cov = 1-gamma(Sigma.h$dist)
Sigma.h = cast( Sigma.h, s1 ~ s2, value="cov" )
Sigma.h$s1 = NULL
Sigma.h = as.matrix(Sigma.h)
det(Sigma.h); eigen(Sigma.h)$values

Sigma.t = data.frame( t1=t )
Sigma.t = merge( Sigma.t, data.frame(t2=t) )
Sigma.t$dist = abs(Sigma.t$t1-Sigma.t$t2)
Sigma.t$cov = 1-gamma(Sigma.t$dist)
Sigma.t = cast( Sigma.t, t1 ~ t2, value="cov" )
Sigma.t$t1 = NULL
Sigma.t = as.matrix(Sigma.t)
det(Sigma.t); eigen(Sigma.t)$values

Sigma.st = data.frame( s1=s )
Sigma.st = merge( Sigma.st, data.frame(t1=t) )
Sigma.st = merge( Sigma.st, data.frame(s2=s) )
Sigma.st = merge( Sigma.st, data.frame(t2=t) )
Sigma.st$h = abs(Sigma.st$s1-Sigma.st$s2)
Sigma.st$u = abs(Sigma.st$t1-Sigma.st$t2)
Sigma.st$cov = 1-gamma(Sigma.st$h)*gamma(Sigma.st$u)
Sigma.st = cast( Sigma.st, s1 + t1 ~ s2 + t2, value="cov" )
Sigma.st$s1 = NULL
Sigma.st$t1 = NULL
Sigma.st = as.matrix(Sigma.st)
det(Sigma.st); eigen(Sigma.st)$values

s = 0:2/2; t=0:2/2
Sigma.st = data.frame( s1=s )
Sigma.st = merge( Sigma.st, data.frame(t1=t) )
Sigma.st = merge( Sigma.st, data.frame(s2=s) )
Sigma.st = merge( Sigma.st, data.frame(t2=t) )
Sigma.st$h = abs(Sigma.st$s1-Sigma.st$s2)
Sigma.st$u = abs(Sigma.st$t1-Sigma.st$t2)
Sigma.st$cov = 1-gamma(Sigma.st$h)*gamma(Sigma.st$u)
Sigma.st = cast( Sigma.st, s1 + t1 ~ s2 + t2, value="cov" )
Sigma.st$s1 = NULL
Sigma.st$t1 = NULL
Sigma.st = as.matrix(Sigma.st)
det(Sigma.st); eigen(Sigma.st)$values
