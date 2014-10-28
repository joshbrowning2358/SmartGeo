setwd("~/Professional Files/Mines/SmartGeo/Queens/")
load("Data/station_pairs.RData")
load("Results/new_sp_variograms.RData")

#First, use a model with a very short range (hence sparse) and a subset of stations:
model = fits[[1]][[2]]
model$space$range = c(0,10)
model$time$range = c(0,10)
filter = as.numeric(as.character(station_pairs$s1)) <=210 &
  as.numeric(as.character(station_pairs$s2)) <=210
station_pairs = station_pairs[filter,]

simulate(model, station_pairs, t=0:1*24, method=2)

covMat = data.frame( station_pairs[,1:3] )
covMat = merge( covMat, data.frame(t1=t) )
covMat = merge( covMat, data.frame(t2=t) )
covMat$cov = covMod( covMat$dist, abs(covMat$t1-covMat$t2) )
covMat = cast( covMat, s1 + t1 ~ s2 + t2, value="cov" )
rownames(covMat) = paste0(covMat[,1], "_", covMat[,2])
covMat = covMat[,-1:-2]
covMat = as.matrix( covMat )
all(t(covMat)==covMat)
eigen(covMat)$values
