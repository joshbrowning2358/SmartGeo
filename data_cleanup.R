library(reshape)
library(plyr)
library(rgdal)
library(bigmemory)

setwd("~/Professional Files/Mines/SmartGeo/Queens/")

#ground contains data for ground deformation by station_ID
ground = read.csv(file="Data/AMTS_Filtered.csv")
colnames(ground)[1] = "Time"
colnames(ground) = gsub("X", "", colnames(ground))
ground[ground==-999] = NA
ground = melt(ground, id.vars="Time")
ground$Time = as.POSIXct("01-01-1900", "%d-%m-%Y", tz="EST") + ground$Time*60*60*24
colnames(ground) = c("Time", "StationID", "Value")

ground$Date = as.Date(ground$Time)
groundSmall = ddply(ground, c("Date", "StationID"), function(df){
  out = apply(df, 2, mean, na.rm=T)
  return(out)
} )

# setwd("Data")
# groundBM = big.matrix(nrow=nrow(ground), ncol=15, backingfile="groundBM"
#         ,descriptorfile="groundBM.desc")
# setwd("..")
# cnames = rep(NA, 15)
# ground[,2] = as.numeric( as.character( ground[,2] ) )
# for(i in 1:3){
#   groundBM[,i] = ground[,i]
#   cnames[i] = colnames(ground)[i]
# }

#station contains the easting and northing for all stations
station = read.csv(file="Data/AMTS Sensor Coordinates.csv")
station = data.frame(t(station))
colnames(station) = c("Easting", "Northing")
station$StationID = gsub("X","",rownames(station))
rownames(station) = NULL
station = cbind(station, project(as.matrix(station[,c("Easting","Northing")])
  ,proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
colnames(station)[4:5] = c("Longitude", "Latitude")


#tunnel contains the location of the tunnels
tunnel = read.csv(file="Data/AMTS Tunnel Alignment.csv", skip=1)
cnames = c("AdvanceFl", "Easting", "Northing", "Elevation")
colnames(tunnel) = c("Time", paste0(cnames, "YL"), paste0(cnames, "A"), paste0(cnames, "D"), paste0(cnames, "BC") )
tunnel$Time = as.POSIXct("01-01-1900", "%d-%m-%Y", tz="EST") + tunnel$Time*60*60*24

distA = data.frame(dist=rep(NA, nrow(ground)))
distA$Time = ground$Time
distA$StationID = ground$StationID
distYL = distA
distD = distA
distBC = distA
for( i in 1:nrow(station)){
  ID = station[i,]$StationID
  loc = station[i,c("Easting", "Northing")]
  distYL$dist[distYL$StationID==ID] = apply(tunnel[,c("EastingYL", "NorthingYL")], 1, function(x){sqrt(sum((x-loc)^2))} )
  distBC$dist[distBC$StationID==ID] = apply(tunnel[,c("EastingBC", "NorthingBC")], 1, function(x){sqrt(sum((x-loc)^2))} )
  distA$dist[distA$StationID==ID] = apply(tunnel[,c("EastingA", "NorthingA")], 1, function(x){sqrt(sum((x-loc)^2))} )
  distD$dist[distD$StationID==ID] = apply(tunnel[,c("EastingD", "NorthingD")], 1, function(x){sqrt(sum((x-loc)^2))} )
  cat("Row",i,"completed.\n")
}

ground = data.frame(ground, A=distA[,1], BC=distBC[,1], D=distD[,1], YL=distYL[,1])
save(ground, file="Data/ground_with_distance.RData")

groundBM[,4] = distA$dist; cnames[4] = "dist.A"
rm(distA)
groundBM[,4] = distBC$dist; cnames[4] = "dist.BC"
rm(distBC)
groundBM[,4] = distD$dist; cnames[4] = "dist.D"
rm(distD)
groundBM[,4] = distYL$dist; cnames[4] = "dist.YL"
rm(distYL)

write.csv(cnames, row.names=F, file="Data/groundBM_cnames.csv")
