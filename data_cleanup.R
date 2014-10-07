library(reshape)
library(MASS)
library(plyr)
library(dplyr)
#library(rgdal)

if(Sys.info()[4]=="JOSH_LAPTOP")
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
if(grepl("ch",Sys.info()[4]))
  setwd("~/SmartGeo/Queens")

#ground contains data for ground deformation by station_ID
ground = read.csv(file="Data/AMTS_Filtered.csv")
colnames(ground)[1] = "Time"
colnames(ground) = gsub("X", "", colnames(ground))
ground[ground==-999] = NA
ground = melt(ground, id.vars="Time")
ground$Time = as.POSIXct("01-01-1900", "%d-%m-%Y", tz="EST") + ground$Time*60*60*24
colnames(ground) = c("Time", "StationID", "Value")

#Create a new variable: differences in the values
ground = arrange(ground, StationID, Time) 
ground$deltaValue = c(NA,diff(ground$Value))
ground$deltaValue[diff(as.numeric(ground$StationID))!=0] = NA #Don't difference across stations

#Search for outliers
station = group_by( ground, StationID )
sds = summarize( station, mu=huber(Value)[[1]], s=huber(Value)[[2]] )
#ggplot( sds, aes(x=StationID, y=s ) ) + geom_point()
ground = merge(ground, sds)
ground$outlier = F
ground$stat = (ground$deltaValue-ground$mu)/ground$s
test1 = ground$stat >= 10
test2 = ground$stat <= -10
#Classify as an outlier if we have a high and then low difference.
ground$outlier[test1 & lead(test2,1)] = T
ground$outlier[test2 & lead(test1,1)] = T
ground$mu = NULL
ground$s = NULL
ground$stat = NULL

range = 15000:17000
qplot( range, ground$Value[range], geom="line" ) + geom_point(aes(color=ground$outlier[range]))
qplot( range, ground$deltaValue[range], color=ground$outlier[range] )

#station contains the easting and northing for all stations
station = read.csv(file="Data/AMTS Sensor Coordinates.csv")
station = data.frame(t(station))
colnames(station) = c("Easting", "Northing")
station$StationID = gsub("X","",rownames(station))
rownames(station) = NULL
#station = cbind(station, project(as.matrix(station[,c("Easting","Northing")])
#  ,proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
#colnames(station)[4:5] = c("Longitude", "Latitude")


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

ground$Date = as.Date(ground$Time)
day_station = group_by( ground, Date, StationID )
temp = summarise(day_station
  ,Value=mean(Value, na.rm=T)
  ,A=mean(A, na.rm=T)
  ,BC=mean(BC, na.rm=T)
  ,D=mean(D, na.rm=T)
  ,YL=mean(YL, na.rm=T)
)
groundSmall = data.frame(temp)
save(groundSmall, file="Data/ground_daily.RData")

