library(ggplot2)
library(ggmap)
library(gridExtra)
library(scales)
library(plyr)
library(reshape)
library(spdep)
library(gstat)
library(bigmemory)

setwd("~/Professional Files/Mines/SmartGeo/Queens/")
load("Data/Cleaned_Data.RData")
groundBM = attach.big.matrix("Data/groundBM.desc")
cnames = read.csv("Data/groundBM_cnames.csv", stringsAsFactors=F)[,1]
source("~/Professional Files/Mines/Research/Wind QC/Code/SNHT functions for Mandy.R")

#Basic statistics: how many obs per station?
length(unique(ground$Time))
length(unique(tunnel$Time)); dim(tunnel)
nrow(ground)/nrow(tunnel)
statNA = ddply(ground, "StationID", function(df){
  data.frame( cnt=sum(!is.na(df$Value))
             ,maxTime=max(df$Time[!is.na(df$Value)])
             ,minTime=min(df$Time[!is.na(df$Value)])
  )
} )
statNA = merge(station, statNA)
summary(statNA)
summary(station)#What?  How can we not have Northing/Easting for some stations!?
statNA$minTimePct = as.numeric(difftime(statNA$minTime,min(tunnel$Time),units="secs")) /
    as.numeric(difftime(max(tunnel$Time),min(tunnel$Time),units="secs"))
statNA$maxTimePct = as.numeric(difftime(max(tunnel$Time),statNA$maxTime,units="secs")) /
    as.numeric(difftime(max(tunnel$Time),min(tunnel$Time),units="secs"))
p = get_map( location=c(mean(statNA$Longitude, na.rm=T)+.0005, mean(statNA$Latitude, na.rm=T)), zoom=17 )
p = ggmap(p)
p1 = p + geom_point(data=statNA, aes(x=Longitude, y=Latitude, color=cnt) ) +
  labs(x="Longitude", y="Latitude", color="Observation Count")
p2 = p + geom_point(data=statNA, aes(x=Longitude, y=Latitude, color=minTimePct) ) +
  labs(x="Longitude", y="Latitude", color="First observation")
p3 = p + geom_point(data=statNA, aes(x=Longitude, y=Latitude, color=maxTimePct) ) +
  labs(x="Longitude", y="Latitude", color="Last observation")
grid.arrange(p1, p2, p3, ncol=3)
mean(statNA$cnt); max(statNA$cnt)
#Ok, so the average station only has 2000 obs. out of the total possible 11000.  The max
#is only 5759.  It also seems that observations on the west were measured more in the
#beginning, and observations on the east were measured more at the end.  The middle tends 
#to have the highest total counts.

ggsave("Results/Missing_Data.png"
  ,ggplot(ground, aes(x=Time, y=factor(StationID), fill=is.na(Value) ) ) +
      geom_tile() + labs(x="", y="Station", fill="Missing?") + scale_y_discrete(breaks=c())
  ,width=12, height=8 )

#Apply the robust SNHT to each station individually, completely ignoring spatial relationships.
temp = ground[ground$StationID==202,]
qplot( temp$Time, temp$Value )
out = snht(data=temp$Value, period=12*20, robust=T)
qplot( temp$Time, temp$Value, color=out$score )
snhtVals = dlply(ground, "StationID", function(df){
  out = snht(df$Value, period=12*20, robust=T)
  return(out$score)
})
save(snhtVals, file="Results/snhtVals.RData")
snhtMat = do.call("cbind", snhtVals)
summary( as.numeric(snhtMat) )
largeBreaks = apply(snhtMat, 2, max, na.rm=T)>1000
for(i in (1:ncol(snhtMat))[largeBreaks]){
  #Grab data from ground for this StationID only:
  g = ground[ground$StationID==as.numeric(colnames(snhtMat)[i]),]
  naFilt = !is.na(snhtMat[,i])
  print( qplot( g$Time[naFilt], g$Value[naFilt], color=snhtMat[naFilt,i]) +
    labs(title=paste("Station", colnames(snhtMat)[i])) )
  readline("Next?")
}

#Look at spatial relationships
p3 + geom_abline(slope=.2, intercept=40.748-.2*(-73.932), col="red", linetype=2 )
#Slope of region is about 20 degrees
start = Sys.time()
temp = ground[ground$Time<=ground$Time[3000],]
temp = merge(temp, station, by="StationID")
temp = temp[!is.na(temp$Longitude),]
coordinates(temp) = c("Easting", "Northing")
is(temp)
dirs = 4
v = variogram(Value ~ 1, data=temp[!is.na(temp$Value),]
    ,alpha=20-180/(2*dirs)+1:dirs*180/dirs)
ggplot(v, aes(x=dist, y=gamma)) + geom_point() +
  facet_wrap( ~ dir.hor )
Sys.time()-start



#Look at trends over time
cnames
toPlot = data.frame( groundBM[groundBM[,2] %in% as.numeric(station$StationID[1:15]),] )
colnames(toPlot) = cnames
toPlot = toPlot[,!is.na(colnames(toPlot))]
#Not sure why 1969-12-31 17:00:00 is the origin.  Maybe 1970 in UTC?
as.POSIXct(toPlot$Time[1:5], origin=as.POSIXct("1969-12-31 17:00:00"))
ground$Time[1:5]
toPlot$Time = as.POSIXct(toPlot$Time, origin=as.POSIXct("1969-12-31 17:00:00"))
pA = ggplot(toPlot, aes(x=Time, y=Value, group=StationID, color=dist.A) ) +
  geom_line(alpha=.3) + coord_cartesian(ylim=c(-2,1))
pD = ggplot(toPlot, aes(x=Time, y=Value, group=StationID, color=dist.A) ) +
  geom_line(alpha=.3) + coord_cartesian(ylim=c(-2,1))

ggplot(ground[ground$StationID==station$StationID[3],]
    ,aes(x=Time, y=Value, group=StationID) ) +
  geom_line()
