library(ggmap)
library(rgdal)

if(Sys.info()[4]=="jb")
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
if(Sys.info()[4]=="JOSH_LAPTOP")
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
load("Data/Cleaned_Data.RData")
source_github("https://raw.githubusercontent.com/rockclimber112358/Random_Code/master/gg_color_hue.R")

#Basic plot of station locations and tunnel movement
toPlot = ground[ground$Time==ground$Time[1],]
toPlot = merge(toPlot, station, by="StationID")
p = get_map( location=c(mean(toPlot$Longitude, na.rm=T), mean(toPlot$Latitude, na.rm=T)), zoom=17 )
p = ggmap(p)
p + geom_point(data=toPlot[is.na(toPlot$Value),], aes(x=Longitude, y=Latitude), color="red", alpha=.2) +
  geom_point(data=toPlot[!is.na(toPlot$Value),], aes(x=Longitude, y=Latitude, color=Value) )
p + geom_point(data=toPlot[!is.na(toPlot$Value),], aes(x=Longitude, y=Latitude), size=.1 )

tunnelYL = unique(tunnel[,c("EastingYL", "NorthingYL")])
tunnelYL = cbind(tunnelYL
  ,project(as.matrix(tunnelYL), proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
colnames(tunnelYL) = c("Easting", "Northing", "Latitude", "Longitude")

tunnelBC = unique(tunnel[,c("EastingBC", "NorthingBC")])
tunnelBC = cbind(tunnelBC
  ,project(as.matrix(tunnelBC), proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
colnames(tunnelBC) = c("Easting", "Northing", "Latitude", "Longitude")

tunnelA = unique(tunnel[,c("EastingA", "NorthingA")])
tunnelA = cbind(tunnelA
  ,project(as.matrix(tunnelA), proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
colnames(tunnelA) = c("Easting", "Northing", "Latitude", "Longitude")

tunnelD = unique(tunnel[,c("EastingD", "NorthingD")])
tunnelD = cbind(tunnelD
  ,project(as.matrix(tunnelD), proj=slot(CRS("+init=ESRI:102718"), "projargs"), inv=T ) )
colnames(tunnelD) = c("Easting", "Northing", "Latitude", "Longitude")

col = gg_color_hue(4)
ggsave("Results/Spatial_Plot.png",
  p + geom_point(data=toPlot[!is.na(toPlot$Value),], aes(x=Longitude, y=Latitude, color="Stations") ) +
    geom_point(data=tunnelYL, aes(x=Latitude, y=Longitude, color="YL"), size=1 ) +
    geom_point(data=tunnelA, aes(x=Latitude, y=Longitude, color="A"), size=1 ) +
    geom_point(data=tunnelD, aes(x=Latitude, y=Longitude, color="D"), size=1 ) +
    geom_point(data=tunnelBC, aes(x=Latitude, y=Longitude, color="BC"), size=1 ) +
    scale_color_manual(breaks=c("Stations", "YL", "A", "BC", "D")
      ,values=c(Stations="#000000", YL=col[1], A=col[2], BC=col[3], D=col[4])) +
    labs(x="Longitude", y="Latitude", color="")
)