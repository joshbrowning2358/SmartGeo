library(RandomFields)
library(ggplot2)
library(plyr)
library(dplyr)
library(spdep)
library(gstat)
library(spacetime)
library(sqldf)
library(CompRandFld)
library(fields)
library(mgcv)

if(Sys.info()[4]=="jb")
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
if(Sys.info()[4]=="JOSH_LAPTOP")
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
load("Data/Cleaned_Data.RData")
rm(tunnel); gc()


#######################################################################
# Basic statistics: how many observations per station?
#######################################################################

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


#######################################################################
#Apply the robust SNHT to each station individually, completely ignoring spatial relationships.
#######################################################################

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

#######################################################################
# Investigate times when tunnels are close to data
#######################################################################

load("Data/Cleaned_Data.RData")
YL = grepl("\\.YL", colnames(tunnel))
A = grepl("\\.A", colnames(tunnel))
D = grepl("\\.D", colnames(tunnel))
BC = grepl("\\.BC", colnames(tunnel))
minYLdist = apply( tunnel[,YL], 1, min, na.rm=T)
minAdist = apply( tunnel[,A], 1, min, na.rm=T)
minDdist = apply( tunnel[,D], 1, min, na.rm=T)
minBCdist = apply( tunnel[,BC], 1, min, na.rm=T)
bnd1 = geom_vline(xintercept=800)
bnd2 = geom_vline(xintercept=4500)
qplot( 1:length( minYLdist ), minYLdist, geom="line" ) + bnd1 + bnd2
qplot( 1:length( minAdist ), minAdist, geom="line" ) + bnd1 + bnd2
qplot( 1:length( minDdist ), minDdist, geom="line" ) + bnd1 + bnd2
qplot( 1:length( minBCdist ), minBCdist, geom="line" ) + bnd1 + bnd2
prd1 = 1:800
prd2 = 801:4500
prd3 = 4501:nrow(tunnel)

#######################################################################
# Investigate slope of region, determine transformation
#######################################################################

angle = function( pt1, pt2 ){
  deltaX = pt1[1]-pt2[1]
  deltaY = pt1[2]-pt2[2]
  atan(deltaY/deltaX)
}
rows = (1:nrow(station))[!is.na(station$Easting)]
angleVec = c()
for(i in 1:length(rows)){
#  sapply( (i+1):length(rows), function(j){
  for( j in (i+1):length(rows)){
    angleVec = c(angleVec, angle(station[rows[i],c("Easting","Northing")]
                                ,station[rows[j],c("Easting","Northing")] ) )
  }
  cat("Finished",i,"\n")
}

slope=.249
ggplot( station, aes(x=Easting, y=Northing) ) + geom_point() +
  geom_abline(intercept=30+211714.3-slope*1003357, slope=slope)
theta = atan(slope)
theta*180/pi #About 14 degrees off of East
#Rotate coordinate system by theta degrees counterclockwise:
station$East2 = station$Easting*cos(theta) + station$Northing*sin(theta)
station$North2 = -station$Easting*sin(theta) + station$Northing*cos(theta)
ggplot( station, aes(x=East2, y=North2) ) + geom_point()
save(station, file="Data/station_new_coords.RData")

#######################################################################
# Play with RFsim from CompRandFld package
#######################################################################

n = 20
nT = 20
x = 1:n
y = 1:n
t = 1:nT

d = RFsim(x, y, t, corrmodel="exp_exp", grid=TRUE
  ,param=list(nugget=0, mean=0,scale_s=2, scale_t=2,sill=1)
)$data
for(i in 1:dim(d)[3]){
  image.plot(1:n, 1:n, d[,,i], zlim=c(min(d),max(d)))
  Sys.sleep(1)
}
mod = FitComposite(data=d[,,1], coordx=1:n, coordy=1:n, corrmodel="exponential", grid=TRUE)
mod = FitComposite(data=d[,,2], coordx=1:n, coordy=1:n, corrmodel="exponential", grid=TRUE)
mod = FitComposite(data=d[,,3], coordx=1:n, coordy=1:n, corrmodel="exponential", grid=TRUE)
mod = FitComposite(data=d[,,4], coordx=1:n, coordy=1:n, corrmodel="exponential", grid=TRUE)
mod = FitComposite(data=d, coordx=1:n, coordy=1:n, coordt=1:nT
  ,corrmodel="exp_exp", grid=TRUE)

#######################################################################
# Model deformation as a function of distance
#######################################################################

source_github("https://raw.githubusercontent.com/rockclimber112358/Random_Code/master/plot_gam.R")
load("Data/ground_with_distance.RData")
load("Data/station_new_coords.RData")
station = station[!is.na(station$East2),]
ground = merge(ground, station[,c("StationID", "East2", "North2")]
        ,by="StationID")
ground = ground[!is.na(ground$Value),] #remove cases without deformation data
fit = gam( Value ~ s(A) + s(BC) + s(D) + s(YL), data=ground )
summary(fit)

#######################################################################
# Fit Spatio-temporal variograms with real data!
#######################################################################

#Slope of region is about 20 degrees
start = Sys.time()
temp = ground[ground$Time<=ground$Time[1000],]
temp = merge(temp, station, by="StationID")
temp = temp[!is.na(temp$Longitude),]
coordinates(temp) = c("Easting", "Northing")
dirs = 4
v = variogram(Value ~ 1, data=temp[!is.na(temp$Value),]
    ,alpha=20-180/(2*dirs)+1:dirs*180/dirs)
ggplot(v, aes(x=dist, y=gamma)) + geom_point() +
  facet_wrap( ~ dir.hor )
Sys.time()-start
save(v, file="Results/Variogram_t_1:1000.RData")

temp = ground[ground$Time<=ground$Time[10],]
temp = merge(temp, station, by="StationID")
temp = temp[!is.na(temp$Longitude),]
coordinates(temp) = c("Easting", "Northing")
v2 = variogram(Value ~ 1, data=temp[!is.na(temp$Value),]
    ,alpha=20-180/(2*dirs)+1:dirs*180/dirs, cloud=TRUE)
ggplot( data.frame(v2), aes(x=dist, y=gamma) ) +
  geom_point(alpha=.002) + facet_wrap( ~ dir.hor) + geom_smooth() +
  coord_cartesian(ylim=c(0,.05))

v2 = variogram(Value ~ 1, data=temp[!is.na(temp$Value),], cloud=TRUE, cressie=TRUE)
ggplot( data.frame(v2), aes(x=dist, y=gamma) ) +
  geom_point(alpha=.002) + facet_wrap( ~ dir.hor) + geom_smooth()

temp = ground[ground$Time<=ground$Time[1000],]
temp = merge(temp, station, by="StationID")
temp = temp[!is.na(temp$Longitude),]
coordinates(temp) = c("Easting", "Northing")
v3 = variogram(Value ~ 1, data=temp[!is.na(temp$Value),], cressie=TRUE)
ggplot( data.frame(v3), aes(x=dist, y=gamma) ) +
  geom_point(alpha=1)

timeVg = dlply( ground[ground$Time<=ground$Time[100],], "Time", function(df){
    df = merge(df, station, by="StationID")
    df = df[!apply(df,1,function(x) any(is.na(x))),]
    coordinates(df) = c("Easting", "Northing")
    out = data.frame(variogram(Value ~ 1, data=df, cloud=T))
    out$Time = df$Time[1]
    return(out)
  } )
#timeVg2 = dlply( ground[ground$Time<=ground$Time[1],], "Time", function(df){
#    df = merge(df, station, by="StationID")
#    df = df[!apply(df,1,function(x) any(is.na(x))),]
#    out = data.frame( dist=as.numeric( dist( df[,c("Easting", "Northing")] ) )
#                     ,gamma=as.numeric( dist(df$Value)^2 ) )/2
#    return( out )
#  } )
timeVg = do.call("rbind", timeVg)

#Variogram appears to be flat, and so it may be ok to assume spatial independence
ggplot(timeVg[sample(1:nrow(timeVg),size=5000),], aes(x=dist, y=gamma) ) +
  geom_point() + geom_smooth() + coord_cartesian( ylim=c(0,.1) )



#################################################################
# Isotropy check
#################################################################

times = unique(ground$Time)
timeVg = NULL
for(i in 6136:length(times)){
  df = ground[ground$Time==times[i],]
  df = merge(df, station, by="StationID")
  df = df[!apply(df,1,function(x) any(is.na(x))),]
  if(nrow(df)==0)
    next;
  coordinates(df) = c("Easting", "Northing")
  out = try(data.frame(variogram(Value ~ 1, data=df, cloud=T, alpha=0:35*10)))
  if(nrow(out)==0)
    next;
  if(is(out,"try-error"))
    next;
  out$Time = df$Time[1]
  timeVg = rbind(timeVg, out)
  cat("Iteration",i,"completed out of",length(times),"\n")
  if(i %% 100 == 0){
    save(timeVg, file=paste0("Results/time_variogram_times_",i-99,"_",i) )
    timeVg = NULL
  }
}
save(timeVg, file=paste0("Results/time_variogram_times_",i) )

files = list.files("Results/")
files = files[grepl("time_variogram",files)]
meanVg = NULL
for(file in files){
  load(paste0("Results/",file))
  timeVg$dist2 = floor(timeVg$dist/10)*10
  meanVgTemp = ddply( timeVg, c("dir.hor", "dist2"), function(df){
    data.frame( gamma=mean(df$gamma), np=nrow(df) )
  } )
  meanVgTemp$file = file
  meanVg = rbind(meanVg, meanVgTemp)
}
meanVg$i = as.numeric( gsub("(time_variogram_times_|_[0-9]*)", "", meanVg$file) )
ggplot( meanVg, aes(x=dist2, y=gamma, color=i, group=i) ) + geom_line() +
  facet_wrap( ~ dir.hor, scale="free" )

#Aggregate over all time points, estimate the overall variogram and sd.
vgTotal = ddply(meanVg, c("dir.hor", "dist2"), function(df){
  out = data.frame(gamma=sum(df$gamma*df$np)/sum(df$np)
            ,np=sum(df$np) )
  out$sd = sqrt(2*out$gamma^2/out$np)
  return(out)
} )
ggsave("Results/anisotropic_investigation.png",
  ggplot( vgTotal[vgTotal$dir.hor<180,], aes(x=dist2) ) + geom_line(aes(y=gamma)) +
    geom_ribbon( aes(ymin=gamma-sd, ymax=gamma+sd), alpha=.2 ) +
    facet_wrap( ~ dir.hor ) + coord_cartesian(ylim=c(0,.2))
)

#Aggregate over all time points, estimate the overall variogram and sd.
vgTotal = ddply(meanVg[meanVg$i<2000,], c("dir.hor", "dist2"), function(df){
  out = data.frame(gamma=sum(df$gamma*df$np)/sum(df$np)
            ,np=sum(df$np) )
  out$sd = sqrt(2*out$gamma^2/out$np)
  return(out)
} )
ggsave("Results/anisotropic_investigation_pretunnel.png",
  ggplot( vgTotal[vgTotal$dir.hor<180,], aes(x=dist2) ) + geom_line(aes(y=gamma)) +
    geom_ribbon( aes(ymin=gamma-sd, ymax=gamma+sd), alpha=.2 ) +
    facet_wrap( ~ dir.hor )# + coord_cartesian(ylim=c(0,.2))
)

meanVg$dir.hor2 = floor(meanVg$dir.hor/45)*45
meanVg$dir.hor2[meanVg$dir.hor2>=180] = meanVg$dir.hor2[meanVg$dir.hor2>=180]-180
meanVgSub = ddply(meanVg, c("dir.hor2", "dist2", "i"), function(df){
  data.frame( gamma=mean(df$gamma, na.rm=T) )
} )
ggplot( meanVgSub, aes(x=dist2, y=gamma, color=i, group=i) ) + geom_line(alpha=.2) +
  facet_wrap( ~ dir.hor2, scale="free" ) + geom_smooth(aes(group=NA)) +
  labs(x="Distance", y="Variogram", color="Time Step")

#################################################################
# Space-time modeling
#################################################################

#Create spatial object

sp = ground[ground$Time==ground$Time[1],]
sp = merge(sp, station, by="StationID")
sp = sp[!is.na(sp$Longitude),]
validStations = as.character(sp$StationID)
sp = sp[,c("Easting", "Northing")]
coordinates(sp) = c("Easting", "Northing")

time = ground$Time[1:10]

data = filter(ground, Time %in% time, StationID %in% validStations)
data = arrange(data, Time, StationID)
stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
plot(stdf)
vst = variogramST( Value ~ 1, data=stdf )
png("Results/spatio-temporal_empirical_variogram_first_10.png")
plot(vst)
dev.off()

time = ground$Time[1:100]

data = filter(ground, Time %in% time, StationID %in% validStations)
data = arrange(data, Time, StationID)
stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
vst = variogramST( Value ~ 1, data=stdf, tlags=1:50 )
png("Results/spatio-temporal_empirical_variogram_first_100.png")
plot(vst)
library(rgl)
x = unique( round( vst$dist ) )
y = unique( vst$timelag )
z = expand.grid(x, y)
z$z = NA
for(i in 1:nrow(z))
  z$z[i] = vst[round(vst$dist)==z$Var1[i] & vst$timelag==z$Var2[i]]
rgl.surface(x=x, y=y, z=vst$gamma)
dev.off()

time = ground$Time[1:1000]

data = filter(ground, Time %in% time, StationID %in% validStations)
data = arrange(data, Time, StationID)
stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
vst = variogramST( Value ~ 1, data=stdf )
png("Results/spatio-temporal_empirical_variogram_first_1000.png")
plot(vst)
dev.off()


for( prd in list(prd1, prd2, prd3) ){

  time = unique(ground$Time[prd])

  data = filter(ground, Time %in% time, StationID %in% validStations)
  data = arrange(data, Time, StationID)
  stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
  vst = variogramST( Value ~ 1, data=stdf, tlags=0:28*9 )
  png(paste0("Results/spatio-temporal_empirical_variogram_prd_",min(prd),".png"))
  plot(vst)
  dev.off()
  save(vst, file=paste0("Results/spatio-temporal_empirical_variogram_prd_",min(prd),".RData"))
}

data = filter(ground, Time %in% time, StationID %in% validStations)
data = arrange(data, Time, StationID)
stdf = STFDF(sp=sp, time=time, data=data[,"Value",drop=F] )
vst = variogramST( Value ~ 1, data=stdf, tlags=0:28*9 )
png("Results/spatio-temporal_empirical_variogram_all.png")
plot(vst)
dev.off()
save(vst, file="Results/spatio-temporal_empirical_variogram_all.RData")
load("/media/storage/Github/SmartGeo/spatio-temporal_empirical_variogram_all.RData")
is(vst)
vst = data.frame(vst)

ggplot( data.frame(vst), aes(x=dist, y=gamma, color=timelag, group=timelag) ) +
  geom_line() + labs(x="Spatial Distance", y="Varigoram", color="Time Lag")
ggplot( data.frame(vst), aes(x=timelag, y=gamma, color=spacelag, group=spacelag) ) +
  geom_line() + labs(x="Spatial Distance", y="Varigoram", color="Time Lag")

mod = vgmST("separable"
         ,space=vgm(psill=.90,"Exp", range=600, nugget=0.4),
         ,time =vgm(psill=.90,"Exp", range=40, nugget=0.4),
         ,sill=0.9)

vstModel = fit.StVariogram(vst, model=mod)
mod = vstModel$space
dist = seq(0,1700,10)
qplot( dist, mod$psill[mod$model=="Nug"] + mod$psill[mod$model=="Exp"]*
         (1-exp(-dist/(mod$range[mod$model=="Exp"]))) ) +
  labs(x="Distance", y="Modeled variogram")

mod = vstModel$time
time = 0:100
qplot( time, mod$psill[mod$model=="Nug"] + mod$psill[mod$model=="Exp"]*
         (1-exp(-time/(mod$range[mod$model=="Exp"]))) ) +
  labs(x="Time Lag", y="Modeled variogram")

time = unique(ground$Time)
grid = as.matrix(station[,c("East2", "North2")])[!is.na(station$East2),]
data = RFsim(coordx=grid, corrmodel="exponential"
    ,grid=FALSE, param=list(nugget=.693, mean=0, sill=1, scale=1e10) )
data$data
qplot( grid[,1], grid[,2], color=data$data )
