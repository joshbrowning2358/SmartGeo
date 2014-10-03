library(RandomFields)
library(ggplot2)
library(plyr)
library(dplyr)
library(spdep)
library(gstat)
library(spacetime)
library(sqldf)
library(CompRandFld)

setwd("~/Professional Files/Mines/SmartGeo/Queens/")
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
rm(tunnel); gc() #Won't be needing it for this analysis

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

#######################################################################
# On to real data!
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
load("Results/spatio-temporal_empirical_variogram_all.RData")
is(vst)


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
RFsimulate( vstModel, x=coordinates(sp)
  ,T=c(min(ground$Time), 2.4, max(ground$Time)) )
