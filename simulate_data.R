library(RandomFields)
library(ggplot2)
library(scales)
library(plyr)
library(dplyr)
library(spdep)
library(gstat)
library(spacetime)
library(sqldf)
library(CompRandFld)
library(fields)
library(mgcv)
library(MASS)
library(snht)

if(Sys.info()[4]=="jb")
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
if(Sys.info()[4]=="JOSH_LAPTOP")
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
source("~/GitHub/SmartGeo/functions.R")
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
#Plot each station individually, look for errors of both types
#######################################################################

load("Data/ground_with_distance.RData")
for(stat in unique(ground$StationID) ){
  temp = ground[ground$StationID==stat & !ground$outlier,]
  temp$closeTunnel = ifelse( temp$A<300, "A"
                          ,ifelse(temp$BC<300, "BC"
                          ,ifelse(temp$D<300, "D"
                          ,ifelse(temp$YL<300, "YL", "None" ) ) ) )
  temp$closeTunnel = factor(temp$closeTunnel, levels=c("YL", "A", "BC", "D", "None") )
  temp$closeTunnel[temp$closeTunnel=="None"] = NA
  temp$deltaValue = c(diff(temp$Value), NA)
  p = ggplot(temp, aes(x=Time) ) + 
    geom_line(aes(y=Value), alpha=.8) +
    coord_cartesian(y=c(min(temp$Value, na.rm=T)-.1, max(temp$Value, na.rm=T)+.1) ) +
    labs(x="", y="Deformation", fill="Tunnel", title=paste("Station",stat)) +
    guides(alpha=F)
  if(sum(!is.na(temp$closeTunnel))>0)
    p = p + geom_tile(data=temp[!is.na(temp$closeTunnel),],
      aes(y=0,fill=closeTunnel, alpha=.1, height=1000))

  xRng = as.POSIXct( range(temp$Time[!is.na(temp$Value)]) )
  ggsave(paste0("Results/Univariate_Time_Series/Station_",stat,".png"),
    p + coord_cartesian(y=c(-0.5,0.5), x=xRng)
  )
  p = ggplot(temp, aes(x=Time) ) + 
    geom_line(aes(y=deltaValue), alpha=.8) +
    coord_cartesian(y=c(min(temp$deltaValue, na.rm=T)-.1, max(temp$deltaValue, na.rm=T)+.1) ) +
    labs(x="", y="Deformation", fill="Tunnel", title=paste("Station",stat,"Deltas")) +
    guides(alpha=F)
  if(sum(!is.na(temp$closeTunnel))>0)
    p = p + geom_tile(data=temp[!is.na(temp$closeTunnel),],
      aes(y=0,fill=closeTunnel, alpha=.1, height=1000))
  ggsave(paste0("Results/Univariate_Time_Series/Station_",stat,"_delta.png"),
    p + coord_cartesian(y=c(-0.5,0.5))
  )
}

#######################################################################
#Apply the robust SNHT to each station individually, completely ignoring spatial relationships.
#######################################################################

temp = ground[ground$StationID==202,]
snhtVals = dlply(ground[!ground$outlier,], "StationID", function(df){
  out = snht(df$Value, period=10*60, robust=T)
  return(out$score)
})
save(snhtVals, file="Results/snhtVals.RData")
#load("Results/snhtVals.RData")
snhtMat = do.call("cbind", snhtVals)
summary( as.numeric(snhtMat) )
largeBreaks = apply(snhtMat, 2, max, na.rm=T)>1000
for(i in (1:ncol(snhtMat))[largeBreaks]){
  #Grab data from ground for this StationID only:
  g = filter( ground, StationID==as.numeric(colnames(snhtMat)[i]) )
  naFilt = !is.na(snhtMat[,i])
  print( qplot( g$Time[naFilt], g$Value[naFilt], color=snhtMat[naFilt,i]) +
    labs(title=paste("Station", colnames(snhtMat)[i])) )
  readline("Next?")
}

thresh = 100
tunnel = group_by(ground, Time)
tunnelDist = summarize(tunnel
  ,A=mean(A, na.rm=T)
  ,BC=mean(BC, na.rm=T)
  ,D=mean(D, na.rm=T)
  ,YL=mean(YL, na.rm=T) )
tunnelDist$closeTunnel = ifelse(tunnelDist$A<500, "A"
                        ,ifelse(tunnelDist$BC<500, "BC"
                        ,ifelse(tunnelDist$D<500, "D"
                        ,ifelse(tunnelDist$YL<500, "YL", NA ) ) ) )
#What % exceed the threshold as a function of time?
qplot(tunnelDist$Time
     ,apply(snhtMat, 1, function(x){ sum(x>thresh, na.rm=T)/sum(!is.na(x)) })
     ,geom="line") +
  geom_tile(data=tunnelDist
      ,aes(x=Time, y=0.5, fill=closeTunnel, alpha=is.na(closeTunnel)), height=1) +
  scale_alpha_manual(breaks=c(T,F), values=c(.5,0)) + guides(alpha=F) +
  coord_cartesian(ylim=c(0,.1)) +
  scale_y_continuous("Stations exceeding threshold", label=percent) 

#TO-DO!  Look for spatio-temporal relationships in snht values
#TO-DO!  Examine observed change points, create function to simulate similar errors.

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

ground = loadGround(timeCnt=800)
colnames(ground) = c("s", "t", "value")
ground$t = as.numeric(ground$t)/(24*60*60)
s = data.frame(loadSp())
colnames(s)[1] = "s"
vst = variogramPts( d=ground, s=s, maxs=0.5, maxt=1, s_brks=0:20*100, t_brks=0:30, angle_brks=0:4*45 )
write.csv(vst, file="variogram.csv")

rm(ground)
ground = loadGround(timeCnt=12000)
fits = list()
for( prd in list(prd1, prd2, prd3) ){
  fit = empVario(data=ground[ground$Time %in% unique(ground$Time)[prd],],tlags=0:30*9)
  fits[[length(fits)+1]] = fit
  print(plot.empVario(fit))
}
save(fits, prd1, prd2, prd3, file="Results/sp_variograms.RData")
#load("Results/sp_variograms.RData")

plot.empVario(fits[[1]])
plot.empVario(fits[[2]])

#Refit model using vgm and fit.StVariogram
vst = fits[[2]][[1]]
mod = vgmST("separable"
      ,space=vgm(psill=0, "Sph", range=100, nugget=1),
      ,time =vgm(psill=1, "Sph", range=70, nugget=0),
      ,sill=max(vst$gamma, na.rm=T))
#vstModel = fit.StVariogram(vst, model=mod
#  ,lower=rep(0,5), upper=c(100, .004, 10, .005, .01), method="L-BFGS-B")
vstModel = fit.StVariogram(vst, model=mod
  ,lower=rep(0,5), upper=c(1000, 0.6, 200, .3, .05), method="L-BFGS-B")
fits[[2]][[2]] = vstModel
plot.empVario(fits[[2]])

#Refit model using fitModelST (my optim function)
mod = fitModelST(vst, initial_t=c(0.1, 50))
mod = list(space=data.frame(model=c("Nug", "Sph")
                           ,psill=c(mod[[1]][1], 1-mod[[1]][1])
                           ,range=c(0,mod[[1]][2]) )
          ,time=data.frame( model=c("Nug", "Sph")
                           ,psill=c(mod[[2]][1], 1-mod[[2]][1])
                           ,range=c(0,mod[[2]][2]) )
          ,sill=mod[[3]])
fits[[2]][[2]] = mod
plot.empVario(fits[[2]])
