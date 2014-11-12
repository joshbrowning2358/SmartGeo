library(dplyr)

getDist = function(statCoord){
  stopifnot(is(statCoord,"data.frame"))
  stopifnot(colnames(statCoord)==c("x","y","s"))

  statDist = merge(statCoord, statCoord, by=NULL, suffixes=c("1","2"))
  statDist = mutate( statDist
              ,dist=sqrt((x1-x2)^2+(y1-y2)^2)
              ,angle=atan2(y2-y1,x2-x1))
  #Only add pi because we want all angles in (0,180)
  statDist$angle = ifelse(statDist$angle<0,statDist$angle+pi,statDist$angle)*180/pi
  statDist$x1 = NULL
  statDist$x2 = NULL
  statDist$y1 = NULL
  statDist$y2 = NULL
  return(statDist)
}

#data: data.frame containing the variables of interest
#statCol: column name of the station id variable in data
#timeCol: column name of the time variable in data
#dataCol: column name of the observed variable in data
#statDist: data.frame with four columns:
#  "s1": id of the first location, must have the same ids as data
#  "s2": id of the first location, must have the same ids as data
#  "dist": the distance between the two stations
#  "angle": the angle between the two stations
#distInt: vector of distances for grouping the distance
#timeInt: vector of times for grouping the time
#angleInt: vector of angles for grouping the angle
variogramPts <- function(data, statCol, timeCol, dataCol, statDist
                ,distInt=quantile(statDist$dist,0:10/20)
                ,timeInt=(max(data[,timeCol])-min(data[,timeCol]))*0:10/20
                ,angleInt=c(0,180)){
  #data
  stopifnot(is(data,"data.frame"))
  stopifnot( all(c(statCol,timeCol,dataCol) %in% colnames(data)) )
  #*Col
  stopifnot(is(statCol,"character"))
  stopifnot(is(timeCol,"character"))
  stopifnot(is(dataCol,"character"))
  #statDist
  stopifnot(all(colnames(statDist)==c("s1", "s2", "dist", "angle")))
  stopifnot(all(statDist[,"s1"] %in% data[,statCol]))
  stopifnot(all(statDist[,"s2"] %in% data[,statCol]))
  stopifnot(all(data[,statCol] %in% statDist[,"s1"]))
  stopifnot(all(data[,statCol] %in% statDist[,"s2"]))
  #*Int
  stopifnot(is(distInt, "numeric"))
  stopifnot(all(diff(distInt)>=0))
  stopifnot(is(timeInt, "numeric"))
  stopifnot(all(diff(timeInt)>=0))
  stopifnot(is(angleInt, "numeric"))
  stopifnot(all(diff(angleInt)>=0))
  stopifnot(max(angleInt)<=180 & min(angleInt)>=0)
  
  maxTimeDiff = max(timeInt)
  colnames(data)[colnames(data)==statCol] = "s"
  colnames(data)[colnames(data)==timeCol] = "t"
  colnames(data)[colnames(data)==dataCol] = "z"
  
  out = NULL
  for(time in unique(data$t)){
    data1 = filter(data, t==time)
    data2 = filter(data, t<=time & t>time-maxTimeDiff)
    pairs = merge(data1, data2, by=NULL, suffixes=1:2)
    pairs = merge(pairs, statDist, by=c("s1", "s2") )
    pairs$deltaT = abs(pairs$t1-pairs$t2)
    pairs$tGrp = findInterval(pairs$deltaT, timeInt)
    pairs$dGrp = findInterval(pairs$dist, distInt)
    pairs$aGrp = findInterval(pairs$angle, angleInt)
    
    outTemp = group_by(pairs, tGrp, dGrp, aGrp)
    outTemp = data.frame( summarize(outTemp
        ,vEst = mean( (z1-z2)^2 )
        ,cnt = n()
        ,meanD = mean(dist)
        ,meanT = mean(deltaT)
        ,meanA = mean(angle)
    ) )
    
    #Add current results into total variogram estimate:
    out = rbind(out, outTemp)
  }
  
  #Collapse out data.frame by grouping the three buckets
  out = group_by(out, tGrp, dGrp, aGrp)
  out = data.frame( summarize(out
      ,vEst = sum(vEst*cnt)/sum(cnt)
      ,cnt = sum(cnt)
      ,meanD = sum(meanD*cnt)/sum(cnt)
      ,meanT = sum(meanT*cnt)/sum(cnt)
      ,meanA = mean(meanA*cnt)/sum(cnt)
  ) )
  
  return(out)
}

#Examples of running the analysis:
# data = data.frame(stat=rep(1:10,each=10), time=rep(1:10,times=10), data=rnorm(100))
# statCoord = data.frame(x=runif(10), y=runif(10), s=1:10)
# statDist = getDist(statCoord)
# v = variogramPts(data, statCol="stat", timeCol="time", dataCol="data", statDist)
# ggplot(v, aes(x=dGrp, y=vEst, color=tGrp, group=tGrp) ) + geom_line()
# 
# load("~/Professional Files/Mines/SmartGeo/Queens/Data/ground_with_nnet.RData")
# load("~/Professional Files/Mines/SmartGeo/Queens/Data/station_pairs.RData")
# colnames(station_pairs) = c("s1", "s2", "dist", "angle")
# ground$Time = as.numeric(ground$Time)
# 
# groundSmall = ground[ground$Time<sort(unique(ground$Time))[800],]
# statDistSmall = station_pairs[station_pairs$s1 %in% groundSmall$StationID,]
# statDistSmall = statDistSmall[statDistSmall$s2 %in% groundSmall$StationID,]
# start = Sys.time()
# v2 = variogramPts(groundSmall
#           ,statDist=statDistSmall
#           ,statCol="StationID", timeCol="Time", dataCol="Value"
#           ,distInt=0:65*10
#           ,timeInt=0:10*(60*60*2.4) #One time unit = 2.4 hours
#           ,angleInt=c(0,180) )
# Sys.time() - start
# ggplot(v2, aes(x=dGrp, y=vEst, color=tGrp, group=tGrp) ) + geom_line()
# ggplot(v2, aes(x=tGrp, y=vEst, color=dGrp, group=dGrp) ) + geom_line()
