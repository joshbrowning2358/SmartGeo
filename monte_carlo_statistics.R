library(mvoutlier)
if(Sys.info()[4]=="jb")
  setwd("/media/storage/Github/SmartGeo/")
if(Sys.info()[4]=="")
  
source("functions.R")
#Assume we have random normal data (probably should actually use simulate()).
#Find the 95% values of the outlying statistic:

sp = loadSp()
sp = data.frame(sp)

temp = matrix(rnorm(nrow(sp)*2),ncol=2)
temp[1,] = c(100,100)
toPlot = locoutNeighbor(temp, X=sp$Easting, Y=sp$Northing
        ,propneighb=beta, usemax=1, npoints=50, variant="dist" )
plot(c(0,50), c(0,1), type="n")
for(i in 1:100)
  lines(1:50, toPlot[i,-1], col=gray(0.7))
lines(1:50,rep(beta,50),col=2)
mvoutlier::locoutNeighbor(temp, X=sp$Easting, Y=sp$Northing
        ,propneighb=beta, usemax=1, npoints=50, variant="dist", indices=1 )

#Get the highest and lowest statistics for each time point for 1000 simulations
hilo = NULL
for(i in 1:1000){
  for(beta in 1:3/10){
    toPlot = locoutNeighbor(matrix(rnorm(nrow(sp)*2),ncol=2), X=sp$Easting, Y=sp$Northing
            ,propneighb=beta, usemax=1, npoints=50, variant="dist" )
    hi = apply(toPlot[,-1], 2, max)
    lo = apply(toPlot[,-1], 2, min)
    hilo = rbind(hilo, data.frame(beta=beta, points=1:50, hi=hi, lo=lo) )
  }
}
save(hilo, file="Results/monte_carlo_isolation.RData")
ggplot( hilo, aes(x=points, fill=factor(beta), group=factor(beta)) ) +
  stat_summary(aes(y=hi, color=factor(beta)), fun.y=max, geom="line") +
  stat_summary(aes(y=lo, color=factor(beta)), fun.y=min, geom="line")
