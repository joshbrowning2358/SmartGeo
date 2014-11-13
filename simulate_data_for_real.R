library(CompRandFld)
if(Sys.info()[4]=="jb"){
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
  source("/media/storage/Github/SmartGeo/functions.R")
}
if(Sys.info()[4]=="JOSH_LAPTOP"){
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
  source("~/Github/SmartGeo/functions.R")
}
if(grepl("ch120", Sys.info()[4])){
  setwd("~/SmartGeo/Queens")
  source("~/Github/SmartGeo/functions.R")
}
load("Results/new_sp_variograms_error.RData")

fits = lapply(fits, function(x){
  x[[1]] = x[[1]][!is.na(x[[1]]$dist),]
  x[[1]] = x[[1]][x[[1]]$dist<500,]
  return(x)
})
load("Results/fitted_sp_variograms_error.RData")
load("Data/station_new_coords.RData")
station = station[!is.na(station$East2),]

plot.empVario(fits[[1]], adj=T, rmAni=T, boundaries=0:65*10, model="exponential")
plot.empVario(fits[[2]], adj=T, rmAni=T, boundaries=0:65*10, model="exponential")
plot.empVario(fits[[3]], adj=T, rmAni=T, boundaries=0:65*10, model="exponential")

fits[[1]][[2]]
sim = RFsim(coordx=as.matrix(station[1:371,c("East2","North2")]), coordt=0:80*24
      ,corrmodel="exp_exp", grid=FALSE
      ,param=list(nugget=0.3, mean=0, scale_s=988.5665, scale_t=4.721552, sill=0.03278) )

start = Sys.time()
sim = RFsim(coordx=as.matrix(station[1:371,c("East2","North2")]), coordt=0:100*24
      ,corrmodel="exp_exp", grid=FALSE
      ,param=list(nugget=0.3, mean=0, scale_s=988.5665, scale_t=4.721552, sill=0.03278) )
Sys.time()-start