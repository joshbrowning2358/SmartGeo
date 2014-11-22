library(gstat)
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
library(GA)

fits = lapply(fits, function(x){
  x[[1]] = x[[1]][!is.na(x[[1]]$dist),]
  x[[1]] = x[[1]][x[[1]]$dist<500,]
  return(x)
})

#################################################################
# Model using fit.StVariogram function (utter failure)
#################################################################

mod = vgmST("separable", 
        space=vgm(0.2,"Exp", 20, 0.1),
        time =vgm(0.9,"Exp", 100, 0.1),
        sill=0.02)
model = fit.StVariogram(fits[[1]][[1]], mod, method="L-BFGS-B"
                ,lower=c(0,0,0,0,0) #space range, space Nug psill, time range, time Nug psill, sill
                ,upper=c(1000,1,1000,1,1), wles=TRUE)
fits[[1]][[2]] = model
plot.empVario(fits[[1]], adj=T, model="exponential", boundaries=0:65*10, rmAni=T)

#################################################################
# Optimization functions
#################################################################

#theta: c(sill, psill for spatial exp, log(range) for spatial, psill for time, log(range) for time)
mod = function(theta, h, u){
  stopifnot(length(theta)==5)
  theta[1]*
    (1-theta[2]+theta[2]*(1-exp(-h/exp(theta[3])))) *
    (1-theta[4]+theta[4]*(1-exp(-u/exp(theta[5]))))
}
WRSS = function(theta, modelFunc, emp){
  #theta
  #modelFunc
  stopifnot(is(modelFunc,"function"))
  #emp
  stopifnot(all(c("np","dist","gamma") %in% colnames(emp)))
  
  if(any(is.na(emp$dist))){
    emp = emp[!is.na(emp$dist),]
    warning("Removed rows with NA dist values in emp")
  }
  
  gam.theta=modelFunc(theta, h=emp$dist, u=emp$timelag)
  filt = gam.theta > 0
  sum(((emp$np/(gam.theta^2))*((emp$gamma-gam.theta)^2))[filt])
}
nWRSS = function(theta, modelFunc, emp) -WRSS(theta, modelFunc, emp)

#################################################################
# Optimization via optim()
#################################################################

initial = c(.02, .7, 3, .01, 3)
sol = optim( par=initial, fn=WRSS, emp=fits[[1]][[1]], modelFunc=mod
#     ,lower=c(0,.1,1,0.0001,1), upper=c(1,1,1000,.3,1000)
     ,lower=c(1e-6,.2,0,0.0001,0), upper=c(.1,.999,7,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par
fits[[1]][[2]]$space[2,2] = sol[2]
fits[[1]][[2]]$space[1,2] = 1-sol[2]
fits[[1]][[2]]$space[2,3] = exp(sol[3])
fits[[1]][[2]]$time[2,2] = sol[4]
fits[[1]][[2]]$time[1,2] = 1-sol[4]
fits[[1]][[2]]$time[2,3] = exp(sol[5])
fits[[1]][[2]]$sill = sol[1]
png("Results/Final Variograms/variogram_model_error_prd1.png")
plot.empVario(fits[[1]], model="exponential", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

#################################################################
# Model using a genetic algorithm to optimize, period 1
#################################################################

set.seed(123)
fits[[1]][[1]]$np = fits[[1]][[1]]$np/1000000
fit = ga(type="real-valued", fitness=nWRSS, emp=fits[[1]][[1]], modelFunc=mod
        ,min=c(1e-6,1e-6,0,1e-6,0), max=c(1,1,7,1,7), maxiter=1000 )
sol = slot(fit, "solution")
fits[[1]][[2]]$space[2,2] = sol[2]
fits[[1]][[2]]$space[1,2] = 1-sol[2]
fits[[1]][[2]]$space[2,3] = exp(sol[3])
fits[[1]][[2]]$time[2,2] = sol[4]
fits[[1]][[2]]$time[1,2] = 1-sol[4]
fits[[1]][[2]]$time[2,3] = exp(sol[5])
fits[[1]][[2]]$sill = sol[1]

png("Results/Final Variograms/variogram_model_error_prd1.png")
plot.empVario(fits[[1]], model="exponential", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

#################################################################
# Model using a genetic algorithm to optimize, period 2
#################################################################

set.seed(321)
fits[[2]][[1]]$np = fits[[2]][[1]]$np/1000000
fit = ga(type="real-valued", fitness=nWRSS, emp=fits[[2]][[1]], modelFunc=mod
        ,min=c(0,0,0,0,0), max=c(1,1,7,1,7), maxiter=1000 )
sol = slot(fit, "solution")
fits[[2]][[2]]$space[2,2] = sol[2]
fits[[2]][[2]]$space[1,2] = 1-sol[2]
fits[[2]][[2]]$space[2,3] = exp(sol[3])
fits[[2]][[2]]$time[2,2] = sol[4]
fits[[2]][[2]]$time[1,2] = 1-sol[4]
fits[[2]][[2]]$time[2,3] = exp(sol[5])
fits[[2]][[2]]$sill = sol[1]

png("Results/Final Variograms/variogram_model_error_prd2.png")
plot.empVario(fits[[2]], model="exponential", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

#################################################################
# Model using a genetic algorithm to optimize, period 3
#################################################################

set.seed(321)
fits[[3]][[1]]$np = fits[[3]][[1]]$np/1000000
fit = ga(type="real-valued", fitness=nWRSS, emp=fits[[3]][[1]], modelFunc=mod
        ,min=c(0,0,0,0,0), max=c(1,1,7,1,7), maxiter=1000 )
sol = slot(fit, "solution")
fits[[3]][[2]]$space[2,2] = sol[2]
fits[[3]][[2]]$space[1,2] = 1-sol[2]
fits[[3]][[2]]$space[2,3] = exp(sol[3])
fits[[3]][[2]]$time[2,2] = sol[4]
fits[[3]][[2]]$time[1,2] = 1-sol[4]
fits[[3]][[2]]$time[2,3] = exp(sol[5])
fits[[3]][[2]]$sill = sol[1]

png("Results/Final Variograms/variogram_model_error_prd3.png")
plot.empVario(fits[[3]], model="exponential", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

save(fits, file="Results/fitted_sp_variograms_error.RData")

#################################################################
# Model using a genetic algorithm and sin model, period 1
#################################################################

set.seed(123)
fits[[1]][[1]]$np = fits[[1]][[1]]$np/1000000
fit = ga(type="real-valued", fitness=nWRSS, emp=fits[[1]][[1]], modelFunc=modSin
        ,min=c(1e-6,0,0.1,0,0), max=c(1,1,10,1,7), maxiter=1000 )
sol = slot(fit, "solution")
fits[[1]][[2]] = sol

png("Results/Final Variograms/variogram_model_error_prd1.png")
plot.empVario(fits[[1]], model="sin_exp", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()


initial = c(.02, .7, 5, .01, 3)
empTemp = fits[[1]][[1]]
empTemp = empTemp[empTemp$dist<300,]
sol = optim( par=initial, fn=WRSS, emp=empTemp, modelFunc=modSin
     ,lower=c(1e-6,.1,1,0.0001,0), upper=c(.1,1-1e-4,1000,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par
fits[[1]][[2]] = sol
png("Results/Final Variograms/variogram_model_error_wave_prd1.png")
plot.empVario(fits[[1]], model="sin_exp", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

#################################################################
# Model using a genetic algorithm and sin+exponential model, period 1
#################################################################

initial = c(.02, .4, 5, 4, .01, 3)
empTemp = fits[[1]][[1]]
empTemp = empTemp[empTemp$dist<300,]
sol = optim( par=initial, fn=WRSS, emp=empTemp, modelFunc=modSinExp
     ,lower=c(1e-6,1e-6,1,0,0.0001,0), upper=c(.1,1-1e-4,1000,7,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par
fits[[1]][[2]] = sol
png("Results/Final Variograms/variogram_model_error_wave_prd1.png")
plot.empVario(fits[[1]], model="sinexp_exp", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

surf.colors <- function(x, col = terrain.colors(20)) {
  # First we drop the 'borders' and average the facet corners
  # we need (nx - 1)(ny - 1) facet colours!
  x.avg <- (x[-1, -1] + x[-1, -(ncol(x) - 1)] +
             x[-(nrow(x) -1), -1] + x[-(nrow(x) -1), -(ncol(x) - 1)]) / 4
  # Now we construct the actual colours matrix
  colors = col[cut(x.avg, breaks = length(col), include.lowest = T)]
  return(colors)
}

h = seq(0, 100, .2)
u = seq(0, 800, 10)
z <- outer(h,u, function(x,y) modSinExp(sol, x, y))
col.ind <- cut(z,100) # colour indices of each point
library(rgl)
pal <- colorRampPalette( c("blue", "green") )(100)
#persp3d(h,u,z, col=pal[col.ind])
png("Results/Final Variograms/3d_spatio_temporal_variogram_prd1.png")
persp(h, u, z, theta=-45, phi=30
    ,col=surf.colors(z, col=colorRampPalette( c("blue", "green") )(100))
    ,main="Spatio-Temporal Variogram", xlab="Distance Lag", ylab="Time Lag", zlab="Gamma")
dev.off()

#################################################################
# Model using a genetic algorithm and sin+exponential model, period 2
#################################################################

initial = c(.02, .4, 5, 0, .01, 3)
empTemp = fits[[2]][[1]]
empTemp = empTemp[empTemp$dist<300,]
sol = optim( par=initial, fn=WRSS, emp=empTemp, modelFunc=modSinExp
     ,lower=c(1e-6,0,1,0,0.0001,0), upper=c(.1,1,1000,7,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par
fits[[2]][[2]] = sol
png("Results/Final Variograms/variogram_model_error_wave_prd2.png")
plot.empVario(fits[[2]], model="sinexp_exp", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

h = seq(0, 100, .2)
u = seq(0, 800, 10)
z <- outer(h,u, function(x,y) modSinExp(sol, x, y))
col.ind <- cut(z,100) # colour indices of each point
library(rgl)
pal <- colorRampPalette( c("blue", "green") )(100)
#persp3d(h,u,z, col=pal[col.ind])
png("Results/Final Variograms/3d_spatio_temporal_variogram_prd2.png")
persp(h, u, z, theta=-45, phi=30
    ,col=surf.colors(z, col=colorRampPalette( c("blue", "green") )(100))
    ,main="Spatio-Temporal Variogram", xlab="Distance Lag", ylab="Time Lag", zlab="Gamma")
dev.off()

#################################################################
# Model using a genetic algorithm and sin+exponential model, period 3
#################################################################

initial = c(.02, .4, 5, 4, .01, 3)
empTemp = fits[[3]][[1]]
empTemp = empTemp[empTemp$dist<300,]
sol = optim( par=initial, fn=WRSS, emp=empTemp, modelFunc=modSinExp
     ,lower=c(1e-6,1e-6,1,0,0.0001,0), upper=c(.1,1-1e-4,1000,7,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par
fits[[3]][[2]] = sol
png("Results/Final Variograms/variogram_model_error_wave_prd3.png")
plot.empVario(fits[[3]], model="sinexp_exp", boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

h = seq(0, 100, .2)
u = seq(0, 800, 10)
z <- outer(h,u, function(x,y) modSinExp(sol, x, y))
col.ind <- cut(z,100) # colour indices of each point
library(rgl)
pal <- colorRampPalette( c("blue", "green") )(100)
#persp3d(h,u,z, col=pal[col.ind])
png("Results/Final Variograms/3d_spatio_temporal_variogram_prd3.png")
persp(h, u, z, theta=-45, phi=30
    ,col=surf.colors(z, col=colorRampPalette( c("blue", "green") )(100))
    ,main="Spatio-Temporal Variogram", xlab="Distance Lag", ylab="Time Lag", zlab="Gamma")
dev.off()

save(fits, file="Results/sin_sp_variograms.RData")
