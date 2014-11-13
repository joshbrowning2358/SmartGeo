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
# Genetic Algorithm optimization functions
#################################################################

#theta: c(sill, psill for spatial exp, log(range) for spatial, psill for time, log(range) for time)
mod = function(theta, h, u){
  stopifnot(length(theta)==5)
  theta[1]*
    (1-theta[2]+theta[2]*(1-exp(-h/exp(theta[3])))) *
    (1-theta[4]+theta[4]*(1-exp(-h/exp(theta[5]))))
}
WRSS = function(theta, modelFunc, emp){
  #theta
  #modelFunc
  stopifnot(is(modelFunc,"function"))
  #emp
  stopifnot(all(c("np","dist","gamma") %in% colnames(emp)))
  
  gam.theta=modelFunc(theta, emp$dist, emp$timelag)
	sum((emp$np/(gam.theta^2))*((emp$gamma-gam.theta)^2))
}
nWRSS = function(theta, modelFunc, emp) -WRSS(theta, modelFunc, emp)

optim( par=initial, fn=WRSS, emp=fits[[1]][[1]], modelFunc=mod
#     ,lower=c(0,.1,1,0.0001,1), upper=c(1,1,1000,.3,1000)
     ,lower=c(0,.2,0,0.0001,0), upper=c(.1,.999,7,.3,7)
     ,method="L-BFGS-B", control=list(trace=T))$par

#################################################################
# Model using a genetic algorithm to optimize, period 1
#################################################################

set.seed(321)
fits[[1]][[1]]$np = fits[[1]][[1]]$np/1000000
fit = ga(type="real-valued", fitness=nWRSS, emp=fits[[1]][[1]], modelFunc=mod
        ,min=c(0,0,0,0,0), max=c(1,1,7,1,7), maxiter=1000 )
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
# Model using a genetic algorithm to optimize, period 2
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
