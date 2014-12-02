if(Sys.info()[4]=="jb"){
  setwd("/media/storage/Professional Files/Mines/SmartGeo/Queens/")
  source("/media/storage/Github/SmartGeo/simulate_data_for_real.R")
}
if(Sys.info()[4]=="JOSH_LAPTOP"){
  setwd("~/Professional Files/Mines/SmartGeo/Queens/")
  source("~/Github/SmartGeo/simulate_data_for_real.R")
}
if(grepl("ch120", Sys.info()[4])){
  setwd("~/SmartGeo/Queens")
  source("~/Github/SmartGeo/simulate_data_for_real.R")
}

load("Data/ground_for_nnet.RData")
load("Data/station_pairs.RData")
load("Data/station_new_coords.RData")
load("Results/nnet_model_all_data.RData")
load("Results/fitted_sp_variograms_bessel.RData")

dir.create("Results/Simulated Data")

nnetMod = fit
nnetMod$coefnames
head(ground,2)
covariate = ground[,c("nnetA", "nnetBC", "nnetD", "nnetYL", "nnetEast2", "nnetNorth2", "nnetTime")]

#Generate a realization from the model
sim = simulate(fits, nnetMod, station_pairs, covariate, ground$StationID, ground$Time
        ,randErrorPct=0, sysErrorThresh=Inf)
sim$Error = NULL
sim$model = NULL
sim$randErrorFl = NULL
sim$sysError = NULL
sim = cbind(sim, covariate)

#Fit a neural network to the simulated data, store the errors
fit = nnet( form=Value ~ nnetA + nnetBC + nnetD + nnetYL + nnetEast2 + nnetNorth2 + nnetTime,
    data=sim, size=10, linout=TRUE, maxit=100000)
sim$Prediction = predict(fit, newdata=sim)
oldCov = sim[,grepl("nnet",colnames(sim))]
sim$oldPrediction = predict(nnetMod, newdata=sim)
ggsave("Results/Simulated Data/nnet_model_comparison.png",
  ggplot(sim, aes(x=Prediction, y=oldPrediction) ) + geom_point(alpha=.1) +
    labs(x="New Prediction", y="Old Prediction")
)
sim$Error = sim$Value - sim$Prediction
sim$Error = sim$Error[,1]
sim$Prediction = NULL
sim$oldPrediction = NULL
sim$A = NULL
sim$BC = NULL
sim$D = NULL
sim$YL = NULL
sim$nnetTime = NULL
sim$nnetA = NULL
sim$nnetBC = NULL
sim$nnetD = NULL
sim$nnetYL = NULL
sim$nnetEast2 = NULL
sim$nnetNorth2 = NULL
gc()

save(sim, file="Results/SIM_simulated_dataset_with_nnet.RData")

#Identify outliers in a very rudimentary way (just to remove influence on variograms)
station = group_by( sim, StationID )
sds = summarize( station, mu=huber(Value)[[1]], s=huber(Value)[[2]] )
#ggplot( sds, aes(x=StationID, y=s ) ) + geom_point()
sim = merge(sim, sds)
sim$outlier = F
sim$stat = (sim$Value-sim$mu)/sim$s
sim$outlier[abs(sim$stat)>6] = T
sim$mu = NULL
sim$s = NULL
sim$stat = NULL

#Compute the empirical variogram
prd1 = 1:800
prd2 = 801:4500
prd3 = 4501:6105
fits = list()
times = sort(unique(sim$Time))
for( prd in list(prd1, prd2, prd3) ){
#for( prd in list(prd2, prd3) ){
  fit = empVario(data=sim[sim$Time %in% times[prd] & !sim$outlier,]
    ,tlags=0:30, alpha=c(0,45,90,135), varName="Error", boundaries=0:65*10, cutoff=650)
  fits[[length(fits)+1]] = fit
  print(plot.empVario(fit, boundaries=0:65*10, model=NULL, adj=TRUE, rmAni=T))
  save(fits, prd1, prd2, prd3, file="Results/Final Variograms/SIM_variograms_error_nocontam.RData")
}

#Fit the empirical variogram with a bessel model
empTemp = fits[[1]][[1]]
empTemp = empTemp[empTemp$dist<300,]
set.seed(123)
fit = ga(type="real-valued", fitness=nWRSS.st, emp=empTemp, modelFunc=modBesExp
      ,min=c(1e-6,0,.0001,1e-4,0.0001,0), max=c(.1,1,1/3,10,.3,7), maxiter=1000 )
fits[[1]][[2]] = slot(fit, "solution")
toPlot = fits[[1]]
toPlot[[1]] = toPlot[[1]][toPlot[[1]]$dist<=500,]
png("Results/Final Variograms/SIM_bessel_exponential_model_GA_prd1.png", width=800, height=800)
plot.empVario(toPlot, model=modBesExp, boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

empTemp = fits[[2]][[1]]
empTemp = empTemp[empTemp$dist<300,]
set.seed(123)
fit = ga(type="real-valued", fitness=nWRSS.st, emp=empTemp, modelFunc=modBesExp
      ,min=c(1e-6,0,.0001,1e-4,0.0001,0), max=c(.1,1,1/3,10,.3,7), maxiter=1000 )
fits[[2]][[2]] = slot(fit, "solution")
toPlot = fits[[2]]
toPlot[[1]] = toPlot[[1]][toPlot[[1]]$dist<=500,]
png("Results/Final Variograms/SIM_bessel_exponential_model_GA_prd2.png", width=800, height=800)
plot.empVario(toPlot, model=modBesExp, boundaries=0:65*10, adj=T, rmAni=T)
dev.off()

empTemp = fits[[3]][[1]]
empTemp = empTemp[empTemp$dist<300,]
set.seed(123)
fit = ga(type="real-valued", fitness=nWRSS.st, emp=empTemp, modelFunc=modBesExp
      ,min=c(1e-6,0,.0001,1e-4,0.0001,0), max=c(.1,1,1/3,10,.3,7), maxiter=1000 )
fits[[3]][[2]] = slot(fit, "solution")
toPlot = fits[[3]]
toPlot[[1]] = toPlot[[1]][toPlot[[1]]$dist<=500,]
png("Results/Final Variograms/SIM_bessel_exponential_model_GA_prd3.png", width=800, height=800)
plot.empVario(toPlot, model=modBesExp, boundaries=0:65*10, adj=T, rmAni=T)
dev.off()
