beta = 0.1
toPlot = locoutNeighbor(dat, X, Y, propneighb=beta, usemax=1, npoints=50, variant="dist")
plot(c(0,50), c(0,1), type="n")
for(i in 1:100)
  lines(1:50, toPlot[i,-1], col=gray(0.7))
lines(1:50,rep(beta,50),col=2)

beta = .3
toPlot = locoutNeighbor(matrix(rnorm(200),ncol=2), X=rnorm(100), Y=rnorm(100), propneighb=beta
        ,usemax=1, npoints=50, variant="dist" )
plot(c(0,50), c(0,1), type="n")
for(i in 1:100)
  lines(1:50, toPlot[i,-1], col=gray(0.7))
lines(1:50,rep(beta,50),col=2)