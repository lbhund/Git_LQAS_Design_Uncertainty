

#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

#pl, pu are lower and upper thresholds
#rhol, rhou are lower and upper icc - either constant, or beta prior parameters, parameterized in terms of mean and overdispersion
#n is total sample size, m is number per cluster, k is number of clusters; must specify 2 of the 3
#nsim is number of simulations to base output on


library(lqasu)
nsim 	<- 100000
pl 	<- .75
pu 	<- .9
rholist 	<- list(.01, .1)
n 	<- 60
m 	<- 10
k 	<- n/m

pe 	<- makepe(pl=pl, pu=pu, rho=rholist, n=n, m=m, nsim = nsim)
q 	<- makeq(pl=pl, pu=pu, rho=rholist, n=n, m=m, nsim = nsim) 

bf 	<- function(x) {log((mean(sapply(1:nsim, function(tt) {dbinom(x, n, pe$pue[tt])}))))-
			    log((mean(sapply(1:nsim, function(tt) {dbinom(x, n, pe$ple[tt])}))))}
x 	<- 0:n
bfout	<- sapply(x, bf)

#output is graph of ple and pue, as well as bayes factor as a function of x
png("bffig.png", height=480, width=2*480)
par(mfrow=c(1,2))
plot(density(pe$ple), xlim=c(.6,1), lty=2, lwd=2, cex=1.7, 
	cex.lab=1.7, cex.axis=1.7, xlab="p", ylab="Density", main="")
lines(density(pe$pue), lwd=2)
plot(x, bfout, type="l", ylab="Log-likelihood ratio", xlab="X", cex=1.7, 
			cex.lab=1.7, cex.axis=1.7, cex.main=1.9)
dev.off()
