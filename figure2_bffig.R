

#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

#pl, pu are lower and upper thresholds
#rhol, rhou are lower and upper icc - either constant, or beta prior parameters, parameterized in terms of mean and overdispersion
#n is total sample size, m is number per cluster, k is number of clusters; must specify 2 of the 3
#nsim is number of simulations to base output on


library(VGAM)
plotQ <- function(pl=NULL, pu=NULL, rhol=NULL, rhou=NULL, n=NULL, m=NULL, k=NULL, nsim=1000) {

	if(is.null(pl) == T | is.null(pu) == T | is.null(rhol) == T | is.null(rhou) == T) 
		stop("Must specify ps and rhos")
	if(is.null(n) == T) n <- m*k
	if(is.null(k) == T) k <- n/m
	if(is.null(m) == T) m <- n/k

	psi	<- function(rho) (m-1)*rho/(m*k-1)

	par(mfrow=c(1,2), mar=c(5, 5, 2, 1) + 0.1)
	
	pvec 	<- seq(0, 1, length=100)

	#calculate Q and pcut
	if(length(rhol) == 1 & length(rhou) == 1){
		fixrho	<- TRUE
		psil 	<- psi(rhol)
		psiu 	<- psi(rhou)

		plx 	<- pbeta(pvec, pl*(1/psil - 1), (1-pl)*(1/psil - 1))
		pux 	<- pbeta(pvec, pu*(1/psiu - 1), (1-pu)*(1/psiu - 1))

		flx 	<- dbeta(pvec, pl*(1/psil - 1), (1-pl)*(1/psil - 1))
		fux 	<- dbeta(pvec, pu*(1/psiu - 1), (1-pu)*(1/psiu - 1))
	}
	if(length(rhol) > 1 & length(rhou) > 1){
		fixrho 	<- FALSE
		rhotempl	<- rbeta(nsim, rhol[1]*(1/rhol[2] -1), (1-rhol[1])*(1/rhol[2] -1))
		rhotempu	<- rbeta(nsim, rhou[1]*(1/rhou[2] -1), (1-rhou[1])*(1/rhou[2] -1))
		psil 	<- psi(rhotempl)
		psiu 	<- psi(rhotempu)

		rlx 	<- rbeta(nsim, pl*(1/psil - 1), (1-pl)*(1/psil - 1))
		rux 	<- rbeta(nsim, pu*(1/psiu - 1), (1-pu)*(1/psiu - 1))

		plx	<- sapply(pvec, function(pp) length(which(rlx < pp))/nsim)
		pux	<- sapply(pvec, function(pp) length(which(rux < pp))/nsim)
	}
	Qvec 	<- 1-plx + pux
	Q	<-round(min(Qvec), digits=3)
	pcut 	<- round(pvec[which(Qvec == min(Qvec))], digits=3)
	#plot pdist
	if(fixrho==T){
		plot(pvec, flx, lwd=2, lty=2, type="l", cex=1.7, cex.main=1.9, 
			cex.lab=1.7, cex.axis=1.7, xlab="p", ylab="Density")
		lines(pvec, fux, lwd=2)
	} else{
		plot(density(rlx), xlim=c(0,1), lty=2, lwd=2, cex=1.7, cex.main=1.9, 
			cex.lab=1.7, cex.axis=1.7, xlab="p", ylab="Density")
		lines(density(rux), lwd=2)
	}
	abline(v=pcut)
	text(.8, 20, bquote(p ==.(pcut)), pos=4, cex=1.5)
	text(.8, 18, bquote(Q==.(Q)), pos=4, cex=1.5)

	#plot BF vs x
	if(fixrho==T)
		bf 	<- function(x) (dbetabinom(x, n, pu, psiu, log=T) - dbetabinom(x, n, pl, psil, log=T))
	if(fixrho==F)
		bf 	<- function(x) mean(sapply(1:nsim, function(tt) (dbetabinom(x, n, pu, psiu[tt], log=T) - dbetabinom(x, n, pl, psil[tt], log=T))))
	x 	<- 0:n
	plot(x, sapply(x, bf), type="l", ylab="BF", xlab="X", cex=1.7, 
			cex.lab=1.7, cex.axis=1.7, cex.main=1.9)
}

#output is graph of ple and pue, as well as bayes factor as a function of x
png("bffig.png", height=480, width=2*480)
plotQ(pl=.75, pu=.9, rhol=.01, rhou=.1, n=60, m=10)
dev.off()
