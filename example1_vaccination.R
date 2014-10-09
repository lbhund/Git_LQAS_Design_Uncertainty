

#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output")

#! READ IN DATA !#
library(foreign)
dat 		<- read.dta("C:/Users/lhund/Dropbox/PROJECTS/LQAS/ClusterPaper/minetti.dta")

library(VGAM)
library(lqasu)

### INPUT DESIGN PARAMETERS ###

nsim 	<- 50000
pl	<- .75
pu	<- .9
m 	<- 10
k 	<- 6
n	<- m*k
d 	<- 50
alpha	<- .1
beta	<- .1

### MAKE RHO DISTRIBUTIONS USING METHOD OF MOMENTS ###

mom 	<- function(rhos){
	xbar 	<- mean(rhos)
	sig 	<- var(rhos)
	a 	<- xbar*(xbar*(1-xbar)/sig - 1)
	b 	<- (1-xbar)*(xbar*(1-xbar)/sig - 1)
	p 	<- a/(a+b)
	r	<- 1/(a+b+1)
	out 		<- list()
	out$sims 	<- rbeta(nsim, a, b)
	out$pr	<- c(p, r)
	out$ab 	<- c(a,b)
	return(out)
}

rho 	<- function(sigma2, p){ sigma2/(p*(1-p))}

rhom		<- mom(dat$icc)
rhol		<- mom(dat$icc[dat$coverage < 85])
rhou  	<- mom(dat$icc[dat$coverage >= 85])

medianm 	<- qbeta(.5, shape1=rhom$ab[1], shape2=rhom$ab[2])
medianl 	<- qbeta(.5, shape1=rhol$ab[1], shape2=rhol$ab[2])
medianu 	<- qbeta(.5, shape1=rhou$ab[1], shape2=rhou$ab[2])

### CALCULATE DESIGNS ###

rholist <- list(
	list(rhol$pr, rhou$pr),
	list(medianl, medianu),
	list(rhom$pr, rhom$pr),
	medianm,
	list(rho(.1^2, pl), rho(.1^2, pu)),
	.1)

selist 	<- list(.77, c(.77,.02), c(.77,.005))
splist 	<- list(.84, c(.83, .01), c(.83, .005))

l 	<- list()
b 	<- list()
e 	<- list()
q 	<- list()

designs <- list()
for(kk in 1:length(rholist)){
	designs[[kk]]	<- list()
	designs[[kk]]$rho <- rholist[[kk]]
	designs[[kk]]$se 	<- 1
	designs[[kk]]$sp 	<- 1
}
for(kk in 1:length(selist)){
	jj 			<- kk + length(rholist)
	designs[[jj]]	<- list()
	designs[[jj]]$rho	<- rholist[[1]]
	designs[[jj]]$se 	<- selist[[kk]]
	designs[[jj]]$sp 	<- splist[[kk]]
}
save(designs, file="designs1.Rda")

#make lqas designs
set.seed(1)
for(kk in 1:length(designs)){
	if(kk!=8) 
		l[[kk]] 	<- lqasu( pl=pl, pu=pu, m = 10, alpha=alpha, beta=beta, nsim=nsim, 
					rho = designs[[kk]]$rho, se=designs[[kk]]$se, sp=designs[[kk]]$sp)
	if(kk==8)
		l[[kk]] 	<- "No rule exists"
}
save(l, file="l1.Rda")

#make q
set.seed(2)
for(kk in 1:length(designs)){
	q[[kk]] <- list(Q=0)
	if(length(designs[[kk]]$se) > 1 | length(designs[[kk]]$sp) > 1){
	q[[kk]] 	<- makeq( pl=pl, pu=pu, nsim=nsim, alpha=alpha, beta=beta,
				se=designs[[kk]]$se, sp=designs[[kk]]$sp)
	}
}
save(q, file="q1.Rda")

#make e
set.seed(3)
for(kk in 1:length(designs)){
	e[[kk]] 	<- errors( pl=pl, pu=pu, m = 10, n=n, d=d, nsim=nsim, 
				rho = designs[[kk]]$rho, se=designs[[kk]]$se, sp=designs[[kk]]$sp)
}
save(e, file="e1.Rda")


#make b
for(kk in 1:length(designs)){
	set.seed(4)
	if(kk!=8) 
		b	 	<- lqasu( pl=pl, pu=pu, m = 10, alpha=alpha, beta=beta, nsim=nsim, bfrule=T,
				rho = designs[[kk]]$rho, se=designs[[kk]]$se, sp=designs[[kk]]$sp)
	if(kk==8)
		b	 	<- "No rule exists"
	save(b, file=paste("b1_", kk, ".Rda", sep=""))
}





