
#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output")

#########################################
############# NNT EXAMPLE ###############
#########################################

library(lqasu)
library(VGAM)

### INPUT DESIGN PARAMETERS ###

nsim 	<- 100000
pl	<- .0005
pu	<- .003
n	<- 2540
d 	<- 2
alpha <- .1
beta 	<- .1

### CALCULATE DESIGNS ###

selist 	<- list(1, .7,   .7, .7, .7)
splist 	<- list(1,  1, .999, c(.999,.001), .99)

l 	<- list()
b 	<- list()
e 	<- list()
q 	<- list()


designs <- list()
for(kk in 1:length(selist)){
	designs[[kk]]	<- list()
	designs[[kk]]$rho <- 0
	designs[[kk]]$se 	<- selist[[kk]]
	designs[[kk]]$sp 	<- splist[[kk]]
}
save(designs, file="designs2.Rda")

#make lqas designs
set.seed(1)
for(kk in 1:length(designs)){
	rhotemp		<- NULL
	if(designs[[kk]]$rho > 0) 
		rhotemp 	<- designs[[kk]]$rho
	if(kk!= 4) l[[kk]] <- lqasu( pl=pl, pu=pu, alpha=alpha, beta=beta, m=10, nsim = nsim,
					se=designs[[kk]]$se, sp=designs[[kk]]$sp, rho=rhotemp)
	if(kk==4) l[[kk]]  <- NULL
}
save(l, file="l2.Rda")

#make q
set.seed(1)
for(kk in 1:length(designs)){
	rhotemp		<- NULL
	if(designs[[kk]]$rho > 0) 
		rhotemp 	<- designs[[kk]]$rho
	q[[kk]] <- 0
	if(is.null(rhotemp)==F | length(designs[[kk]]$se) > 1 | length(designs[[kk]]$sp) > 1)	
		q[[kk]] <- makeq(pl=pl, pu=pu, alpha=.1, beta=.1, nsim = nsim, m=10, n=n,
					se=designs[[kk]]$se, sp=designs[[kk]]$sp, rho=rhotemp)
}
save(q, file="q2.Rda")

#make e
set.seed(1)
for(kk in 1:length(designs)){
	rhotemp		<- NULL
	if(designs[[kk]]$rho > 0) 
		rhotemp 	<- designs[[kk]]$rho
	e[[kk]] 		<- errors(pl=pl, pu=pu, n=n, d=d, m=10, nsim = nsim,
					se=designs[[kk]]$se, sp=designs[[kk]]$sp, rho=rhotemp)
}
save(e, file="e2.Rda")
