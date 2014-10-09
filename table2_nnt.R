
library(xtable)

#! LOAD IN OUTPUT FROM EXAMPLE 2 - NNT !#
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/designs2.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/l2.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/e2.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/q2.Rda")

#! FILL IN BF DESIGNS IN SITUATION WHERE EQUIVALENT TO LQAS FOR EFFICIENCY !#

upper 	<- length(designs)
btemp 	<- list()
for(kk in 1:upper){
	btemp[[kk]] <- list(n = NA, rule = NA)
	if(kk!=4){
		pl		<- l[[kk]]$pl
		pu		<- l[[kk]]$pu
		sptemp 	<- designs[[kk]]$sp
		setemp 	<- designs[[kk]]$se
		ple 		<- setemp*pl + (1-sptemp)*(1-pl)
		pue 		<- setemp*pu + (1-sptemp)*(1-pu)
		dtemp 	<- l[[kk]]$rule
		ntemp 	<- l[[kk]]$n
		logbf	 	<- dbinom(dtemp, ntemp, pue, log = T) - dbinom(dtemp, ntemp, ple, log = T)
		bf 		<- exp(logbf)
		btemp[[kk]] <- list(n = ntemp, rule = bf)
	}
}
b 	<- btemp
rm(btemp)

# Plug in NAs for design 4 where rule does not exist
l[[4]] <- list(n = NA, rule=NA)

for(pp in 1:upper){
	if(length(q[[pp]]) > 1)
		q[[pp]] <- q[[pp]]$Q
}

# MAKE BF RULES #
logbf <- NULL
for(kk in 1:length(designs)){
	if(kk!=4){
		pl <- l[[kk]]$pl*l[[kk]]$se + (1-l[[kk]]$pl)*(1- l[[kk]]$sp)
		pu <- l[[kk]]$pu*l[[kk]]$se + (1-l[[kk]]$pu)*(1- l[[kk]]$sp)
		n  <- l[[kk]]$n
		d <-  l[[kk]]$rule 
		logbf[kk] <- dbinom(d, n, pu, log=T) - dbinom(d, n, pl, log=T)
	}
	if(kk==4) logbf[kk] <- NA
}
bf <- exp(logbf)

# Make matrix of results
mat 	<- NULL
for(kk in 1:upper){
	mat 	<- rbind(mat, c(q[[kk]], e[[kk]], 
		 	l[[kk]]$n, l[[kk]]$rule,bf[kk]))
}
mat 	<- as.data.frame(mat)
colnames(mat) 	<- c("Q", "alpha", "beta", "n_lqas", "d", "k")
xtable(mat, digits=c(0, 2, 2, 2, 0, 0, 2))