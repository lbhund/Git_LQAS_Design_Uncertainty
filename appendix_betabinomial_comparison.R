
#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output")

library(VGAM)
set.seed(4)

### CONVOLUTION FUNCTION ###
makedensity <- function(k, m, p, rho){
	K 	<-  rep(m, k)
	d1 	<-  dbetabinom(0:K[1],K[1],p,rho)
	for(i in 2:length(K)){ 
    		d2 <- c( dbetabinom(0:K[i],K[i],p,rho), rep(0,length(d1)-1) )
    		d1 <- c(d1,rep(0,K[i]))
    		zz <- Re( fft( fft(d1)*fft(d2), inverse=TRUE))/length(d1)
    		d1 <- zz
	}	
	ax 	<- 0:(length(d1)-1)
	return(cbind(ax, d1))
}  


### INPUT PARAMETERS ###

out 	<- matrix(NA, nrow=4680, ncol=5)
ii 	<- 0
for(k in 5:10){
for(m in 5:10){
for(p in c(.01, .05, seq(.1, .9, by=.1), .95, .99)){
for(rho in seq(from=.01, to = .1, by = .01)){
	ii 	<- ii + 1
	n 	<- k*m

psi 		<- (m-1)*rho/(m*k-1)

cdf1 	<- NULL
cdf2 	<- NULL
for(cut in 0:n){
	cdf1[cut +1] 	<- sum(makedensity(k, m, p, rho)[1:(cut+1),2])
	cdf2[cut +1] 	<- pbetabinom(cut, n, p, psi)
}

out[ii,] <- c(k, m, p, rho, max(abs(cdf1 - cdf2)))

}}}}

save(out, file="appendix_betabinsim.Rda")
