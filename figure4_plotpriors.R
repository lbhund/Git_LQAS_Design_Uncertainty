
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

### VACCINATION COVERAGE ###
pvec 	<- seq(0, 1, .0001)

ab 	<- function(vec){
	mu 	<- vec[1]
	psi 	<- vec[2]
	a 	<- mu*(1/psi - 1)
	b 	<- (1-mu)*(1/psi - 1)
	return(c(a, b))
}

rho	=list(c(.059, .047), c(.021, .009), c(.033, .039))
se1	=c(.77, .02)
sp1	=c(.83, .01)
se2	=c(.77, .005)
sp2	=c(.83, .005)
m 	= 10
k 	= 6
	
sedist1 	<- dbeta(pvec, shape1=ab(se1)[1], shape2 = ab(se1)[2])
spdist1 	<- dbeta(pvec, shape1=ab(sp1)[1], shape2 = ab(sp1)[2])
sedist2 	<- dbeta(pvec, shape1=ab(se2)[1], shape2 = ab(se2)[2])
spdist2 	<- dbeta(pvec, shape1=ab(sp2)[1], shape2 = ab(sp2)[2])
rholdist	<- dbeta(pvec, shape1=ab(rho[[1]])[1], shape2 = ab(rho[[1]])[2])
rhoudist 	<- dbeta(pvec, shape1=ab(rho[[2]])[1], shape2 = ab(rho[[2]])[2])
rhodist 	<- dbeta(pvec, shape1=ab(rho[[3]])[1], shape2 = ab(rho[[3]])[2])


png("priors_v.png", width=480*2, height=480*2)
par(mfrow=c(2,2),mar=c(5, 3, 1, 1))

plot(pvec, sedist2, xlim=c(.55,1), xlab=expression(s[e]), ylab="", type="l", lwd=2, cex.lab=2.5, cex.axis=1.5)
lines(pvec, sedist1, lty=2, lwd=2)
legend("topright", lty=c(2,1), c("Beta(.77, .02)", "Beta(.77, .005)"), lwd=2, cex=1.5)

plot(pvec, spdist2, xlim=c(.7,1), xlab=expression(s[p]), ylab="", type="l", lwd=2, cex.lab=2.5, cex.axis=1.5)
lines(pvec, spdist1, lty=2, lwd=2)
legend("topright", lty=c(2,1), c("Beta(.83, .01)", "Beta(.83, .005)"), lwd=2, cex=1.5)

plot(pvec, rhodist, xlim=c(0,.3), xlab=expression(rho[l]), ylab="", type="l", lwd=2, cex.lab=2.5, cex.axis=1.5)
lines(pvec, rholdist, lty=2, lwd=2)
legend("topright", lty=c(2,1), c("Beta(.059, .047)", "Beta(.033, .039)"), lwd=2, cex=1.5)

plot(pvec, rhodist, xlim=c(0,.4), xlab=expression(rho[u]), ylab="", type="l", lwd=2, cex.lab=2.5, cex.axis=1.5)
lines(pvec, rhoudist, lty=2, lwd=2)
legend("topright", lty=c(2,1), c("Beta(.021, .009)", "Beta(.033, .039)"), lwd=2, cex=1.5)
dev.off()

### NNT ###

sp	=c(.999, .001)
spdist 	<- dbeta(pvec, shape1=ab(sp)[1], shape2 = ab(sp)[2])
png("priors_n.png", width=480, height=480)
par(mar=c(5, 3, 1, 1))
plot(pvec, spdist, xlim=c(.99,1), xlab=expression(s[p]), ylab="", type="l", lwd=2, cex.lab=2.5, cex.axis=1.5)
legend("topleft", lty=1, c("Beta(.999, .001)"), lwd=2, cex=1.5)
dev.off()
