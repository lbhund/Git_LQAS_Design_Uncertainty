
#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

set.seed(3)
library(lqasu)

pl 	<- .75
pu 	<- .9
alpha <- .1
beta 	<- .1
se 	<- c(.77, .005)
sp 	<- c(.83, .005)
nsim 	<- 50000

q <- makeq(pl=pl, pu=pu, se=se, sp=sp, alpha=alpha, beta=beta, nsim=nsim)

pedist	<- makepe(pl=pl, pu=pu, se=se, sp=sp, nsim=nsim)
xlim 		<- as.numeric(range(unlist(pedist)))
ylim 		<- range(c(density(pedist$ple)$y,density(pedist$pue)$y))

pvec		<- seq(0, 1, by=min(.005, pl/1000))
plx		<- sapply(pvec, function(pp) length(which(pedist$ple < pp))/nsim)
pux		<- sapply(pvec, function(pp) length(which(pedist$pue < pp))/nsim)

plset 	<- which(1-plx <= beta)
puset 	<- which(pux <= alpha)

plemin 	<- range(pvec[plset])[1]
puemax 	<- range(pvec[puset])[2]


amin 		<- min(pux[which(1 - plx <= beta)])
bmin 		<- min(1 - plx[which(pux <= alpha)])
aline 	<- pvec[which(pux==amin)]
bline 	<- pvec[which(1-plx==bmin)] 


file 	<- "q_illustration.png"
png(filename=file, width=480*2, height=480)
par(mar=c(5, 5, 2, 2) + 0.1)
plot(density(pedist$ple), xlim=c(.55, .8), ylim=ylim, main="", xlab="p", cex.axis=1.9, cex.lab=1.9)
lines(density(pedist$pue))
#abline(v=puemax, lty=2)
#abline(v=plemin, lty=2)
abline(v=aline, lty=2)
abline(v=bline, lty=2)
abline(v=q$p)
#legend("topright", paste("Q=",round(q$Q, digits=3)), cex=1.9)
text(.62, 10, expression(f[l]), cex=1.9)
text(.72, 10, expression(f[u]), cex=1.9)
dev.off()

q$p
q$Q
aline
q$alphamin
bline
q$betamin



