

#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

#! READ IN DATA !#
library(foreign)
dat 		<- read.dta("C:/Users/lhund/Dropbox/PROJECTS/LQAS/ClusterPaper/minetti.dta")

### MAKE FIGURE DISPLAYING ICC FOR ALL AREAS ###

png("minetti.png")
plot(dat$coverage, dat$icc, xlab="Coverage", ylab=expression(rho), pch=16, cex.lab=1.5, cex.axis=1.5)
dev.off()
