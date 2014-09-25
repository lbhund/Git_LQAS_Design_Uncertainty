
#! CHANGE DIRECTORY TO WHEREEVER FILE SHOULD BE WRITTEN, e.g. C:\My Documents\Dropbox\... !#
setwd("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/pics")

load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/appendix_betabinsim.Rda")

file 	<- "bb_cdf.png"
png(filename=file, width=480*2, height=480)
#tiff(filename=file, width=8, height=4, units="in", res=300, compression="lzw")
par(mfrow=c(1,2))
plot(out[,3], out[,5], xlab="p", ylab="Maximum distance between cdfs")
plot(out[,4], out[,5], xlab=expression(rho), ylab="Maximum distance between cdfs")
dev.off()

