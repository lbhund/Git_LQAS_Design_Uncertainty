
# CHECK WHETHER BF AND LQAS RULES ARE THE SAME FOR VACCINE EXAMPLE #

load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/l1.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/designs1.Rda")

btemp <- list()
for(kk in 1:length(designs)){
	filename <- paste("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/b1_",
			kk, ".Rda", sep="")
	load(filename)
	btemp[[kk]] <- b
}
b <- btemp
rm(btemp)

out <- rep(F, length(designs))
for(kk in 1:length(designs)){
if(kk!=8) out[kk] <- all.equal(b[[kk]]$drules[,3], 
		         b[[kk]]$drules[order(-b[[kk]]$drules[,3]),3])
}


out