
#! LOAD IN OUTPUT FROM EXAMPLE 1 - VACCINATION !#
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/designs1.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/l1.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/e1.Rda")
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/q1.Rda")

btemp <- list()
for(kk in 1:length(designs)){
	filename <- paste("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/b1_",
			kk, ".Rda", sep="")
	load(filename)
	btemp[[kk]] <- b
}
b <- btemp
rm(btemp)

#FILL IN NAs FOR THE DESIGNS WITH NO RULE#
l[[8]] 	<- list(n=NA, rule=NA)
b[[8]]	<- list(n=NA, rule=NA)

for(pp in 1:length(designs)){
	if(length(q[[pp]]) > 1)
		q[[pp]] <- q[[pp]]$Q
}

library(xtable)
mat 	<- NULL
for(kk in 1:length(designs)){
	mat 	<- rbind(mat, 
			c(q[[kk]], e[[kk]], 
		 	l[[kk]]$n, l[[kk]]$rule, b[[kk]]$n, b[[kk]]$rule))
}
mat 	<- as.data.frame(mat)
colnames(mat) 	<- c("Q", "alpha", "beta", "n_lqas", "d", "n_bf", "k")
xtable(mat, digits=c(0, 2, 2, 2, 0, 0, 0, 2))
