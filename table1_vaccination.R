
#! LOAD IN OUTPUT FROM EXAMPLE 1 - VACCINATION !#
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/output_example1.Rda")

library(xtable)
mat 	<- NULL
for(kk in 1:length(e)){
	mat 	<- rbind(mat, c(e[[kk]], q[[kk]]$Q,
		 	m[[kk]]$n, m[[kk]]$rule, b[[kk]]$n, b[[kk]]$rule))
}
mat 	<- as.data.frame(mat)
colnames(mat) 	<- c("alpha", "beta", "Q", "n_lqas", "d", "n_bf", "k")
xtable(mat, digits=c(0, 2, 2, 2, 0, 0, 0, 2))
