
#! LOAD IN OUTPUT FROM EXAMPLE 2 - NNT !#
load("C:/Users/lhund/Dropbox/PROJECTS/LQAS/SensSpec/output/output_example2.Rda")

library(xtable)

# Plug in NAs for design 4 where rule does not exist
m[[4]]$rule <- NA
m[[4]]$n	<- NA
b[[4]]$rule <- NA
b[[4]]$n	<- NA

# Make matrix of results
mat 	<- NULL
for(kk in 1:length(e)){
	mat 	<- rbind(mat, c(e[[kk]], q[[kk]]$Q,
		 	m[[kk]]$n, m[[kk]]$rule, b[[kk]]$n, b[[kk]]$rule))
}
mat 	<- as.data.frame(mat)
colnames(mat) 	<- c("alpha", "beta", "Q", "n_lqas", "d", "n_bf", "k")
xtable(mat, digits=c(0, 2, 2, 2, 0, 0, 0, 2))