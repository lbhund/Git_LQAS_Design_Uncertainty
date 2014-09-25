
library(lqasu)

pl	<- .0005
pu	<- .003
p 	<- seq(from=.995, to=1, by=.0001)
fp 	<- dbeta(p, .999*(1/.001-1), (1-.999)*(1/.001-1))
se	<- .7
sp	<- c(.999, .001)

file 	<- "nnt.png"
png(filename=file, width=480*2, height=480)
par(mfrow=c(1,2))
plot(p, fp, type="l", ylab="Density")
plotpe(pl, pu, se=se, sp=sp, addq=T, position="topright")
dev.off()

