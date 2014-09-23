

### MAKE EFFECTIVE COVERAGE DISTRIBUTIONS ###
makepe <- function(pl, pu, n=NULL, m=NULL, k=NULL, rho=NULL, se=NULL, sp=NULL, nsim=5000){
	### PARAMETER DEFINITIONS ###
	if(is.null(m) == T) m <- n/k
	if(is.null(k) == T) k <- n/m
	if(is.null(n) == T) n <- m*k
	if(is.null(se) == T & is.null(sp) == T){
		ple 	<- pl	
		pue 	<- pu
	}
	if(is.null(se) == T & is.null(sp) == F){
		se 	<- 1
	}
	if(is.null(se) == F & is.null(sp) == T){
		sp 	<- 1
	}
	if(is.null(rho) == F){if(is.list(rho) == F){if(length(rho) == 1){
		if(rho <= 0 | rho >= 1) stop("rho must be between 0 and 1")}}}

	### MAKE DISTRIBUTIONS ###
	if(is.null(se) == F | is.null(sp) == F){
		if(length(se) == 1 & length(sp) == 1){
			ple 	<- se*pl + (1-sp)*(1-pl)
			pue 	<- se*pu + (1-sp)*(1-pu)
		}
		if(length(se) == 2 & length(sp) == 2){
			sea 	<- se[1]*(1/se[2] - 1)
			seb 	<- (1-se[1])*(1/se[2] - 1)
			spa 	<- sp[1]*(1/sp[2] - 1)
			spb 	<- (1-sp[1])*(1/sp[2] - 1)
			sesamp 	<- rbeta(nsim, sea, seb)
			spsamp	<- rbeta(nsim, spa, spb)
			ple		<- sesamp*pl + (1-spsamp)*(1-pl)
			pue		<- sesamp*pu + (1-spsamp)*(1-pu)
		}
		if(length(se) == 1 & length(sp) == 2){
			spa 	<- sp[1]*(1/sp[2] - 1)
			spb 	<- (1-sp[1])*(1/sp[2] - 1)
			spsamp	<- rbeta(nsim, spa, spb)
			ple		<- se*pl + (1-spsamp)*(1-pl)
			pue		<- se*pu + (1-spsamp)*(1-pu)
		}
		if(length(se) == 2 & length(sp) == 1){
			sea 	<- se[1]*(1/se[2] - 1)
			seb 	<- (1-se[1])*(1/se[2] - 1)
			sesamp 	<- rbeta(nsim, sea, seb)
			ple		<- sesamp*pl + (1-sp)*(1-pl)
			pue		<- sesamp*pu + (1-sp)*(1-pu)
		}
		if(length(se) > 2 | length(sp) > 2)
			stop("Se and Sp not correctly specified")
	} 
	if(is.null(rho)==F & is.list(rho) == F){
		if(length(rho)==1)
			psil	<- psiu	<- (m-1)*rho/(m*k-1)
		if(length(rho)>1){
			psil	<- (m-1)*rho[1]/(m*k-1)
			psiu	<- (m-1)*rho[2]/(m*k-1)
		}
		plea		<- ple*(1/psil - 1)
		pleb		<- (1-ple)*(1/psil - 1)
		puea		<- pue*(1/psiu - 1)
		pueb		<- (1-pue)*(1/psiu - 1)
		pledist 	<- rbeta(nsim, plea, pleb)
		puedist 	<- rbeta(nsim, puea, pueb)
	}
	if(is.list(rho)==T){
		rhola <- rho[[1]][1]*(1/rho[[1]][2] - 1)
		rholb <- (1-rho[[1]][1])*(1/rho[[1]][2] - 1)
		rhoua <- rho[[2]][1]*(1/rho[[2]][2] - 1)
		rhoub <- (1-rho[[2]][1])*(1/rho[[2]][2] - 1)
		rhol 	<- rbeta(nsim, rhola, rholb)
		rhou 	<- rbeta(nsim, rhoua, rhoub)
		psil	<- (m-1)*rhol/(m*k-1)
		psiu	<- (m-1)*rhou/(m*k-1)
		plea		<- ple*(1/psil - 1)
		pleb		<- (1-ple)*(1/psil - 1)
		puea		<- pue*(1/psiu - 1)
		pueb		<- (1-pue)*(1/psiu - 1)
		pledist 	<- rbeta(nsim, plea, pleb)
		puedist 	<- rbeta(nsim, puea, pueb)
	}
	if(is.null(rho)==T){
		pledist 	<- ple
		puedist 	<- pue
	}
	out <- list()
	out$ple <- pledist
	out$pue <- puedist
	return(out)
}

### CALCULATE ERRORS ###
errors <- function(pl, pu, n=NULL, d, m=NULL, k=NULL, rho=NULL, se=NULL, sp=NULL, nsim=5000){
	if(is.null(n) == T) n <- m*k
	pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)
	lower_vec 	<- rbinom(nsim, n, pedist$ple)
	upper_vec 	<- rbinom(nsim, n, pedist$pue)
	alpha	<- length(which(upper_vec <= d))/nsim
	beta  <- length(which(lower_vec > d))/nsim
	return(c(alpha, beta))
}

### CALCULATE LQAS RULE ###
lqasu <- function(pl, pu, rho=NULL, se=NULL, sp=NULL, alpha = 0.1, beta = 0.1, m = NULL, 
			k = NULL, miter = F, nsim = 3000, add = 1){

	#define pe dist first if independent of n for computational efficiency
	if(is.null(rho) == T)
		pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)

	# fix clusters or samples within cluster
    	stop 	<- F
    	fixk 	<- TRUE
    	if (is.null(k) == T) 
    	    fixk 	<- FALSE
    	if (is.null(m) == T & fixk == FALSE) 
		m <- 1
    	#	stop("Must specify m or k") 
    	if (fixk == TRUE) 
    	    m 	<- 1
    	if (fixk == FALSE) 
    	    k 	<- 1
 	while (stop == F) {
		if (fixk == TRUE) 
    	            m 	<- m + add
    	      if (fixk == FALSE) 
    	            k 	<- k + add
    	      n 	<- m * k
		if(is.null(rho) == F)
			pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)
    	      d 	<- (round(0.5 * min(mean(pedist$ple)) * n)):(round(2 * max(mean(pedist$pue)) * n))


		lower_vec 	<- rbinom(nsim, n, pedist$ple)
		upper_vec 	<- rbinom(nsim, n, pedist$pue)
		alpha_vec	<- sapply(d, function(dd) {
			return(length(which(upper_vec <= dd))/nsim)})
		beta_vec	<- sapply(d, function(dd) {
			return(length(which(lower_vec > dd))/nsim)})

    	      ind 	<- rep(0, length(d))
      	ind[alpha_vec <= alpha & beta_vec <= beta] <- 1
      	if (sum(ind) > 0) 
      		stop 	<- T
      	rule 	<- d[which(ind == 1)]
      }
	out 		<- NULL
      out$rule 	<- rule
      out$n	 	<- n
      out$alpha 	<- alpha_vec[which(ind == 1)]
      out$beta 	<- beta_vec[which(ind == 1)]
      out$alpha_max 	<- alpha
      out$beta_max 	<- beta
      out$pl 	<- pl
      out$pu 	<- pu
	out$imperfect <- F
	if(is.null(se) == F | is.null(sp) == F){
		out$imperfect <- T
		out$sp	<- sp
		out$se	<- se
	}	
      out$m 	<- m
      out$k 	<- k
      out$rho 	<- rho
      out$family 	<- "cluster"
      class(out) 	<- "lqasu"
     	return(out)
}


### SUMMARIZE LQAS DESIGN ###
summary.lqasu <-function (object) {

	if(length(object$rule) > 1) 
		cat("NOTE: More than one decision rule found. \n")

    	cat("The LQAS design parameters are: \n\n")
    	cat("Sample Size:", object$n, "\n")
    	if(length(object$rule) ==1)
  		cat("Decision Rule:", object$rule, "\n\n")
    	if(length(object$rule) > 1)
		cat("Decision Rules:", object$rule, "\n\n")

    	cat("Classify as high if X > d.", "\n\n", sep = "")
    	cat("The true error levels for this design are: \n alpha =", 
        	round(object$alpha, digits = 4), "\n  beta =", round(object$beta, digits = 4), "\n\n")

    	cat("Design Parameters:", "\n")
    	cat("pl = ", object$pl, "and pu = ", object$pu, "\n")
    	cat("alpha = ", object$alpha_max, "and beta = ", object$beta_max, "\n\n")

	if(object$imperfect == T){
		cat("Imperfect instrument with mean sensitivity ", object$se[1], " and specificity ", object$sp[1], ".\n", sep="")
		if(length(object$se > 1) | length(object$sp) > 1) cat("Sensitivity and specificity assumed unknown. \n")
		#("Revised mean thresholds are pl = ", round(object$ple, digits=4), " and pu = ", round(object$pue, digits=4), ".\n\n", sep="")
	}

	if(object$family=="cluster"){
		cat("Probabilities were calculated assuming cluster sampling was used.  \n")
		if(length(object$rho)==1)
			cat("Sample ", object$k, " clusters and ", object$m, " individuals per cluster, assuming an ICC of ", object$rho, ".\n", sep="")
		if(length(object$rho)==2 & is.list(object$rho) == F)
			cat("Sample ", object$k, " clusters and ", object$m, " individuals per cluster, assuming an ICC of \n", 
			object$rho[1], " for pl = ", object$pl, " and ", object$rho[2], " for pu = ", object$pu, ".\n", sep="")
		if(length(object$rho)==2 & is.list(object$rho) == T)
			cat("Uncertainty in ICC incorporated by design. \n
			Sample ", object$k, " clusters and ", object$m, " individuals per cluster, assuming a mean ICC of \n", 
			round(object$rho[[1]][1], digits=4), " for pl = ", object$pl, " and ", round(object$rho[[2]][1], digits=4), " for pu = ", object$pu, ".\n", sep="")

	}
}


### MAKE BF RULE ###

bfrule <- function(pl, pu, rho=NULL, se=NULL, sp=NULL, alpha = 0.1, beta = 0.1, m = NULL, 
			k = NULL, miter = F, nsim = 3000, add = 1){

	#define pe dist first if independent of n for computational efficiency
	if(is.null(rho) == T)
		pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)

	# fix clusters or samples within cluster
    	stop 	<- F
    	fixk 	<- TRUE
    	if (is.null(k) == T) 
    	    fixk 	<- FALSE
    	if (is.null(m) == T & fixk == FALSE) 
		m <- 1
    	#	stop("Must specify m or k") 
    	if (fixk == TRUE) 
    	    m 	<- 1
    	if (fixk == FALSE) 
    	    k 	<- 1
 	while (stop == F) {
		if (fixk == TRUE) 
    	            m 	<- m + add
    	      if (fixk == FALSE) 
    	            k 	<- k + add
    	      n 	<- m * k
		if(is.null(rho) == F)
			pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)
    	      d 	<- (round(0.5 * min(mean(pedist$ple)) * n)):(round(2 * max(mean(pedist$pue)) * n))
		d 	<- d[which(d <= n)]

		bf 		<- function(x) mean(sapply(1:nsim, function(tt) 
			(dbinom(x, n, pedist$pue[tt], log=T) - dbinom(x, n, pedist$ple[tt], log=T))))
		xu_dist	<- table(rbinom(nsim, n, pedist$pue))/nsim
		xl_dist	<- table(rbinom(nsim, n, pedist$ple))/nsim
		bf_u 		<- sapply(as.numeric(row.names(xu_dist)), function(tt) bf(tt))
		bf_l		<- sapply(as.numeric(row.names(xl_dist)), function(tt) bf(tt))	

		alpha_vec 	<- sapply(d, function(dd) sum(0, xu_dist[which(bf_u <= bf(dd))]))
		beta_vec  	<- sapply(d, function(dd) sum(0, xl_dist[which(bf_l >  bf(dd))]))

	    	ind 		<- rep(0, length(d))
	      ind[alpha_vec <= alpha & beta_vec <= beta] <- 1
	      if (sum(ind) > 0){
      		stop 	<- T
	      	rule 		<- d[which(ind == 1)]
			bfrule	<- bf(rule)
			drulesmat 	<- t(sapply(0:n, function(tt) c(tt, aa <- bf(tt), ifelse(aa <= bfrule, 1, 0))))
			drules	<- drulesmat[which(drulesmat[,2]==1),1]
		}
	}
	out 		<- NULL
      out$rule 	<- exp(bf(rule))
	out$drules	<- drulesmat
      out$n	 	<- n
      out$alpha 	<- alpha_vec[which(ind == 1)]
      out$beta 	<- beta_vec[which(ind == 1)]
      out$alpha_max 	<- alpha
      out$beta_max 	<- beta
      out$pl 	<- pl
      out$pu 	<- pu
	out$imperfect <- F
	if(is.null(se) == F | is.null(sp) == F){
		out$imperfect <- T
		out$sp	<- sp
		out$se	<- se
	}	
      out$m 	<- m
      out$k 	<- k
      out$rho 	<- rho
      out$family 	<- "cluster"
      class(out) 	<- "lqasu"
     	return(out)
}



### MAKE Q ###

makeq <- function(pl, pu, n=NULL, m=NULL, k=NULL, rho=NULL, se=NULL, sp=NULL, nsim=5000){
	
	if(is.null(rho) == T & length(sp) == 1 & length(se)==1) 
		stop("No design uncertainty - must specify rho or measurement error prior")
	if(is.null(n) == T) n <- m*k
	pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)
	pvec	<- seq(0, 1, by=.005)
	plx	<- sapply(pvec, function(pp) length(which(pedist$ple < pp))/nsim)
	pux	<- sapply(pvec, function(pp) length(which(pedist$pue < pp))/nsim)

	Qvec 	<- 1-plx + pux
	Q	<- min(Qvec)
	pcut 	<- pvec[which(Qvec == min(Qvec))]
	out 	<- list()
	out$Q	<- Q
	out$p	<- pcut
	return(out)
}


### PLOT EFFECTIVE COVERAGE ###

plotpe <- function(pl, pu, n=NULL, m=NULL, k=NULL, rho=NULL, se=NULL, sp=NULL, nsim=5000,addq=F,...){
	pedist	<- makepe(pl=pl, pu=pu, n=n, m=m, k=k, rho=rho, se=se, sp=sp, nsim=nsim)
	xlim <- as.numeric(range(unlist(pedist)))
	ylim <- range(c(density(pedist$ple)$y,density(pedist$pue)$y))
	plot(density(pedist$ple), xlim=xlim, ylim=ylim,...)
	lines(density(pedist$pue))
	if(addq == T){
		pvec	<- seq(0, 1, by=.005)
		plx	<- sapply(pvec, function(pp) length(which(pedist$ple < pp))/nsim)
		pux	<- sapply(pvec, function(pp) length(which(pedist$pue < pp))/nsim)

		Qvec 	<- 1-plx + pux
		Q	<- min(Qvec)
		pcut 	<- pvec[which(Qvec == min(Qvec))]
		legend("topleft", paste("Q=",round(Q, digits=3)))
		abline(v=pcut)
	}
}

