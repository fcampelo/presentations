# Implements stop criteria for Evolution Strategies
# This routine accesses the calling environment, which means that changes made 
# in the variables contained in \code{env} WILL change the original values. 
# DO NOT change anything in \code{env} unless you're absolutely sure of what 
# you're doing.
#
# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

check_stop_criteria <- function(){
	
	# Get calling environment
	env 			<- parent.frame()
	crits 		<- env$stopcrit$names
	keep.running 	<- 1
	
	for (crit in crits){
		keep.running <- keep.running * !(do.call(crit,
								     args = list()))
	}
	
	return(as.logical(keep.running))
}

# Stop criterion: maximum number of iterations
stop_maxiter <- function(){
	env <- parent.frame(n = 2)
	return(env$t >= env$stopcrit$maxiter)
}

# Stop criterion: maximum number of objective function calls
stop_maxeval <- function(){
	env <- parent.frame(n = 2)
	return(env$nfe >= env$stopcrit$maxevals)
}

# Stop criterion: statistical test of convergence to epsilon-neighborhood of
# (local) optimum.
# Rnorm is the desired R/sigma ("normalized R") distance to the optimum
stop_stat <- function(){
	env <- parent.frame(n = 2)
	stopifnot(all(c("Rnorm", "sig.level") %in% names(env$stopcrit)),
		    is.numeric(env$stopcrit$Rnorm),
		    env$stopcrit$Rnorm > 0,
		    is.numeric(env$stopcrit$sig.level),
		    env$stopcrit$sig.level > 0 && env$stopcrit$sig.level < 1)
	
	# Probability of successful mutation if the current point is distant from
	# the optimum of a (locally-) spherical function by R (=Rnorm*sigma)
	# Calculate once and store in the main environment of run_ES
	if (!any("mstar" == names(env))) {
		env$mstar <- calc_mstar(alpha  = env$stopcrit$sig.level,
						Rnorm  = env$stopcrit$Rnorm,
						n      = env$n,
						lambda = env$lambda)
	}

	# Echo current state to the console
	if(env$no.improv == 0){
		if(!any("mstar" == names(env))){
			env$last.improv <- 0
		}
		env$stabspan <- env$t - env$last.improv
		env$last.improv <- env$t
		cat("\nIteration ", env$t - 1, ":", 
	    		env$stabspan, "/", env$mstar, " without improvement.")}
	
	# Evaluate and return stop criterion
	return(env$no.improv >= env$mstar)
}


calc_mstar <- function(alpha, Rnorm, n, lambda){
	ps <- calc_ps(n, 1, Rnorm)
	if (lambda == 1){
		return(ceiling(log(alpha / 2) / log(1 - ps)))
	} else {
		za 	<- qnorm(alpha)
		a  	<- (2 * ps / za) ^ 2
		b     <- 4 * ps * (3 * ps - 1)
		c 	<- - (za ^ 2) * (8 * ps + 1)
		d 	<- 4 * ps * (za ^ 2) * (3 * ps + (za ^ 2) * (ps - 1))
		return(ceiling(max(Re(polyroot(c(d, c, b, a)))) / lambda))
	}
}