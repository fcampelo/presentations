# Calculate the probability of a successful mutation in a (locally-) spherical 
# function, in terms of :
# - n 	: dimension of the problem
# - sigma	: standard deviation of the isotropic Gaussian mutation 
# - R		: distance from the optimum
#
# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

calc_ps <- function(n, sigma, R){
	
	# ===== Error catching
	stopifnot(
		is.numeric(c(n, sigma, R)), # all inputs must be numeric
		n == floor(n),              # n must be an integer
		n > 1,                      # n must be greater than 1
		sigma > 0,                  # sigma must be positive
		R > 0)                      # R must be positive
	
	# ===== Calculate Phi(n)
	if (n == 2){
		Phi <- 1 / pi
	} else {
		jvec <- 0:(n-3)
		intvals <- unlist(lapply(jvec,
						 FUN  = int_psine,
						 linf = 0,
						 lsup = pi))
		Phi <- (sqrt(pi) ^ (-n)) * prod(intvals)
	}
	
	# ===== Calculate probability of success, p_s(n, sigma, R)
	
	# Integration argument (function)
	funarg <- function(phi1, n, sigma, R){
		# Lower incomplete gamma function
		ligf <- pgamma(2 * (R * cos(phi1) / sigma) ^ 2, n / 2) * gamma(n / 2)
		
		# multiply by power sine
		ligf * sin(phi1) ^ (n - 2)
	}
	
	# Phi(n) * int_{0}^{pi/2} [funarg(phi1, n, sigma, R)] d(phi1)
	Phi * integrate(f     = funarg, 
			    lower = 0, 
			    upper = pi / 2,
			    n     = n,
			    sigma = sigma,
			    R     = R)$value
}


# integral of a power sine function: 
# \int\limits_{linf}^{lsup} [(sin(t)) ^ n] dt
int_psine <- function(n, linf, lsup){
	integrate(f     = function(x, n) (sin(x)) ^ n,
		    lower = linf,
		    upper = lsup,
		    n     = n)$value
}