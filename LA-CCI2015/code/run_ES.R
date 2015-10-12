# Simple (1+lambda)-ES
# WARNING:  all operations are performed with all variables standardized to the 
# 		interval [0, 1]. Adjust sigma accordingly.
#
# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

run_ES <- function(lambda, sigma, probpars, stopcrit, seed = NULL){
	
	# ===== Error catching
	stopifnot(
		is.numeric(lambda),
		lambda == floor(lambda),
		lambda > 0,
		is.numeric(sigma) && sigma > 0,
		all(c("name", "xmin", "xmax") %in% names(probpars)),
		length(probpars$xmin) == length(probpars$xmax),
		length(probpars$xmin) > 0,
		is.numeric(c(probpars$xmin, probpars$xmax)),
		all(probpars$xmin < probpars$xmax),
		is.null(seed) || (is.numeric(seed) && seed > 0 && seed == floor(seed)))
	
	if(is.null(seed))  seed  <- as.numeric(Sys.time())
	
	
	# ===== Initialize algorithm
	set.seed(seed)				# Set PRNG
	n   <- length(probpars$xmin)		# problem dimension
	t   <- 0      				# iterations counter
	x   <- runif(n)				# initialize random candidate solution
	f   <- evaluate_point(probpars, x)	# evaluate x
	nfe <- 1 					# function evaluations counter
	keep.running <- TRUE 			# iteration control flag
	no.improv    <- 0 			# lack of improvement counter
	
	lhist    <- 1000
	fhist    <- numeric(lhist)			# initialize history vector
	fhist[1] <- f
	
	# ===== Iterative cycle
	while(keep.running){
		# Generate trial vectors
		Xc <- matrix(x,
				 nrow = lambda,
				 ncol = n,
				 byrow = TRUE) + 
			matrix(rnorm(n * lambda, 
					 mean = 0, 
					 sd = sigma),
				 ncol = n,
				 nrow = lambda)
		
		# Truncate to limits
		Xc <- pmax(0*Xc, pmin(0*Xc + 1, Xc))

		# evaluate trial vectors
		fc <- evaluate_population(probpars, Xc)
		nfe <- nfe + lambda
		
		# Selection
		if (any(fc <= f)){
			f <- min(fc)
			x <- Xc[which.min(fc), ]
			no.improv <- 0
		} else{
			no.improv <- no.improv + 1
		}
		
		# Update iteration counter
		t <- t + 1
		
		# Update history vector
		fhist[t + 1] <- f
		
		# Preallocate more space if needed for fhist
		if(t + 1 == lhist){
			fhist <- c(fhist, numeric(lhist))
			lhist <- 2*lhist
		}
		
		# Check stop criteria
		keep.running <- check_stop_criteria()
		
	}
	
	cat("\n\nOptimization finished at iteration", t)
	
	# Return output
	return(list(x     = probpars$xmin + x * (probpars$xmax - probpars$xmin),
			f     = f,
			fhist = fhist[1:(t+1)],
			nfe   = nfe,
			iter  = t))
}