# If you use this code, please cite:
# F. Campelo, "Towards Statistical Convergence Criteria for Mutation-based
# 	Evolutionary Algorithms", 2nd Latin American Congress on Evolutionary 
#	Computation - LA-CCI, Curitiba, Brazil, 2015.

library(ggplot2)

# Define server logic
shinyServer(function(input, output) {
	# Calculate m*
	output$mstar <- renderText(
		paste("Required iterations without improvement:",
			calc_mstar(alpha  = input$alpha, 
				     Rnorm  = input$Rnorm,
				     n      = input$N,
				     lambda = input$lambda)))
	
	# Run optimization after user clicks RUN 
	observeEvent(input$run,{
		# Optimization run - setup
		probpars <- list(name  = "sphere", 
				     xmin = rep(-0.5, input$N), 
				     xmax = rep(0.5,  input$N))
		stopcrit <- list(names = c("stop_stat", "stop_maxiter"), 
				     Rnorm = input$Rnorm, 
				     sig.level = input$alpha,
				     maxiter = ifelse(input$maxiter==0,
				     		     Inf,
				     		     input$maxiter))
		
		# Run optimization
		out <- run_ES(input$lambda, input$sigma, probpars, stopcrit)
		
		output$viewplot <- renderPlot({
			ggplot(data = data.frame(t = seq_along(out$fhist),
							 f = out$fhist),
				 aes(x = t, y = f)) + 
				geom_line() + geom_point() +
				ggtitle("Algorithm progress") + 
				xlab("Iteration") + ylab ("Function value") +
				theme(axis.title = element_text(size = 16),
					axis.text  = element_text(size = 14),
					plot.title = element_text(size = 18)) + 
				annotate("text", 
					   x = 0.1*out$iter, y = out$fhist[1], 
					   label = paste("Final value:", signif(out$f,6), 
					   		  "\nTotal iters:", out$iter), 
					   hjust = 0) + 
				scale_x_log10()
		})
	})
	
	output$counter <- 
	  renderText({
	    if (!file.exists("counter.Rdata")) counter <- 0
	    else load("counter.Rdata")
	    counter  <- counter + 1
	    save(counter, 
	         file = "counter.Rdata")     
	    paste("Hits: ", counter)
	  })
})




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




evaluate_point <- function (probpars, x){
  
  # Denormalize point
  x <- probpars$xmin + x * (probpars$xmax - probpars$xmin)
  
  # Evaluate candidate solution
  y <- do.call(probpars$name,
               args = list(x))
  
  # Return evaluated value
  return(y)
}

# If we have a population
evaluate_population <- function (probpars, X){
  
  # Denormalize population
  X <- denormalize_population(probpars, X)
  
  # Evaluate candidate solutions
  Z <- do.call(probpars$name,
               args = list(X))
  
  # Return evaluated values
  return (Z)
}

# Denormalize population
denormalize_population <- function(probpars, X){
  # Denormalize population
  LL <- matrix(rep(probpars$xmin, nrow(X)),
               ncol = ncol(X),
               byrow = TRUE)
  UL <- matrix(rep(probpars$xmax, nrow(X)),
               ncol = ncol(X),
               byrow = TRUE)
  return(LL + X * (UL - LL))
}




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




sphere <- function(X){
  if(is.matrix(X)) {
    return(sqrt(rowSums(X ^ 2)))
  }
  else {
    return(sqrt(sum(X ^ 2)))
  }
}