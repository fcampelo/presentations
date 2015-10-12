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