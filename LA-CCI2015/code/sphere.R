sphere <- function(X){
	if(is.matrix(X)) {
		return(sqrt(rowSums(X ^ 2)))
	}
	else {
		return(sqrt(sum(X ^ 2)))
	}
}