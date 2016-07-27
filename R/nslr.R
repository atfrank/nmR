nslr <- function(X){
	#' Sum of Logarithmic Ranks Function
	#'
	#' This function allows you to compute the normalized sum of logarithmic ranks
	#' @param X vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)
	#' @export
	#' @examples
	#' random_nslr(sample(c(rep(0,100),rep(1,10))))
  
  ri <- which(X==1)
  N <- length(X)
  i <-  1:length(ri)
  SLRmax <- -sum(log(i/N))
  return(-sum(log(ri/N))/SLRmax)
} 
