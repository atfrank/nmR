enrichment_factor <- function(X, fraction=0.05){
	#' Enrichment Factor Function
	#'
	#' This function computes the enrichment factor
	#' @param X vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)
	#' @export
	#' @examples
	#' enrichment_factor(sample(c(rep(0,100),rep(1,10))))
  
  ri <- which(X==1)
  N <- length(X)
  n_actives <-  length(ri)    
  B <- n_actives/N
  
  X <- X[1:ceiling(fraction*N)]
  ri <- which(X==1)
  N <- length(X)
  n_actives <-  length(ri)
  A  <- n_actives/N  
  return(A/B)
} 

relative_enrichment_factor <- function(X, fraction=0.10){
	#' Relative Enrichment Factor Function
	#'
	#' This function computes the relative enrichment factor
	#' @param X vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)
	#' @export
	#' @examples
	#' relative_enrichment_factor(sample(c(rep(0,100),rep(1,10))))
  N <- length(X)
  n <-  length(which(X==1))      
  X <- X[1:ceiling(fraction*N)]
  na <-  length(which(X==1))      
  N <- length(X)
  return(na/min(c(N, n)))
} 
