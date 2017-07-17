get_ensemble_average <- function(w, X){
  #' Ensemble Averaging Function
  #'
  #' This function compute ensemble-averaged of X using w
  #' @author Aaron T. Frank
  #' @param w weights vector 
  #' @param X data matrix to be averaged
  #' @export
  #' @examples
  #' get_ensemble_average(w, X)  
  MX = nrow(X)
  NX = ncol(X)
  .rowSums(matrix(c(w),byrow=T,nrow=MX,ncol=NX)*X,MX,NX)
}

fitness <- function(p, mask = NULL, weights = 1, alpha = 10){
  #' GA Fitness Function
  #'
  #' This function compute ensemble-averaged of X using w
  #' @author Aaron T. Frank
  #' @param p GA parameters (vector)
  #' @param mask mask certain parameters (vector)
  #' @param alpha coefficient to L1 regularization term
  #' @export
  #' @examples
  #' fitness(p)  
  p <- p/sum(p)
  if (!is.null(mask)){p[mask] <- 0}
  ensemble_averaged <- get_ensemble_average(p, ensemble)
  return(-mean(weights*abs(ensemble_averaged-target))-alpha*(sum(p))) # Note that this is a L1 regularized optimization
}

run_ga_selection <- function(target, ensemble, cycles = 100, population_size = 100, seed = 12345, binary = TRUE, weights = 1, monitor = FALSE){
  require(GA)
  #' GA Optimization Function
  #'
  #' This function runs GA optimizations
  #' @author Aaron T. Frank
  #' @param target Actual data (vector)
  #' @param ensemble Matrix of predicted data. Number of columns correspond to the number of ensemble members and rows the number sample points in the target data.
  #' @param alpha Coefficient to L1 regularization term
  #' @param cycles GA optimization cycles (integer)
  #' @param population_size GA population size (integer)
  #' @param seed Seed for random number generator
  #' @param binary Run binary selection mode?
  #' @param weights Weights used to individual data points
  #' @param monitor Print progress?
  #' @export
  #' @examples
  #' fitness(p)  
  ensemble_size <- ncol(ensemble)
  if (binary){
    GA <- ga(parallel = FALSE, type = "binary", nBits = ensemble_size, fitness = fitness, monitor = monitor, popSize = population_size, maxiter = cycles, weights = weights)
  } else {
    GA <- ga(parallel = FALSE, type = "real-valued", min = rep(0, ensemble_size), max = rep(1, ensemble_size), fitness = fitness, monitor = monitor, popSize = population_size, maxiter = cycles, weights = weights)
  }
  return(GA)
}
