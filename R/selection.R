# (0) set working directory
setwd("~/GitSoftware/nmR/")

compute_chemical_shifts <- function(w,X){
  # Helper function to compute the chemical shifts
  # Args:
  #   w: weights -- vector 
  #   X: data matrix of features -- data matrix
  # Returns:
  #   returns: the computed chemical shifts shifts as weighted sum
  MX = nrow(X)
  NX = ncol(X)
  .rowSums(matrix(c(w),byrow=T,nrow=MX,ncol=NX)*X,MX,NX)
}

fitness <- function(w, mask = NULL, weights = 1){
  # Helper function to the overall fitness for GA
  # Args:
  #   w: weights -- vector
  # Returns:
  #   returns: correlation coefficient between target and ensemble-averaged chemical shifts
  w <- w/sum(w)
  if (!is.null(mask)){w[mask] <- 0}
  ensemble_averaged <- compute_chemical_shifts(w, ensemble)
  return(-mean(weights*abs(ensemble_averaged-target))-10*(sum(w))) # Note that this is a L1 regularized optimization
}

runGA <- function(target, ensemble, cycles = 100, population_size = 100, seed = 12345, binary = TRUE, weights = 1){
  require(GA)
  # Helper function that runs GA regression to general refined parameters for model
  # Args:
  #   target: actual data -- vector
  #   ensemble: matrix of predicted data -- number of columns correspond to the number of ensemble members and rows the number sample points in the target data
  #   cycles (optional): optimization cycles -- integer
  #   populationSize (optional): population or number of solutions to evolve -- integer
  #   seed: random number seed -- double  
  # Returns:
  #   returns: resulting model parameters -- data frame
  ensemble_size <- ncol(ensemble)
  if (binary){
    GA <- ga(parallel = FALSE, type = "binary", nBits = ensemble_size, fitness = fitness, monitor=TRUE, popSize = population_size, maxiter=cycles, weights = weights)
  } else {
    GA <- ga(parallel = FALSE, type = "real-valued", min = rep(0,25), max = rep(1,25), fitness = fitness, monitor=TRUE, popSize = population_size, maxiter = cycles, weights = weights)
  }
  return(GA)
}

# (1) load library
library(nmR)
target <- as.matrix.data.frame(read.table("data/observed_vector.txt"))
ensemble <- as.matrix.data.frame(read.table("data/predicted_matrix.txt"))
weights <- as.matrix.data.frame(read.table("data/weights_vector.txt"))

GA <- runGA(target, ensemble, binary = FALSE, cycles = 5000, weights = weights)
rmsd <- read.table("data/1SCL.txt")
rmsd$sel <- as.vector(GA@solution)/max(as.vector(GA@solution))
