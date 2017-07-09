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
  .rowSums(matrix(c(w),byrow=T,nrow=MX,ncol=NX)*X,MX,NX)/NX
}

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

fitness <- function(w){
  # Helper function to the overall fitness for GA
  # Args:
  #   w: weights -- vector
  # Returns:
  #   returns: correlation coefficient between target and ensemble-averaged chemical shifts
  w <- w/sum(w)
  ensemble_averaged <- compute_chemical_shifts(w, ensemble) 
  #return(cor(target, ensemble_averaged, method = "kendall")-1+sum(w))
  return(-mean(abs(ensemble_averaged-target))-10*(sum(w)))
}

runGA <- function(target, ensemble, cycles = 100, population_size = 25, seed = 12345){
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
  GA <- ga(parallel = FALSE, type = "binary", nBits = ensemble_size, fitness = fitness, monitor=TRUE, seed = seed, popSize = population_size, maxiter = cycles, keepBest = TRUE)
  return(GA)
}


# (1) load library
library(nmR)
target <- as.matrix.data.frame(read.table("data/observed_vector.txt"))
ensemble <- as.matrix.data.frame(read.table("data/predicted_matrix.txt"))
rmsd <- read.table("data/1SCL.txt")

#GA <- ga(parallel = FALSE, type = "binary", nBits = ensemble_size, fitness = fitness, monitor=TRUE, popSize = 100, maxiter=100)
GA <- ga(parallel = FALSE, type = "real-valued", min = rep(0,25), max = rep(1,25), fitness = fitness, monitor=TRUE, popSize = 25, maxiter=20000)


rmsd$sel <- as.vector(GA@solution)/max(as.vector(GA@solution))
