detect_referencing_error <- function(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based = FALSE, verbose = FALSE){
  #' Automated Referencing Error Detection Using Bayesian Regression
  #'
  #' This function allows user to detect potential referencing errors by comparing observed and computed chemical shifts
  #' @param observed_chemical_shifts observed chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS.
  #' @param computed_chemical_shifts observed chemical shift dataframe. Should contain field: resname, resid, nucleus, predCS
  #' @param residue_and_nucleus_based if TRUE, both residue type and nucleus are considered when detecting errors. Default considers only nucleus type
  #' @param verbose if TRUE, print progress log from MCMCpack. Default is FALSE
  #' @export
  #' @examples
  #' detect_referencing_error(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based=TRUE)
  #' detect_referencing_error(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based=FALSE)
  require(MCMCpack)
  cs <- merge(observed_chemical_shifts, computed_chemical_shifts)
  cs <- subset(cs, expCS!=0)
  cs$secondary <- cs$expCS - cs$predCS
  if (!residue_and_nucleus_based){
    f <- MCMCregress(secondary~nucleus+0, data = cs, burnin = 1000, mcmc = 10000, thin = 1, verbose = verbose)
  } else {
    f <- MCMCregress(secondary~nucleus+resname+0, data = cs, burnin = 1000, mcmc = 10000, thin = 1, verbose = verbose)
  }
  return(summary(f))
}

accuracy_estimate <- function(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based = FALSE, verbose = FALSE){
  #' Automated Referencing Error Detection Using Bayesian Regression
  #'
  #' This function allows user to detect potential referencing errors by comparing observed and computed chemical shifts
  #' @param observed_chemical_shifts observed chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS.
  #' @param computed_chemical_shifts computed chemical shift dataframe. Should contain field: resname, resid, nucleus, predCS
  #' @param residue_and_nucleus_based if TRUE, both residue type and nucleus are considered when detecting errors. Default considers only nucleus type
  #' @param verbose if TRUE, print progress log from MCMCpack. Default is FALSE
  #' @export
  #' @examples
  #' detect_referencing_error(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based=TRUE)
  #' detect_referencing_error(observed_chemical_shifts, computed_chemical_shifts, residue_and_nucleus_based=FALSE)
  require(MCMCpack)
  cs <- merge(observed_chemical_shifts, computed_chemical_shifts)
  cs <- subset(cs, expCS!=0)
  cs$secondary <- abs(cs$expCS - cs$predCS)
  if (!residue_and_nucleus_based){
    f <- MCMCregress(secondary~nucleus, data = cs, burnin = 1000, mcmc = 10000, thin = 1, verbose = verbose)
  } else {
    f <- MCMCregress(secondary~nucleus+resname, data = cs, burnin = 1000, mcmc = 10000, thin = 1, verbose = verbose)
  }
  return(summary(f))
}

correct_referencing_error <- function(observed_chemical_shifts, computed_chemical_shifts, ratio_mean_sd = 3, threshold = -2, verbose = FALSE){
  #' Correct Referencing Errors
  #'
  #' This function allows user to detect potential referencing errors by comparing observed and computed chemical shifts
  #' @param observed_chemical_shifts observed chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS
  #' @param computed_chemical_shifts observed chemical shift dataframe. Should contain field: resname, resid, nucleus, predCS
  #' @param ratio_mean_sd ratio of mean to sd of the error. 
  #' @param threshold cutoff to use for identifying systematic errors.
  #' @param verbose if TRUE, print progress log from MCMCpack. Default is FALSE
  #' @export
  #' @examples
  #' correct_referencing_error(observed_chemical_shifts, computed_chemical_shifts)

  err <- detect_referencing_error(observed_chemical_shifts, computed_chemical_shifts, FALSE, verbose)
  err <- as.data.frame(err$statistics)
  err$ratio <- abs(err$Mean/err$SD)
  err <- err[((err$ratio > ratio_mean_sd) & (err$Mean < threshold)),]
  err <- err[!(rownames(err) %in% "sigma2"),]
  err$nucleus <- substr(rownames(err), 8, 12)
  err$Mean <- round(err$Mean, 3)
  err <- err[,(colnames(err) %in% c("nucleus", "Mean"))]
  observed_chemical_shifts <- merge(observed_chemical_shifts, err, all = TRUE)
  observed_chemical_shifts$Mean[is.na(observed_chemical_shifts$Mean)] <- 0.0
  observed_chemical_shifts$Mean <- abs(observed_chemical_shifts$Mean)
  observed_chemical_shifts$cs_corrected <- observed_chemical_shifts$expCS + abs(observed_chemical_shifts$Mean)
  expnames <- c("resname","resid","nucleus","cs_corrected","error")
  observed_chemical_shifts <- observed_chemical_shifts[,expnames]
  return(list(shifts=observed_chemical_shifts, errors=err))
}
