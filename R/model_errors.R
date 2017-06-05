model_erorrs <- function(cs_input, errorType){		
	#' Model Errors Function
	#'
	#' This function allows you to computes the error for each indiviudal model (or conformer)
	#' @param cs_input chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS, weight.
	#' @param n_total total number of samples. Defaults to 1000
	#' @param cycles number of random ranking to generate. Default to 10000.
	#' @export
	#' @examples
	#' random_nslr(10, n_total = 100, cycles = 1000)

  # Computes the error for each indiviudal model (or conformer)
  # Args:
  #   x: dataframe with columns expCS and predCS -- data frame
  # Returns:
  #		data frame with errors for each model (or conformer) 
  require(plyr)
  errorType <- tolower(errorType)
	if (errorType == "cor_p"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=score_pearson)
	}
	if (errorType == "cor_s"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=score_spearman)
	}
	if (errorType == "cor_k"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=score_kendall)
	}
	if (errorType == "mae"){
			means <- ddply(cs_input, c("model"), .fun=score_mae)
	}
	if (errorType == "rmse"){
			means <- ddply(cs_input, c("model"), .fun=score_rmse)
	}

	if (errorType == "geo_mae"){
			means <- ddply(cs_input, c("model"), .fun=score_geo_mae)
	}
  if (errorType == "prob"){
    means <- ddply(cs_input, c("model"), .fun=score_probability)
  }
  
	if (errorType == "flatchi2_0.25"){
			means <- ddply(cs_input, c("model"), .fun=score_flat_chi2, chi2_c = 0.25)
	}	
	if (errorType == "flatchi2_0.5"){
			means <- ddply(cs_input, c("model"), .fun=score_flat_chi2, chi2_c = 0.5)
	}
	if (errorType == "flatchi2_1"){
			means <- ddply(cs_input, c("model"), .fun=score_flat_chi2, chi2_c = 1)
	}
	if (errorType == "flatchi2_2"){
			means <- ddply(cs_input, c("model"), .fun=score_flat_chi2, chi2_c = 2)
	}
	if (errorType == "flatChi2_3"){
			means <- ddply(cs_input, c("model"), .fun=score_flat_chi2, chi2_c = 3)
	}
  
  colnames(means) <- c("model", "error")
	return(means)
  
}
