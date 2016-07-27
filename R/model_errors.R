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
	if (errorType == "COR_p"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=getCOR_p)
	}
	if (errorType == "COR_s"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=getCOR_s)
	}
	if (errorType == "COR_k"){
			means <- ddply(.data=cs_input, .var=c("model"),.fun=getCOR_k)
	}
	if (errorType == "RMSE"){
			means <- ddply(cs_input, c("model"), .fun=getRMSE)
	}
	if (errorType == "MAE"){
			means <- ddply(cs_input, c("model"), .fun=getMAE)
	}
	if (errorType == "rRMSE"){
			means <- ddply(cs_input, c("model"), .fun=getrRMSE)
	}
	if (errorType == "rMAE"){
			means <- ddply(cs_input, c("model"), .fun=getrMAE)
	}
	if (errorType == "LOGIC"){
			means <- ddply(cs_input, c("model"), .fun=getLOGIC)
	}
	if (errorType == "rMAE_COR_k"){
			means <- ddply(cs_input, c("model"), .fun=getrMAE_COR_k)
	}
	if (errorType == "rRMSE_COR_k"){
			means <- ddply(cs_input, c("model"), .fun=getrRMSE_COR_k)
	}
	if (errorType == "flatChi2_3_COR_k"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2_3_COR_k)
	}	

	if (errorType == "rMAE_COR_p"){
			means <- ddply(cs_input, c("model"), .fun=getrMAE_COR_p)
	}
	if (errorType == "rRMSE_COR_p"){
			means <- ddply(cs_input, c("model"), .fun=getrRMSE_COR_p)
	}
	if (errorType == "flatChi2_3_COR_p"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2_3_COR_p)
	}	
	if (errorType == "rMAE_COR_s"){
			means <- ddply(cs_input, c("model"), .fun=getrMAE_COR_s)
	}
	if (errorType == "rRMSE_COR_s"){
			means <- ddply(cs_input, c("model"), .fun=getrRMSE_COR_s)
	}
	if (errorType == "flatChi2_3_COR_s"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2_3_COR_s)
	}	

	if (errorType == "flatChi2_0.125"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 0.125)
	}

	if (errorType == "flatChi2_0.25"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 0.25)
	}	
	if (errorType == "flatChi2_0.5"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 0.5)
	}
	if (errorType == "flatChi2_1"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 1)
	}
	if (errorType == "flatChi2_2"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 2)
	}
	if (errorType == "flatChi2_3"){
			means <- ddply(cs_input, c("model"), .fun=getflatChi2, chi2_c = 3)
	}
	return(means)
}
