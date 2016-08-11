kendall_score <- function(x){
	#' Kendall Correlation Scoring Function
	#'
	#' This function computes the 1 - tau
	#' @param x input dataframe. Should contain field: predCS and expCS.
	#' @export
	#' @examples
	#' kendall_score(x)
	return(1-cor(x$expCS,x$pred,method="kendall"))
}

pearson_score <- function(x){
	#' Pearson Correlation Scoring Function
	#'
	#' This function computes the 1 - R
	#' @param x input dataframe. Should contain field: predCS and expCS.
	#' @export
	#' @examples
	#' pearson_score(x)
	return(1-cor(x$expCS,x$pred,method="pearson"))
}

spearman_score <- function(x){
	#' Spearman Correlation Scoring Function
	#'
	#' This function computes the 1 - rho
	#' @param x input dataframe. Should contain field: predCS and expCS.
	#' @export
	#' @examples
	#' spearman_score(x)
	return(1-cor(x$expCS,x$pred,method="spearman"))
}

score_mae <- function(x){
	#' Weighted Mean Absolute Error Function
	#'
	#' This function computes the weighted (or reduced) mean-absolute-error (MAE)
	#' @param x input dataframe. Should contain field: weight, predCS and expCS.
	#' @export
	#' @examples
	#' score_mae(x)
	x$error <- x$expCS - x$predCS
	return(c(mean(abs(x$weight*x$error))))
}

score_rmse <- function(x){
	#' Weighted Root Mean Squared Error Function
	#'
	#' This function computes the weighted (or reduced) root-mean-squared-error (RMSE)
	#' @param x input dataframe. Should contain field: weight, predCS and expCS.
	#' @export
	#' @examples
	#' score_rmse(x)
  x$error <- x$expCS - x$predCS
	return(sqrt(mean(x$weight*x$weight*x$error*x$error)))
}

score_huber <- function(x){
	#' Pseudo-Huber Error Function
	#'
	#' This function a weighted (or reduced) root-mean-squared-error (RMSE)
	#' @param x input dataframe. Should contain field: w, predCS and expCS.
	#' @export
	#' @examples
	#' score_rmse(x)
	delta <- 1/x$weight
	delta2 <- delta*delta
  a <- x$expCS - x$predCS
  
	return(mean(delta2*(sqrt(1+(a/delta)^2)-1)))
}

score_flat_chi2 <- function(x, chi2_c = 1){
	#' Weighted Flat Chi^2 Error Function
	#'
	#' This function computes  the (or reduced) flat Chi^2 error (return 1 - flat_chi2).
	#' @param x input dataframe. Should contain field: weight, predCS and expCS.
	#' @param chi2_c double. The scaling factor for the error in the exponent.
	#' @export
	#' @examples
	#' score_flat_chi2(x)
	x$error <- x$expCS - x$predCS	
	return(mean(1-exp(-1.0*(x$weight*x$error/chi2_c)**2)))	
}

score_probability <- function(x){
	#' Probability Error Function (similar to the Quality Score)
	#'
	#' This function computes the joint probability of observing the exhibited errors between measured and predicted chemical shifts
	#' @param x input dataframe. Should contain field: weight, predCS and expCS.
	#' @export
	#' @examples
	#' score_flat_chi2(x)
	x$error <- x$expCS - x$predCS
	x$sigma <- 1/x$weight	
	return(1-prod((1/sqrt(2*x$sigma*x$sigma*pi))*exp(-1.0*(x$error/(sqrt(2)*x$sigma)**2))))}
