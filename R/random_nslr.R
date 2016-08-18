random_nslr <- function(n_actives, n_total = 1000, cycles = 10000){
	#' Randomized Sum of Logarithmic Ranks Function
	#'
	#' This function allows you to compute a randomize normalized sum of logarithmic ranks
	#' @param n_actives number of positive samples.
	#' @param n_total total number of samples. Defaults to 1000
	#' @param cycles number of random ranking to generate. Default to 10000.
	#' @export
	#' @examples
	#' random_nslr(10, n_total = 100, cycles = 1000)
	
	data <- data.frame(row=1:n_total)
	data$status <- 0
	data$status[1:n_actives] <- 1

	slrs <- NULL
	for (i in 1:cycles){
		slrs <- c(slrs, nslr(data[sample(1:nrow(data)),]$status))
	}
	mean(slrs)
}
