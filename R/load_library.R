load_library <- function(package){
 	#' Load Package Function
	#'
	#' This function loads a package if previous installed or install the package and then loads it.
	#' @param package name of package to load.
	if(package %in% rownames(installed.packages()) == FALSE) {install.packages(package,dep=TRUE,repos='http://cran.mtu.edu/')}
	library(package)
}