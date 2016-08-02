load_library <- function(package_name, repo='http://cran.mtu.edu/'){
	#' Package Loader
	#'
	#' This function load a package if installed, and install then loads if not
	#' @param package_name name (string) of package	
  #' @param repo repository (string) to use when installing the package	
	#' @export
	#' @examples
	#' load_library(package_name='plyr')

   if(package_name %in% rownames(installed.packages()) == FALSE) {install.packages(package_name,dep=TRUE,repos=repo)}
   library(package_name)
}

update_package <- function(package_path='/Users/atfrank/GitSoftware/nmR', package_name='nmR'){
	#' Package Documentation Update Function
	#'
	#' This function updates the package specified by package_path and package_name
	#' @param package_path path (string) to package
	#' @param package_name name (string) of package	
	#' @export
	#' @examples
	#' update_package(package_path='/Users/atfrank/GitSoftware/nmR', package_name='nmR')

	library('roxygen2')
	library('devtools')
	setwd(package_path)
	document()
	setwd("..")
	install(package_name)
}