update_package <- function(package_path='/Users/atfrank/GitSoftware/nmR', package_name='nmR'){
	#' Packahe Documentation Update Function
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