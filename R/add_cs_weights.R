add_cs_weights <- function(cs_input, weights, atomBasedWeights=FALSE){	
	#' Differential Weighting Function
	#'
	#' This function allows you to add weights chemical shift data 
	#' @param cs_input chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS.
	#' @param weights chemical shift weights dataframe. Should contain field: nucleus and weight.
	#' @param atomBasedWeights should weights by added based on nucleus type and resname. If TRUE, weights dataframe should also contain a resname field. Defaults to FALSE.
	#' @export
	#' @examples
	#' add_cs_weights(cs, weights, atomBasedWeights=FALSE)
	#' add_cs_weights(cs, weights, atomBasedWeights=TRUE)
	
	cs_input$type <- "carbon"
	cs_input$type[grep("H",cs_input$nucleus)] <- "proton"
	cs_input$type[grep("N",cs_input$nucleus)] <- "nitrogen"
	cs_input$weight <- 1.0
	
	if(atomBasedWeights){
		for (res in c("GUA","ADE","CYT","URA")){
			for (n in unique(weights$nucleus)){
				cs_input$weight[cs_input$nucleus==n & cs_input$resname==res] <- mean(subset(weights,(nucleus==n&resname==res))$weight)
			}
		}
	} else {
		for (n in weights$nucleus){
			cs_input$weight[cs_input$nucleus==n] <- mean(subset(weights,(nucleus==n))$weight)
		}
	}	
	cs_input$weight <- 1/(cs_input$weight)
	return(cs_input)
}
