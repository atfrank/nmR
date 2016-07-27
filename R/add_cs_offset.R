add_cs_offset <- function(cs_input, cs_offset){
	#' Chemical Shift Referencing Function
	#'
	#' This function allows you to add an offset to chemical shift based on nucleus type and resname (i.e., residue name)
	#' @param cs_input chemical shift dataframe. Should contain field: resname, resid, nucleus, expCS.
	#' @param cs_offset chemical shift offset dataframe. Should contain field: resname, resid, nucleus, offset.
	#' @export
	#' @examples
	#' add_cs_offset(cs, offsets)
	for (res in unique(cs_offset$resname)){
		for (n in unique(cs_offset$nucleus)){
		  offset <- mean(subset(cs_offset,(nucleus==n&resname==res))$offset)		  
		  if (!is.na(offset) && offset > 1.0){
		  	cat(sprintf("%s %s %s %s\n", id, res, n, offset))
				cs_input$expCS[cs_input$nucleus==n & cs_input$resname==res] <- cs_input$expCS[cs$nucleus==n & cs_input$resname==res] +  offset
			}
		}
	}
	return(cs_input)
}
