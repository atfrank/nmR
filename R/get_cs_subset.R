 get_cs_subset <- function(cs_input, nuc="all"){
	#' A Chemical Shift Selection Function
	#'
	#' This function allows you to select a subset of the chemical shift data based on nucleus type
	#' @param cs_input chemical shift dataframe
	#' @param nuc name of nucleus type to select
	#' @export
	#' @examples
	#' get_cs_subset(cs, "H1'")  
	#' get_cs_subset(cs, "all")  
	#' get_cs_subset(cs, "baseCarbon")  
	
	nucleiGroups <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")  
	if((nuc %in% nucleiGroups)){
		cs <- subset(cs_input,nucleus==nuc)
	} else {
		if(nuc=="all"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")
		}
		if(nuc=="noH1primes"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")
		}
		if(nuc=="noH5primes"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H2","H5","H6","H8")
		}
		if(nuc=="proton_baseCarbon"){
			sele <- cs_input$nucleus %in% c("C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")
		}
		if(nuc=="sugar_baseCarbon"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''")
		}
		if(nuc=="baseProton_sugarCarbon"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","H2","H5","H6","H8")
		}
		if(nuc=="sugarProton_baseCarbon"){
			sele <- cs_input$nucleus %in% c("C2","C5","C6","C8","H1'","H2'","H3'","H4'","H5'","H5''")
		}
		if(nuc=="base"){
			sele <- cs_input$nucleus %in% c("C2","C5","C6","C8","H2","H5","H6","H8")
		}
		if(nuc=="sugar"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","H1'","H2'","H3'","H4'","H5'","H5''")
		}
		if(nuc=="carbon"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")
		}
		if(nuc=="proton"){
			sele <- cs_input$nucleus %in% c("H1'","H2'","H3'","H4'","H5'","H5''","H2","H5","H6","H8")
		}
		if(nuc=="baseCarbon"){
			sele <- cs_input$nucleus %in% c("C2","C5","C6","C8")
		}
		if(nuc=="baseProton"){
			sele <- cs_input$nucleus %in% c("H2","H5","H6","H8")
		}
		if(nuc=="sugarCarbon"){
			sele <- cs_input$nucleus %in% c("C1'","C2'","C3'","C4'","C5'")
		}
		if(nuc=="sugarProton"){
			sele <- cs_input$nucleus %in% c("H1'","H2'","H3'","H4'","H5'","H5''")
		}
		cs <- cs_input[sele,]
	}	
	return(cs)
}
