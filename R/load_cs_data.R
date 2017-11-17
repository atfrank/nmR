 load_cs_data <- function(csfile, remove_assigned=FALSE, accuracyFile=NULL, cs_offset=NULL, residues=NULL, atomBasedWeights=FALSE, names=c("processor","model","resid","resname","nucleus","predCS","expCS","randCS","id")){
	#' Chemical Shift Data Loading Function
	#'
	#' This function allows you to load chemical shift data from file
	#' @param csfile input file name
	#' @param remove_assigned remove data with zero values? Defaults to FALSE
	#' @param accuracyFile expected accuracy file. Defaults to NULL
	#' @param cs_offset file containing offsets. Defaults to NULL
	#' @param residues vector of residues to retain data for. Defaults to NULL
	#' @param atomBasedWeights vector of specifying the residues for to data will be retained. Defaults to NULL
	#' @param names vector of column names. Defaults to c("processor","model","resid","resname","nucleus","predCS","expCS","randCS","id")
	#' @export
	#' @examples
	#' load_cs_data("data/predicted_shifts.txt")

	# read in raw chemical shift data
	cs <- read.table(csfile,col.names=names)

	# remove unassigned chemical shifts
	if(remove_assigned){cs <- subset(cs,expCS!=0)}

	# retain chemical shifts for specified residues
	if (!is.null(residues)) {
			cs <- cs[cs$resid %in% residues,]  
	}

	# retain only protonated chemical shifts
	cs <- subset(cs,!((resname=="ADE")&(nucleus=="C5")))
	cs <- subset(cs,!((resname=="ADE")&(nucleus=="C6")))
	cs <- subset(cs,!((resname=="GUA")&(nucleus=="C5")))
	cs <- subset(cs,!((resname=="GUA")&(nucleus=="C6")))
	cs <- subset(cs,!((resname=="URA")&(nucleus=="C2")))
	cs <- subset(cs,!((resname=="CYT")&(nucleus=="C2")))	
	
	
  # add differential error weights to carbon and proton nuclei
  if(!is.null(accuracyFile)){
		# read in accuracy stats
		accuracyFile <- read.table(accuracyFile, header = TRUE)  
  	cs <- add_cs_weights(cs, accuracyFile, atomBasedWeights)
  }  
	
  # add offets
  if(!is.null(cs_offset)){
		# read in accuracy stats
		cs_offset <- read.table(cs_offset, col.names = c("resname", "nucleus", "offset"))  
		cs <- add_cs_offset(cs, cs_offset) 
  }
	return(cs)
}
