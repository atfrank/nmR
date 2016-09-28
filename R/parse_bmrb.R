get_bmrb_files <- function(ID,rnaID="TEST",outpath="data/BMRB/",goal="shifts",listID=NULL,offset=0,coffset=FALSE,str=FALSE,nuclei=NULL,residues=NULL,verbose=FALSE){  
  # Main function that uses helper functions to download a BMRB entry and then format the output as
  # **resname resid nucleus chemical_shift chemical_shift_error
  # Args:
  #   ID: BMRB entry ID
  #   rnaID: ID tag is the name of the formatted output file
  #   outpath: path pointing to desired location for the formatted output file
  #   goal: shifts,summary,listcheck
  #   listID: listID to extract from BMRB
  #   rnaID: string appended to 'measured_shifts_' to uniquely identify chemical shifts dataset
  #   offset: renumber residue using by this offset 
  #   coffset: check if 13C chemical shifts have systematic referencing error [default %default]"),
  #   str: use str format for output
  #   nuclei: restrict output to these nuclei
  #   residues: restrict output to these residues
  #   verbose progress information
  #
  # Returns:
  #   no R objects
  #   writes a formatted output file
  
  # source user needed user functions
  source("scripts/identify_cs_offsets.R",local=TRUE)	  
  # get BMRB entry
  y <- download_and_parse_bmrb(ID) 
  print(y$shifts)
  if(!is.null(y$shifts)){
    # check number of chemical shift list in BMRB entry
    if (goal=="listcheck"){
      get_bmrb_list_info(y$shifts)
    }	
    # summarize info in BMRB entry
    if (goal=="summary"){
      get_bmrb_summary(ID, y)
    }	
    # write out chemical shifts
    if (goal=="shifts"){
      write_out_shifts(y,listID=listID,offset=offset,rnaID=rnaID,outpath=outpath,residues=residues,nuclei=nuclei,str=str,verbose=verbose)		
      # check for 13C chemical shifts offsets
      if (coffset){
        referencing_errors_detector(y)
      }
    }
  }
}

download_and_parse_bmrb <- function(ID,verbose=FALSE){
  #' Download and Parse Function
  #'
  #' download the BMRB file specified by ID and generate a chemical shift data frame
  #' @param ID BMRB database ID
  #' @param verbose logical indicating print info to screen
  #' @export
  
  # make temporary directory  
  system("mkdir -p /tmp/bmrb/")		
  
  # download file
  if (verbose){
    download.file(url=paste("http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/",ID,sep=""), destfile="/tmp/bmrb/test.bmrb",quiet=F)
  } else {
    download.file(url=paste("http://rest.bmrb.wisc.edu/bmrb/NMR-STAR3/",ID,sep=""), destfile="/tmp/bmrb/test.bmrb",quiet=T)
  }
  A <- readLines(con="/tmp/bmrb/test.bmrb")
  A <- sub("^[ \t]+","",A) #remove leading whitespace
  A <- sub("[ \t]+$","",A) #remove trailing whitespace
  B <- A[count.fields("/tmp/bmrb/test.bmrb",blank.lines.skip=F)==24|count.fields("/tmp/bmrb/test.bmrb",blank.lines.skip=F)==23]
  B <- B[!grepl("citation",B)]
  B <- B[!grepl("direct",B)]
  B <- B[!grepl("RNA",B)]
  B <- B[!grepl("entity_1",B)]
  B <- B[complete.cases(B)]
  if(length(B)>0){
    y <- read.table(textConnection(B))
    
    if(ncol(y)==24){
      colnames(y) <- c("ID","Assembly_atom_ID","Entity_assembly_ID","Entity_ID","Comp_index_ID","Seq_ID","Comp_ID","Atom_ID","Atom_type","Atom_isotope_number","Val","Val_err","Assign_fig_of_merit","Ambiguity_code","Occupancy","Resonance_ID","Auth_entity_assembly_ID","Auth_asym_ID","Auth_seq_ID","Auth_comp_ID","Auth_atom_ID","Details","Entry_ID","Assigned_chem_shift_list_ID")
    }
    if(ncol(y)==23){
      colnames(y) <- c("ID","Assembly_atom_ID","Entity_assembly_ID","Entity_ID","Comp_index_ID","Seq_ID","Comp_ID","Atom_ID","Atom_type","Atom_isotope_number","Val","Val_err","Assign_fig_of_merit","Ambiguity_code","Occupancy","Resonance_ID","Auth_entity_assembly_ID","Auth_seq_ID","Auth_comp_ID","Auth_atom_ID","Details","Entry_ID","Assigned_chem_shift_list_ID")
    }
    
    fields <- c("Comp_ID","Seq_ID","Atom_ID","Val","Val_err","Assigned_chem_shift_list_ID") 
    y$Seq_ID <- as.numeric(as.character(y$Seq_ID))
    polymers <- A[grepl("_Entity.Polymer_type",A)]
    polymers <- polymers[!grepl("_Entity.Polymer_type_details",polymers)]
    polymers <- read.table(textConnection(polymers),col.names=c("polymer","polymerID"))$polymerID
    
    temperature <- pH <- 9999.9
    conditionsT <- A[grepl("temperature",A)]
    conditionsT <- conditionsT[grepl(ID,conditionsT)]
    conditionsp <- A[grepl("pH",A)]
    conditionsp <- conditionsp[grepl(ID,conditionsp)]
    if (length(conditionsp)==1 && length(conditionsT)==1){
      temperature <- read.table(textConnection(conditionsT),col.names=c("var","T", "unit","error","BMRB","list"))$T
      pH <- read.table(textConnection(conditionsp),col.names=c("var","pH","error","unit","BMRB","list"))$pH
    }
    return(list(shifts=y,poly=polymers,info=A,temp=temperature,pH=pH,fields=fields))
  } else {
    return(list(shifts=NULL))
  }
}

get_bmrb_summary <- function(ID,y){
  #' Summary Information Function
  #'
  #' prints of general info about the BMRB entry (e.g., temperature and pH)
  #' @param ID BMRB database ID
  #' @param y dataframe returned by download_and_parse_bmrb
  #' @export
  polymers <- y$poly
  info <- y$info
  y <- y$shifts
  temperature <- pH <- 9999.9
  conditionsT <- info[grepl("temperature",info)]
  conditionsT <- conditionsT[grepl(ID,conditionsT)]
  conditionsp <- info[grepl("pH",info)]
  conditionsp <- conditionsp[grepl(ID,conditionsp)]
  if (length(conditionsp)==1 && length(conditionsT)==1){
    temperature <- read.table(textConnection(conditionsT),col.names=c("var","T", "unit","error","BMRB","list"))$T
    pH <- read.table(textConnection(conditionsp),col.names=c("var","pH","error","unit","BMRB","list"))$pH
  }
  cat("Assembly Information:\n")
  list_ID <- unique(y$Assigned_chem_shift_list_ID)
  for (i in seq_along(list_ID)){
    r <- range(subset(y,Assigned_chem_shift_list_ID==list_ID[i])[,"Seq_ID"])
    res <- unique(subset(y,Assigned_chem_shift_list_ID==list_ID[i])[,"Comp_ID"])
    cat(sprintf("Assigned_chem_shift_list_ID %s; reside number range %s-%s; %s;\n",i,r[1],r[2],polymers[i]))
  }       
  cat("\nSample Information:\n")
  s <- read.table(textConnection(info[grepl("_Sample_condition_list.ID",info)]),col.names=c("sample","sampleID"))$sampleID
  cat(sprintf("Sample_conditions:\n"))
  cat(sprintf("%s %s\n",temperature,pH))
  cat(sprintf("%s\n",conditionsT))
  cat(sprintf("%s\n",conditionsp))		
  cat("\nAssembly Information:\n")
  cat(sprintf("%s\n",info[grepl("_Assembly.Name",info)]))			
}

get_bmrb_list_info <- function(y){
  #' Count Datasets Function
  #'
  #' counting chemcial shift lists
  #' @param y dataframe returned by download_and_parse_bmrb
  #' @export
  
	list_ID <- unique(y$Assigned_chem_shift_list_ID)
	cat(sprintf("%s\n",length(list_ID)))
	return(list_ID)
}

referencing_errors_detector <- function (input,first=NULL,last=NULL) {
  #' Referencing Error Detection Function
  #'
  #' Identifies whether 13C chemical shift data is systematically shifted. Use the chemical shifts of the a terminal G:C base pair and compare to refernce values
  #' @param input chemical shift data frame (expects columns to correspond to "resname","resid","nucleus","expCS","error" or "resname","resid","nucleus","expCS"
  #' @param first: integer specifying residue number of 5' terminal residue in chemical shift list file
  #' @param last:  integer specifying residue number of 3' terminal residue in chemical shift list file
  #' @export
  
  cs.table <- y
  ifelse(ncol(cs.table) == 5, colnames(cs.table) <- c("resname", "resid", "nucleus", "expCS", "error"), colnames(cs.table) <- c("resname","resid", "nucleus", "expCS"))
  cs.table$type <- "carbon"
  cs.table$type[grep("H",cs.table$nucleus)] <- "proton"
  
  if (nrow(subset(cs.table,type=="carbon"))<1) {
    return("Error: No 13C chemical shifts in dataset")
  }
  
  first.resid <- min(cs.table$resid)
  last.resid <- max(cs.table$resid)
  if (!is.null(first))
    first.resid <- first
  if (!is.null(last))
    last.resid <- last      
  first.resname <- subset(cs.table, resid == first.resid)[1, "resname"]
  last.resname <- subset(cs.table, resid == last.resid)[1, "resname"]
  if (first.resname != "GUA" || last.resname != "CYT") {
    return("Error: No terminal GUA:CYT")
  }
  else {
    comp.first.C8 <- subset(cs.table, resid == first.resid & nucleus == "C8")[, "expCS"]
    comp.last.C1p <- subset(cs.table, resid == last.resid & nucleus == "C1'")[, "expCS"]
    comp.last.C3p <- subset(cs.table, resid == last.resid & nucleus == "C3'")[, "expCS"]
    comp.last.C5 <- subset(cs.table, resid == last.resid & nucleus == "C5")[, "expCS"]
    comp <- c(comp.first.C8, comp.last.C1p, comp.last.C3p, 
              comp.last.C5)
    if (length(comp.first.C8) != 1) {
      comp.first.C8 <- NA
    }
    if (length(comp.last.C1p) != 1) {
      comp.last.C1p <- NA
    }
    if (length(comp.last.C3p) != 1) {
      comp.last.C3p <- NA
    }
    if (length(comp.last.C5) != 1) {
      comp.last.C5 <- NA
    }
    ref.first.C8 <- mean(c(138.7, 139.7))
    ref.last.C1p <- mean(c(92.5, 93.4))
    ref.last.C3p <- mean(c(69.4, 70.4))
    ref.last.C5 <-  mean(c(97.4, 98.8))
    ref <- c(ref.first.C8, ref.last.C1p, ref.last.C3p, ref.last.C5)
    offset <- mean(c(ref.first.C8, ref.last.C1p, ref.last.C3p,ref.last.C5) - c(comp.first.C8, comp.last.C1p, comp.last.C3p, comp.last.C5), na.rm = T)
    if(is.na(offset)){
      return("check this")
    } else {
      return(round(offset,3))
    }        
  }
}

write_out_shifts <- function(y,rnaOnly=TRUE,listID=NULL,offset=offset,rnaID="1234",outpath="",residues=NULL,nuclei=NULL,str=FALSE,verbose=TRUE){
  #' Chemical Shift Output Function
  #'
  #' writes out chemcial shift dataframe to file
  #' @param y dataframe returned by download_and_parse_bmrb
  #' @param rnaOnly logical, write out only RNA chemical shifts?
  #' @param listID writing chemical shifts corresponding to BMRB list specified by listID
  #' @param offset add offset residue numbers
  #' @param rnaID id used in the file name to which chemical shifts are written
  #' @param outpath path pointing to location to save file
  #' @param residues if vector not NULL, will only write out chemical shifts specified residues
  #' @param nuclei if vector not NULL, will only write out chemical shifts specified nuclei
  #' @param str should chemical shifts be written in STAR format?
  #' @param verbose logical containing whether info will be printed to screen
  #' @export
  
  temperature=y$temp
  pH=y$pH
  fields=y$fields
  y <- y$shifts
  if (!is.null(residues)) {
    y <- y[y$Seq_ID %in% residues,]  
  }
  if (!is.null(nuclei)) {
    y <- y[y$Atom_ID %in% nuclei,]  
  }
  if(!is.null(listID)){
		y <- subset(y,Assigned_chem_shift_list_ID==listID)
  } 
  y$Seq_ID <- y$Seq_ID + offset
  y$Comp_ID <- as.character(y$Comp_ID)
  y$Comp_ID[y$Comp_ID=="G"] <- "GUA"
  y$Comp_ID[y$Comp_ID=="A"] <- "ADE"
  y$Comp_ID[y$Comp_ID=="C"] <- "CYT"
  y$Comp_ID[y$Comp_ID=="U"] <- "URA"
  y$Comp_ID[y$Comp_ID=="T"] <- "THY"
  y$temperature <- temperature
  y$pH <- pH
  
  if(rnaOnly)
    y <- y[y$Comp_ID %in% c("GUA","ADE","CYT","URA"),]
  if(nrow(y) > 0){
    if(verbose){
      print(y[,fields])
    }
    if(str){
      y <- data.frame(a=1:nrow(y),b=y$Seq_ID,c=y$Seq_ID,d=substr(as.character(y$Comp_ID),1,1),e=y$Atom_ID,f=substr(as.character(y$Atom_ID),1,1),g=y$Val,h=".",i=".")
      write.table(y,file=paste(outpath,"measured_shifts_",rnaID,".str",sep=""),col.names=F,row.names=F,quote=F)
    } else {
      write.table(y[,c("Comp_ID","Seq_ID","Atom_ID","Val","Val_err")],file=paste(outpath,"measured_shifts_",rnaID,".dat",sep=""),col.names=F,row.names=F,quote=F)
    }
  }
}

