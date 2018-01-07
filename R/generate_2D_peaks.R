create_residue_hmqc_peaks <- function(cs, protons=c("H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), carbons=c("C1'", "C2'", "C3'", "C4'", "C5'", "C5'", "C2", "C5", "C6", "C8")){
  #' HMQC Peaks Generation Helper Function
  #'
  #' This function allows you to convert chemical shift list to chemical shift peak table
  #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
  #' @param protons vector of proton nuclei.
  #' @param carbons vector of carbon nuclei.
  #' @export
  #' @examples
  #' create_residue_hmqc_peaks(cs)

  if(length(protons)!=length(carbons)){stop("list of protons and carbons should be the same")}
  if(length(unique(cs$model))!=1 || length(unique(cs$resid))!=1 ){stop("this function works only for a single model and a single resid. Try using create_peaks")}

  peak_H <- NULL
  peak_C <- NULL
  type_H <- NULL
  type_C <- NULL
  weight_H <- NULL
  weight_C <- NULL

  for (i in 1:length(protons)){
    if(length(cs$predCS[cs$nucleus==protons[i]])!=0){
      peak_H <- c(peak_H, cs$predCS[cs$nucleus==protons[i]])
      peak_C <- c(peak_C, cs$predCS[cs$nucleus==carbons[i]])
      type_H <- c(type_H, protons[i])
      type_C <- c(type_C, carbons[i])
      weight_H <- c(weight_H, cs$weight[cs$nucleus==protons[i]])
      weight_C <- c(weight_C, cs$weight[cs$nucleus==carbons[i]])
    }
  }
  return(data.frame(type_H=type_H, type_C=type_C, peak_H=peak_H, peak_C=peak_C, weight_H=weight_H, weight_C=weight_C))
}

create_residue_tocsy_peaks <- function(cs){
  #' TOCSY Peaks Generation Helper Function
  #'
  #' This function allows you to convert chemical shift list to chemical shift peak table
  #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
  #' @export
  #' @examples
  #' create_residue_tocsy_peaks(cs)

  # generate a simulated 2D TOCSY spectrum from a list of assigned peaks
  nuc <- c("H1'","H2'","H3'","H4'","H5'","H2","H5","H6","H8")
  nuc_name <- c("h1p","h2p","h3p","h4p","h5p","h2","h5","h6","h8")
  for (i in seq_along(nuc_name)){
    if(nrow(subset(cs, nucleus==nuc[i]))<1){
      assign(nuc_name[i], 9999)
    }
    else {
      assign(nuc_name[i],subset(cs,nucleus==nuc[i])$predCS[1])
    }
  }
  p1 <- p2 <- p1_nam <- p2_nam <- NULL
  # "H1'" correlations
  p1 <- c(p1,h1p,h1p,h1p,h1p,h2p,h3p,h4p,h5p)
  p2 <- c(p2,h2p,h3p,h4p,h5p,h1p,h1p,h1p,h1p)
  p1_nam <- c(p1_nam,"H1'","H1'","H1'","H1'","H2'","H3'","H4'","H5'")
  p2_nam <- c(p2_nam,"H2'","H3'","H4'","H5'","H1'","H1'","H1'","H1'")

  # "H2'" correlations
  p1 <- c(p1,h2p,h2p,h2p,h2p,h1p,h3p,h4p,h5p)
  p2 <- c(p2,h1p,h3p,h4p,h5p,h2p,h2p,h2p,h2p)
  p1_nam <- c(p1_nam,"H2'","H2'","H2'","H2'","H1'","H3'","H4'","H5'")
  p2_nam <- c(p2_nam,"H1'","H3'","H4'","H5'","H2'","H2'","H2'","H2'")

  # "H3'" correlations
  p1 <- c(p1,h3p,h3p,h3p,h3p,h1p,h2p,h4p,h5p)
  p2 <- c(p2,h1p,h2p,h4p,h5p,h3p,h3p,h3p,h3p)
  p1_nam <- c(p1_nam,"H3'","H3'","H3'","H3'","H1'","H2'","H4'","H5'")
  p2_nam <- c(p2_nam,"H1'","H2'","H4'","H5'","H3'","H3'","H3'","H3'")

  # "H4'" correlations
  p1 <- c(p1,h4p,h4p,h4p,h4p,h1p,h2p,h3p,h5p)
  p2 <- c(p2,h1p,h2p,h3p,h5p,h4p,h4p,h4p,h4p)
  p1_nam <- c(p1_nam,"H4'","H4'","H4'","H4'","H1'","H2'","H3'","H5'")
  p2_nam <- c(p2_nam,"H1'","H2'","H3'","H5'","H4'","H4'","H4'","H4'")

  # "H5'" correlations
  p1 <- c(p1,h5p,h5p,h5p,h5p,h1p,h2p,h3p,h4p)
  p2 <- c(p2,h1p,h2p,h3p,h4p,h5p,h5p,h5p,h5p)
  p1_nam <- c(p1_nam,"H5'","H5'","H5'","H5'","H1'","H2'","H3'","H4'")
  p2_nam <- c(p2_nam,"H1'","H2'","H3'","H4'","H5'","H5'","H5'","H5'")

  resname <- unique(x$resname)
  if (resname =="URA" || resname == "CYT"){
    p1 <- c(p1,h5,h6)
    p2 <- c(p2,h6,h5)

    p1_nam <- c(p1_nam,"H5","H6")
    p2_nam <- c(p2_nam,"H6","H5")
  }

  if ( resname == "ADE" ){
    p1 <- c(p1,h2,h6,h8,h6)
    p2 <- c(p2,h6,h2,h6,h8)

    p1_nam <- c(p1_nam,"H2","H6","H8","H6")
    p2_nam <- c(p2_nam,"H6","H2","H6","H8")
  }

  spectrum <- data.frame(pair=paste(p1_nam,p2_nam,sep=":"),cs1=p1,cs2=p2)
  spectrum_H <- subset(spectrum,cs1!=9999 & cs2!=9999)

  # generate a simulated 2D TOCSY spectrum from a list of assigned peaks
  nuc <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")
  nuc_name <- c("c1p","c2p","c3p","c4p","c5p","c2","c5","c6","c8")
  for (i in seq_along(nuc_name)){
    if(nrow(subset(x,nucleus==nuc[i]))<1){
      assign(nuc_name[i], 9999)
    }
    else {
      assign(nuc_name[i],subset(x,nucleus==nuc[i])$predCS[1])
    }
  }
  p1 <- p2 <- p1_nam <- p2_nam <- NULL
  # "C1'" correlations
  p1 <- c(p1,c1p,c1p,c1p,c1p,c2p,c3p,c4p,c5p)
  p2 <- c(p2,c2p,c3p,c4p,c5p,c1p,c1p,c1p,c1p)
  p1_nam <- c(p1_nam,"C1'","C1'","C1'","C1'","C2'","C3'","C4'","C5'")
  p2_nam <- c(p2_nam,"C2'","C3'","C4'","C5'","C1'","C1'","C1'","C1'")

  # "C2'" correlations
  p1 <- c(p1,c2p,c2p,c2p,c2p,c1p,c3p,c4p,c5p)
  p2 <- c(p2,c1p,c3p,c4p,c5p,c2p,c2p,c2p,c2p)
  p1_nam <- c(p1_nam,"C2'","C2'","C2'","C2'","C1'","C3'","C4'","C5'")
  p2_nam <- c(p2_nam,"C1'","C3'","C4'","C5'","C2'","C2'","C2'","C2'")

  # "C3'" correlations
  p1 <- c(p1,c3p,c3p,c3p,c3p,c1p,c2p,c4p,c5p)
  p2 <- c(p2,c1p,c2p,c4p,c5p,c3p,c3p,c3p,c3p)
  p1_nam <- c(p1_nam,"C3'","C3'","C3'","C3'","C1'","C2'","C4'","C5'")
  p2_nam <- c(p2_nam,"C1'","C2'","C4'","C5'","C3'","C3'","C3'","C3'")

  # "C4'" correlations
  p1 <- c(p1,c4p,c4p,c4p,c4p,c1p,c2p,c3p,c5p)
  p2 <- c(p2,c1p,c2p,c3p,c5p,c4p,c4p,c4p,c4p)
  p1_nam <- c(p1_nam,"C4'","C4'","C4'","C4'","C1'","C2'","C3'","C5'")
  p2_nam <- c(p2_nam,"C1'","C2'","C3'","C5'","C4'","C4'","C4'","C4'")

  # "C5'" correlations
  p1 <- c(p1,c5p,c5p,c5p,c5p,c1p,c2p,c3p,c4p)
  p2 <- c(p2,c1p,c2p,c3p,c4p,c5p,c5p,c5p,c5p)
  p1_nam <- c(p1_nam,"C5'","C5'","C5'","C5'","C1'","C2'","C3'","C4'")
  p2_nam <- c(p2_nam,"C1'","C2'","C3'","C4'","C5'","C5'","C5'","C5'")

  resname <- unique(x$resname)
  if (resname =="URA" || resname == "CYT"){
    p1 <- c(p1,c5,c6)
    p2 <- c(p2,c6,c5)

    p1_nam <- c(p1_nam,"C5","C6")
    p2_nam <- c(p2_nam,"C6","C5")
  }

  if ( resname == "ADE" ){
    p1 <- c(p1,c2,c6,c8,c6)
    p2 <- c(p2,c6,c2,c6,c8)

    p1_nam <- c(p1_nam,"C2","C6","C8","C6")
    p2_nam <- c(p2_nam,"C6","C2","C6","C8")
  }

  spectrum <- data.frame(pair=paste(p1_nam,p2_nam,sep=":"),cs1=p1,cs2=p2)
  spectrum_C <- subset(spectrum,cs1!=9999 & cs2!=9999)
  return(rbind(spectrum_C, spectrum_H))
}

create_residue_cosy_peaks<- function(cs){
  #' COSY Peaks Generation Helper Function
  #'
  #' This function allows you to convert chemical shift list to chemical shift peak table
  #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
  #' @export
  #' @examples
  #' create_residue_cosy_peaks(cs)

  # generate a simulated 2D COSY spectrum from a list of assigned peaks
  nuc <- c("H1'","H2'","H3'","H4'","H5'","H2","H5","H6","H8")
  nuc_name <- c("h1p","h2p","h3p","h4p","h5p","h2","h5","h6","h8")
  for (i in seq_along(nuc_name)){
    if(nrow(subset(cs,nucleus==nuc[i]))<1){
      assign(nuc_name[i], 9999)
    }
    else {
      assign(nuc_name[i],subset(cs,nucleus==nuc[i])$predCS[1])
    }
  }
  p1 <- p2 <- NULL
  p1 <- c(h1p,h2p,h2p,h3p,h3p,h4p,h4p,h5p)
  p2 <- c(h2p,h1p,h3p,h2p,h4p,h3p,h5p,h4p)

  p1_nam <- c("H1'","H2'","H2'","H3'","H3'","H4'","H4'","H5'")
  p2_nam <- c("H2'","H1'","H3'","H2'","H4'","H3'","H5'","H4'")

  resname <- unique(x$resname)
  if (resname =="URA" || resname == "CYT"){
    p1 <- c(p1,h5,h6)
    p2 <- c(p2,h6,h5)

    p1_nam <- c(p1_nam,"H5","H6")
    p2_nam <- c(p2_nam,"H6","H5")
  }

  spectrum <- data.frame(pair=paste(p1_nam,p2_nam,sep=":"),cs1=p1,cs2=p2)
  spectrum_H <- subset(spectrum,cs1!=9999 & cs2!=9999)

  # generate a simulated 2D COSY spectrum from a list of assigned peaks
  nuc <- c("C1'","C2'","C3'","C4'","C5'","C2","C5","C6","C8")
  nuc_name <- c("c1p","c2p","c3p","c4p","c5p","c2","c5","c6","c8")
  for (i in seq_along(nuc_name)){
    if(nrow(subset(x,nucleus==nuc[i]))<1){
      assign(nuc_name[i], 9999)
    }
    else {
      assign(nuc_name[i],subset(x,nucleus==nuc[i])$predCS[1])
    }
  }
  p1 <- p2 <- NULL
  p1 <- c(c1p,c2p,c2p,c3p,c3p,c4p,c4p,c5p)
  p2 <- c(c2p,c1p,c3p,c2p,c4p,c3p,c5p,c4p)

  p1_nam <- c("C1'","C2'","C2'","C3'","C3'","C4'","C4'","C5'")
  p2_nam <- c("C2'","C1'","C3'","C2'","C4'","C3'","C5'","C4'")

  resname <- unique(x$resname)
  if (resname =="URA" || resname == "CYT"){
    p1 <- c(p1,c5,c6)
    p2 <- c(p2,c6,c5)

    p1_nam <- c(p1_nam,"C5","C6")
    p2_nam <- c(p2_nam,"C6","C5")
  }

  spectrum <- data.frame(pair=paste(p1_nam,p2_nam,sep=":"),cs1=p1,cs2=p2)
  spectrum_C <- subset(spectrum,cs1!=9999 & cs2!=9999)
  return(rbind(spectrum_C, spectrum_H))
}

create_peaks <- function(cs, type = "hmqc", grouping = c("model", "resid", "resname"), protons=c("H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), carbons=c("C1'", "C2'", "C3'", "C4'", "C5'", "C5'", "C2", "C5", "C6", "C8")){
  #'  Peaks Generation Function
  #'
  #' This function allows you to convert chemical shift list to chemical shift peak table
  #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
  #' @param type type of 2D spectra.
  #' @param grouping variables used to group data.
  #' @param protons vector of proton nuclei.
  #' @param carbons vector of carbon nuclei.

  #' @export
  #' @examples
  #' create_hmqc_peaks(cs)
  require(plyr)
  if(type == "hmqc"){
    peaks <- plyr::ddply(.data = cs, .variables = grouping, .fun = create_residue_hmqc_peaks, protons, carbons)
  }

  if(type == "cosy"){
    peaks <- plyr::ddply(.data = cs, .variables = grouping, .fun = create_residue_cosy_peaks)
  }

  if(type == "tocsy"){
    peaks <- plyr::ddply(.data = cs, .variables = grouping, .fun = create_residue_tocsy_peaks)
  }
  return(peaks)
}

compass_score <- function (Q, P) {
  #' COMPASS Scoring Function
  #'
  #' This function allows you to compare two multi-dimensional spectra
  #' See: 10.1016/j.str.2015.07.019: Experimental Protein Structure Verification by Scoring with a Single, Unassigned NMR Spectrum
  #' @param Q reference spectrum (experimental)
  #' @param P comparison spectrum (simulated)
  #' @export
  #' @examples
  #' compass_score(cs)

  # computes the COMPASS score between points in 2D space
  stopifnot(is.numeric(P), is.numeric(Q))
  if (is.vector(P))
    P <- matrix(P, ncol = 1)
  if (is.vector(Q))
    Q <- matrix(Q, ncol = 1)
  if (ncol(P) != ncol(Q))
    stop("'P' and 'Q' must have the same number of columns.")
  D <- pracma::distmat(Q, P)
  return(list(scores=median(apply(D, 1, min)), indices=apply(D, 1, which.min)))
}

