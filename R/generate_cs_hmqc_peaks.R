 create_residue_peaks <- function(cs, protons=c("H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), carbons=c("C1'", "C2'", "C3'", "C4'", "C5'", "C5'", "C2", "C5", "C6", "C8")){
   #' HMQC Peaks Generation Helper Function
   #'
   #' This function allows you to convert chemical shift list to chemical shift peak table
   #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
   #' @param protons vector of proton nuclei.
   #' @param carbons vector of carbon nuclei.
   #' @export
   #' @examples
   #' create_peaks(cs)

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


 create_peaks <- function(cs, protons=c("H1'", "H2'", "H3'", "H4'", "H5'", "H5''", "H2", "H5", "H6", "H8"), carbons=c("C1'", "C2'", "C3'", "C4'", "C5'", "C5'", "C2", "C5", "C6", "C8"), grouping=c("model", "resid")){
   #' HMQC Peaks Generation Function
   #'
   #' This function allows you to convert chemical shift list to chemical shift peak table
   #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
   #' @param protons vector of proton nuclei.
   #' @param carbons vector of carbon nuclei.
   #' @param grouping variables used to group data
   #' @export
   #' @examples
   #' create_peaks(cs)

   plyr::ddply(.data = cs, .variables = grouping, .fun = create_residue_peaks, protons, carbons)
 }
