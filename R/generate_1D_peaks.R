create_1D_specturm <- function(cs, bw = 0.01, proton = TRUE){
  #' 1D Spectrum Generation Function
  #'
  #' This function generate a 1D NMR spectrum for chemical shift
  #' @param cs input chemical shift dataframe. Should contain field: model, resid, nucleus, weight, and predCS
  #' @param bw bandwidth for density calculation
  #' @param proton generate spectrum for proton nuclei; if not, will generate spectrum for carbon.
  #' @export
  #' @examples
  #' create_1D_specturm(cs, proton = FALSE)
  if(proton){
    cs <- cs[grepl("H", cs$nucleus),]$predCS
    xlim <- c(15, 2)
  } else {
    cs <- cs[grepl("C", cs$nucleus),]$predCS
    xlim <- c(160, 70)
  }
  density <- density(cs, bw = bw)
  plot(density, xlim = xlim)
  return(density)
}
