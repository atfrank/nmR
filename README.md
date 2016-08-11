# nmR
Functions for the Analysis and Comparison of NMR Data. As of now, the focus is NMR-derived chemical shifts.

# Create (thanks to: https://hilaryparker.com/2014/04/29/writing-an-r-package-from-scratch/)
library(roxygen2)
library("devtools")
setwd("/Users/atfrank/Documents/GitSoftware")
create("nmR")
setwd("./nmR")
document()
setwd("..")
install("nmR")

# Update
library(roxygen2)
library("devtools")
setwd("/Users/atfrank/Documents/GitSoftware")
setwd("./nmR")
document()
setwd("..")
install("nmR")

# Install
setwd("/Users/atfrank/Documents/GitSoftware")
install("nmR")