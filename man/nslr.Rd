% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nslr.R
\name{nslr}
\alias{nslr}
\title{Sum of Logarithmic Ranks Function}
\usage{
nslr(X)
}
\arguments{
\item{X}{vector of 0 (inactives) and 1 (actives) that was sorted based on some scores (e.g., agreement between measured and predicted shifts)}
}
\description{
This function allows you to compute the normalized sum of logarithmic ranks
}
\examples{
random_nslr(sample(c(rep(0,100),rep(1,10))))
}
