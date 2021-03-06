% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{score_probability}
\alias{score_probability}
\title{Probability Error Function (similar to the Quality Score)}
\usage{
score_probability(x)
}
\arguments{
\item{x}{input dataframe. Should contain field: weight, predCS and expCS.}
}
\description{
This function computes the joint probability of observing the exhibited errors between measured and predicted chemical shifts
}
\examples{
score_flat_chi2(x)
}
