% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{score_kendall}
\alias{score_kendall}
\title{Kendall Correlation Scoring Function}
\usage{
score_kendall(x)
}
\arguments{
\item{x}{input dataframe. Should contain field: predCS and expCS.}
}
\description{
This function computes the 1 - tau
}
\examples{
score_kendall(x)
}
