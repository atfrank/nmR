% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/selection.R
\name{fitness}
\alias{fitness}
\title{GA Fitness Function}
\usage{
fitness(p, mask = NULL, weights = 1, alpha = 10)
}
\arguments{
\item{p}{GA parameters (vector)}

\item{mask}{mask certain parameters (vector)}

\item{alpha}{coefficient to L1 regularization term}
}
\description{
This function compute ensemble-averaged of X using w
}
\examples{
fitness(p)
}
\author{
Aaron T. Frank
}
