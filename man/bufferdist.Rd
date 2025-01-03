% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reconnect_functions.R
\name{bufferdist}
\alias{bufferdist}
\title{identify distance of low dispersal probability (<= 0.1 \%) for a given dispersal capacity (alpha) and a distance probability function (dispfun)}
\usage{
bufferdist(
  alpha = 500,
  dispfun = function(d, alpha) {
     exp(-(1/alpha) * abs(d))
 },
  prob = 0.001
)
}
\arguments{
\item{alpha}{numeric vector of average dispersal distance(s) in m}

\item{dispfun}{function for dispersal probability with a given distance (d) and dispersal probability (alpha)}
}
\value{
a numeric vector of buffer distances to use for given alphas
}
\description{
function to find the optimal distance (or buffer radius) to consider where dispersal probability <= 0.1 \% for a given alpha and dispersal probability function
}
