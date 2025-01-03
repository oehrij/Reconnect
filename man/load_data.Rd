% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Reconnect_functions.R
\name{load_data}
\alias{load_data}
\title{load csv data files from source}
\usage{
load_data(path = NULL, pattern = NULL)
}
\arguments{
\item{path}{character string describing path where csv data are stored}

\item{pattern}{character string describing regex pattern of csv files that should be loaded}
}
\value{
a numeric vector of buffer distances to use for given alphas
}
\description{
import separate csv files as one coherent dataset
}
