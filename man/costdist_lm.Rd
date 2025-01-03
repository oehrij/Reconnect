% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CONN_functions.R
\name{costdist_lm}
\alias{costdist_lm}
\title{Resistance distance to habitat in landscape}
\usage{
costdist_lm(
  x = NULL,
  resi = NULL,
  scale = 1,
  maxiter = 100,
  dist_type = "edge",
  cost_type = "least-cost"
)
}
\arguments{
\item{x}{A raster object where habitat patches are set to 1 and matrix area is set to NA OR a shapefile object containing habitat patches}

\item{resi}{resistance raster file, should have the same crs and extent as sf}

\item{scale}{numeric value: scale factor in the terra::costDist function - the cost distance is divided by this number}

\item{maxiter}{numeric value: maximum number of iterations used by the terra::costDist function}

\item{dist_type}{character, distance type: either "edge" or "centroid" (default: "edge").}

\item{cost_type}{character, gdistance algorithm applied to calculate distances: at the moment, only "least-cost" (default) available.}
}
\value{
A list with an sf object of the habitat patches (sfpa) and a corresponding square matrix (cosdi) where pairwise cost distances among habitat patches are described according to the above defined criteria
}
\description{
this is a wrapper function for the terra::costDist function
}
