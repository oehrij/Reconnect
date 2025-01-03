% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CONN_functions.R
\name{FragStats}
\alias{FragStats}
\title{Make a selection of FRAGSTATS landscape statistical metrics}
\usage{
FragStats(
  x = NULL,
  metrics = c("lsm_c_np", "lsm_c_pd", "lsm_c_ai", "lsm_c_clumpy", "lsm_c_cohesion",
    "lsm_c_nlsi", "lsm_c_mesh", "lsm_c_pladj"),
  directions = 4,
  count_boundary = FALSE,
  edgedepth = 2
)
}
\arguments{
\item{x}{habitat Raster*Layer file, should be binary if raster: 1= habitat and 0= non-habitat}

\item{metrics}{character, abbreviation of metrics, cf. landscapemetrics \href{https://r-spatialecology.github.io/landscapemetrics/}{R-package}}

\item{directions}{numeric, directions argument for directions in which cells are connected (4 or 8)}

\item{count_boundary}{A species-specific dispersal distance (e.g. average gap-crossing distance, in the units of the mdist values, default = 500)}

\item{edgedepth}{numeric, distance (measured in nr of cells) a cell has to be away from the patch boundary to be considered as core cell.}
}
\value{
A named vector of landscape statistical metrics: default: number of patches (lsm_c_np),patch density (lsm_c_pd), aggregation index (lsm_c_ai), clumpiness index (lsm_c_clumpy), patch cohesion index (lsm_c_cohesion), normalized landscape shape index (lsm_c_nlsi), effective mesh size (lsm_c_mesh), percentage of like adjacencies (lsm_c_pladj)
}
\description{
calculate basic FRAGSTATS metrics, such as patch cohesion using the landscapemetrics \href{https://r-spatialecology.github.io/landscapemetrics/}{R-package}
}
