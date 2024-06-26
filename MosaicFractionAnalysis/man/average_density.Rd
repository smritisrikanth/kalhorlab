% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MF_util.R
\name{average_density}
\alias{average_density}
\title{calculates average density across several densities}
\usage{
average_density(density_list, colname = "count")
}
\arguments{
\item{density_list}{list of density objects with the same set of bins}

\item{colname}{name of column to average}
}
\value{
a single density representing the average across the list
}
\description{
calculates average density across several densities
}
