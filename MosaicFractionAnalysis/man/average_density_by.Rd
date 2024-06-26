% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MF_Functions.R
\name{average_density_by}
\alias{average_density_by}
\title{calculates average density for each unique combination of a given set of variables}
\usage{
average_density_by(nested_data, by_var, normalized = FALSE, unfiltered = FALSE)
}
\arguments{
\item{nested_data}{nested_data object including the following columns: density, all columns in by_var}

\item{by_var}{columns by which to nest}

\item{normalized}{logical value indicating whether to average normalized counts instead of counts}

\item{unfiltered}{logical value indicating whether to average unfiltered counts instead of counts}
}
\value{
nested_density object with average_density column
}
\description{
calculates average density for each unique combination of a given set of variables
}
