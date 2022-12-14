% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_pop.R
\name{combine_simpops}
\alias{combine_simpops}
\title{Concatenates timestamps and population size}
\usage{
combine_simpops(simpop1, simpop2)
}
\arguments{
\item{simpop1}{simpop from initial run}

\item{simpop2}{simpop from next run}
}
\value{
simpop/phylo - with branch lengths defined as above (sometimes stochastic)
}
\description{
Concatenates timestamps and population size and drivercount

DEPRECATED. This function is now automatically called by sim_pop
}
