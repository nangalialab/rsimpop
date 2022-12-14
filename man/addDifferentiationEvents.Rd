% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim_pop.R
\name{addDifferentiationEvents}
\alias{addDifferentiationEvents}
\title{Adds a differentiation events to a simpop. Can either specify relevant tips or randomly generate.}
\usage{
addDifferentiationEvents(
  tree,
  cfg,
  newCellType,
  idx = NULL,
  nEvent = -1,
  currentCompartment = -1
)
}
\arguments{
\item{tree}{\itemize{
\item simpop object.
}}

\item{cfg}{\itemize{
\item List. simpop config.
}}

\item{newCellType}{\itemize{
\item character
}}

\item{idx}{\itemize{
\item integer/vector. Index of the tip where to add the differentiation event.
}}

\item{currentCompartment}{\itemize{
\item integer. ID (value) of compartment where the differentiation event will be added to.
}}

\item{nEvents}{\itemize{
\item integer. the number of differentiation events to be added
}}
}
\value{
simpop object
}
\description{
Adds a differentiation events to a simpop. Can either specify relevant tips or randomly generate.
}
