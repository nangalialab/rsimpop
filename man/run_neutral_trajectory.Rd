% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/wrapped_sims.R
\name{run_neutral_trajectory}
\alias{run_neutral_trajectory}
\title{Runs a simple single compartment neutral simulation with a specified trajectory}
\usage{
run_neutral_trajectory(
  simpop,
  initial_division_rate,
  trajectory,
  b_verbose = 1
)
}
\arguments{
\item{simpop}{\itemize{
\item Rate of symmetric cell division during development
}}

\item{initial_division_rate}{\itemize{
\item Rate of symmetric cell division once population equilibrium is reached.
}}

\item{trajectory}{\itemize{
\item data.frame - with fields ts(timestamp in days),target_pop_size,division_rate, compartment.
}}

\item{b_verbose}{\itemize{
\item FALSE or 0 for silence. 1 for progress per year. 2 to report all changes to target population size
}}

\item{nyears}{\itemize{
\item Total number of years to run the simulation
}}
}
\value{
simpop object.
}
\description{
Runs a simple single compartment neutral simulation with a specified trajectory
}
\examples{
trajectory=data.frame(ts=365*(1:80),target_pop_size=5e4+100*(1:80),division_rate=1/(2*190))
trajectory$target_pop_size[5:10]=2*trajectory$target_pop_size[5:10]
trajectory$target_pop_size[11:15]=0.2*trajectory$target_pop_size[11:15]
sp=run_neutral_trajectory(NULL,0.5,trajectory)
plot(sp)
}
