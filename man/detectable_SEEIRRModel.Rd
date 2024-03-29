% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence_functions.R
\name{detectable_SEEIRRModel}
\alias{detectable_SEEIRRModel}
\title{Detectable SEEIRR model}
\usage{
detectable_SEEIRRModel(pars, times)
}
\arguments{
\item{pars}{Vector of labeled SEEIRR model parameters.}

\item{times}{Time points for prevalence calculation.}
}
\value{
Returns detectable prevalence values for given points in time.
}
\description{
This function returns a vector of detectable prevalence
values for given times according to SEEIRR model parameters.
It utilizes C_SEEIRR_model_rlsoda, a function included
in the virosolver package that is written in C (code can be found in /src/)
and rlsoda, a solver for ordinary differential equations.
}
\examples{
FIX ME

}
\seealso{
Other incidence functions: 
\code{\link{SEEIRRModel}()},
\code{\link{detectable_SEIRModel}()},
\code{\link{exponential_growth_model}()},
\code{\link{gaussian_process_model}()},
\code{\link{reverse_gp_model}()},
\code{\link{solveSEEIRRModel_rlsoda_wrapper}()},
\code{\link{solveSEEIRRModel_rlsoda}()},
\code{\link{solveSEIRModel_lsoda_wrapper}()},
\code{\link{solveSEIRModel_lsoda}()},
\code{\link{solveSEIRModel_rlsoda_wrapper}()},
\code{\link{solveSEIRModel_rlsoda}()},
\code{\link{solveSEIRswitch_rlsoda_wrapper}()},
\code{\link{solveSEIRswitch_rlsoda}()}
}
\author{
James Hay, \email{jameshay218@gmail.com}
}
\concept{incidence functions}
