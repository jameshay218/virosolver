% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence_functions.R
\name{solveSEEIRRModel_rlsoda_wrapper}
\alias{solveSEEIRRModel_rlsoda_wrapper}
\title{Wrapper function for the solveSEIRRModel_rlsoda code, which is written in C.}
\usage{
solveSEEIRRModel_rlsoda_wrapper(pars, times)
}
\arguments{
\item{pars}{Vector of labeled SEEIRR model parameters.}

\item{times}{Time points for incidence calculation.}
}
\value{
Returns incidence values for given points in time.
}
\description{
This function returns a vector of incidence values for given times
according to SEIRR model parameters. It utilizes rlsoda, a solver for ordinary
differential equations but is not available on CRAN.
}
\examples{
FIX ME

}
\seealso{
Other incidence functions: 
\code{\link{SEEIRRModel}()},
\code{\link{detectable_SEEIRRModel}()},
\code{\link{detectable_SEIRModel}()},
\code{\link{exponential_growth_model}()},
\code{\link{gaussian_process_model}()},
\code{\link{reverse_gp_model}()},
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
