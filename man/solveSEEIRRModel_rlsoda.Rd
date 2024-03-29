% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/incidence_functions.R
\name{solveSEEIRRModel_rlsoda}
\alias{solveSEEIRRModel_rlsoda}
\title{Solve SEEIRR model with rlsoda}
\usage{
solveSEEIRRModel_rlsoda(ts, init, pars, compatible = FALSE)
}
\arguments{
\item{ts}{Vector of time points for incidence calculation.}

\item{init}{Vector of labeled initial conditions for system of equations.}

\item{pars}{Vector of labeledSEEIRR model parameters.}

\item{compatible}{Boolean variable indicating whether
results should be transposed in deSolve compatible format.
Set to FALSE by default.}
}
\value{
Returns matrix of class deSolve containing solutions for SEEIRR ODEs
for each time in ts.
}
\description{
This function solves the SEEIRR ODEs for a series of
time points and returns a matrix of the solved values. rlsoda
must be manually downloaded from GitHub.
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
\code{\link{solveSEEIRRModel_rlsoda_wrapper}()},
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
