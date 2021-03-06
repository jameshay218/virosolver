% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/viral_load_functions.R
\name{pred_dist_wrapper}
\alias{pred_dist_wrapper}
\title{This function calculates predicted Ct value densities for each observed
time and returns a tibble of Ct values, predicted Ct value densities,
and observed times.}
\usage{
pred_dist_wrapper(
  test_cts,
  obs_times,
  ages,
  pars,
  prob_infection,
  symptom_surveillance = FALSE
)
}
\arguments{
\item{test_cts}{Vector of Ct values.}

\item{obs_times}{Vector of times at which data was collected (cross sections).}

\item{ages}{Vector of ages (time since infection).}

\item{pars}{Model parameters.}

\item{prob_infection}{Vector of the probability of infection.}

\item{symptom_surveillance}{Boolean if TRUE, then takes arguments from pars and generates the Ct distribution assuming symptom-based surveillance}
}
\value{
Tibble containing Ct values, densities, and times
}
\description{
This function calculates predicted Ct value densities for each observed
time and returns a tibble of Ct values, predicted Ct value densities,
and observed times.
}
\examples{
FIX ME

}
\seealso{
Other viral load functions: 
\code{\link{simulate_viral_loads_example}()},
\code{\link{viral_load_func}()}
}
\author{
James Hay, \email{jhay@hsph.harvard.edu}
}
\concept{viral load functions}
