% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/likelihoods.R
\name{prob_detectable_curve}
\alias{prob_detectable_curve}
\title{Probability curve of detectable Ct values.}
\usage{
prob_detectable_curve(pars, ages)
}
\arguments{
\item{pars}{Model parameters.}

\item{ages}{Vector of times since infection.}
}
\value{
Returns a gamma-distributed curve of the probabilities of detectable
viral loads.
}
\description{
This function takes the model parameters and returns a gamma-distributed curve of
the probabilities of detectable viral loads with parameters pars and ages.
}
\seealso{
Other likelihood functions: 
\code{\link{likelihood_R}()},
\code{\link{likelihood_cpp_wrapper}()},
\code{\link{likelihood_detectable}()},
\code{\link{likelihood_pos_R}()},
\code{\link{p_a}()},
\code{\link{prop_detectable_single}()},
\code{\link{prop_detectable}()}
}
\author{
James Hay, \email{jhay@hsph.harvard.edu}
}
\concept{likelihood functions}
