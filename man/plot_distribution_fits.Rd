% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots.R
\name{plot_distribution_fits}
\alias{plot_distribution_fits}
\title{Plot distribution fits}
\usage{
plot_distribution_fits(
  chain,
  obs_dat,
  MODEL_FUNC,
  nsamps = 100,
  pos_only = TRUE
)
}
\arguments{
\item{chain}{A dataframe containing the MCMC samples}

\item{obs_dat}{A dataframe containing observed Ct values and time of sample
collection. NULL by default.}

\item{MODEL_FUNC}{Function that expects a vector of model parameters with names
corresponding to the parameter control table and returns a single log posterior probability}

\item{nsamps}{Number of samples. Defaults to 100.}

\item{pos_only}{pos_only flag uses only Ct values below the limit of detection. Defaults to TRUE.}
}
\value{
Returns two stacked ggplots.
}
\description{
Plot predicted Ct distribution fits from model
}
\examples{
\dontrun {
model_func_gp <- create_posterior_func(par_tab,example_ct_data,NULL,incidence_function,"model")
p_distribution_fit_gp <- plot_distribution_fits(chain_comb, example_ct_data, model_func_gp,100,pos_only=FALSE)
}

}
\seealso{
Other plots: 
\code{\link{plot_prob_infection}()},
\code{\link{predicted_distribution_fits}()}
}
\author{
James Hay, \email{jhay@hsph.harvard.edu}
}
\concept{plots}
