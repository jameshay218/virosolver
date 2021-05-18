
[![DOI](https://zenodo.org/badge/301812162.svg)](https://zenodo.org/badge/latestdoi/301812162)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)

# A method to infer epidemic dynamics using viral load data: virosolver

This README provides instructions for setup and installation of the
`virosolver` package. A vignette is provided with a simple case study
using simulated data to re-estimate the true underlying incidence curve.

Documentation for this package is a work in progress, but all code is
working correctly. This particular version accompanies the git
repository at <https://github.com/jameshay218/virosolver_paper>, where
use cases and instructions are provided.

# Setup

This package relies on the `lazymcmc` R package, which is used for the
MCMC procedure. This is easy to do with
`devtools::install_github("jameshay218/lazymcmc")`. *However*, for many
analyses where multi-modal posteriors are suspected, a separate branch
implementing parallel tempering is needed. I’d recommend you set this up
as follows: - Install the `lazymcmc` base package using
`devtools::install_github` as above. Any time this version is used,
`library(lazymcmc)` is called. - Clone the `parallel_tempering` branch
from
[here](https://github.com/jameshay218/lazymcmc/tree/parallel_tempering).
Whenever this version is needed, then `devtools::load_all("PATH TO
LAZYMCMC PARALLEL TEMPERING REPO")` is called instead.

It is possible to install the parallel tempering branch directly with
`devtools::load_all("jameshay218/lazymcmc",ref="parallel_tempering")`,
but I prefer to load the package locally until I merge the parallel
tempering and master branches.

A number of generic R packages are also used throughout:

``` r
c("tidyverse","ggthemes","ggpubr","ggsci","data.table","patchwork",
"fitdistrplus","deSolve","lazymcmc","odin","doParallel","coda")
```

Finally, this code uses compiled code with Rcpp, so you’ll need a C++
compiler (Rtools on Windows, Xcode on Mac). See
[here](http://adv-r.had.co.nz/Rcpp.html).

# Case study vignette - simulated data

``` r
library(virosolver)
library(tidyverse)
## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.0 ──
## ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
## ✓ tibble  3.0.3     ✓ dplyr   1.0.2
## ✓ tidyr   1.1.2     ✓ stringr 1.4.0
## ✓ readr   1.3.1     ✓ forcats 0.5.0
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## x dplyr::filter() masks stats::filter()
## x dplyr::lag()    masks stats::lag()
library(patchwork)
```

## Rationale

`virosolver` takes an input data frame of Ct values with associated
sample collection dates from quantitative reverse transcription PCR
(RT-qPCR) testing, and reconstructs the incidence curve that gave rise
to those measurements. The logic is as follows:

  - There is an unobserved incidence curve that describes the generation
    of new infections over time;
  - Susceptible individuals may become infected at some point in time
    with probability equal to the per-capita incidence on each day;
  - Following infection, viral loads in the infected individual follow
    some set of predictable kinetics (ie. viral loads go up then down);
  - Individuals are sampled at random from the population (a random
    cross-sectional sample), and thus an individual’s viral load is
    measured at some unknown point (as a Ct value) in their infection
    course;
  - If we mostly measure low viral loads, then most of the individuals
    we sampled were in the late stage of their infection and incidence
    was likely declining, but if we mostly measure high viral loads,
    then most individuals were early on in their infection course and
    incidence was likely growing.

By capturing this model in a mathematical model, we can obtain a
probabilistic estimate of the underlying incidence curve having observed
a set of Ct values at some point in time.

That was a bit of a whirlwind explanation – the key idea is that you can
estimate incidence based on cross-sections of observed Ct values. This
case study explores the application of `virosolver` to SARS-CoV-2. A
full explanation can be found in the accompanying
[paper](https://doi.org/10.1101/2020.10.08.20204222).

## Data

`virosolver` expects a data frame as input data in long format, where
each row corresponds to one tested sample. There should be one column
labeled as `t`, giving the time in days a sample was taken, and one
column labeled `ct`, giving the Ct value of that tested sample. Note
that Ct values are semi-quantitative and their scale depends on the
platform/instrument used. It is assumed that all Ct values within the
data frame are on the same scale and therefore internally consistent
across time points.

``` r
data(example_ct_data)
print(head(example_ct_data %>% filter(ct < 40)))
## # A tibble: 6 x 2
##       t    ct
##   <int> <dbl>
## 1    55  21.3
## 2    55  30.2
## 3    55  20.3
## 4    55  22.7
## 5    55  24.3
## 6    55  26.4
```

``` r
## Plot only detectable Ct values
p_ct_data <- ggplot(example_ct_data %>% filter(ct < 40)) + 
  geom_violin(aes(x=t,group=t,y=ct),scale="width",fill="grey70",draw_quantiles=c(0.025,0.5,0.975)) + 
  geom_jitter(aes(x=t,y=ct),size=0.1,width=2,height=0) + 
  scale_y_continuous(trans="reverse") +
  theme_bw() +
  ylab("Ct value") +
  xlab("Observation time") +
  ggtitle("Observed Ct values < 40 (the limit of detection) over time")
p_ct_data
```

![](man/figures/unnamed-chunk-4-1.png)<!-- -->

``` r
p_detectable_data <- example_ct_data %>% 
  mutate(detect=ct < 40) %>% 
  group_by(t) %>% 
  summarize(prev=sum(detect)/n()) %>% 
  ggplot() + geom_point(aes(x=t,y=prev)) + 
  theme_bw() + 
  scale_y_continuous(limits=c(0,0.5)) + 
  ylab("Percent detectable") +
  ggtitle("Proportion of samples with Ct values < 40") +
  xlab("Observation time") 
## `summarise()` ungrouping output (override with `.groups` argument)
p_detectable_data
```

![](man/figures/unnamed-chunk-5-1.png)<!-- -->

A script to generate these data can be found in the `extdata` folder,
which can be found
[here](https://github.com/jameshay218/virosolver/tree/master/man/inst/extdata).

## Ct model

A key part of the model is the assumed viral kinetics curve. This
describes the mode and variation of Ct values on each day post
infection. This is the population-level distribution – it does not track
individual-level viral kinetics curve, so the variation about the mode
captures *all* variation arising from sampling variation,
individual-level heterogeneity etc. It is *CRUCIAL* to check the
parameters underpinning this model when applying `virosolver` to a new
dataset, as this calibration will be entirely dependent on the
population being tested and the PCR instrument used. We assume the
following Ct model for this simulation:

``` r
data(example_gp_partab)
pars <- example_gp_partab$values
names(pars) <- example_gp_partab$names

## Solve the Ct model over a range of times since infection (referred to as "ages")
test_ages <- seq(1,50,by=1)

## This gives the modal Ct value
cts <- viral_load_func(pars, test_ages)

p_ct_model <- ggplot(data.frame(ct=c(40,cts),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=ct)) + 
  scale_y_continuous(trans="reverse",
                     limits=c(40,10)) +
  theme_bw() +
  ylab("Modal Ct value") +
  xlab("Days since infection")

## Note that this model does not solve for t=0, 
## as it is always assumed that no one is detectable 0 days post infection
prop_detect <- prop_detectable(test_ages,pars, cts)
p_ct_model_detectable <- ggplot(data.frame(p=c(0,prop_detect),t=c(0,test_ages))) + 
  geom_line(aes(x=t,y=p)) + 
  theme_bw() +
  ylab("Proportion of infections still detectable") +
  xlab("Days since infection")
p_ct_model/p_ct_model_detectable
```

![](man/figures/unnamed-chunk-6-1.png)<!-- -->

An intuitive way to look at this curve is to simulate observations from
it and plot the simulated Ct values ordered by days since infection.
From this, we can see the substantial amount of variation in
measurements on each day post infection, but also that there is still
some information in the Ct values. Note also that under this model,
individuals become truly undetectable (ie. fully recovered) under some
daily Bernoulli process. Those who remain detectable for a long time
demonstrate a “narrowing” of observed Ct values, where the variation
about the mode decreases. This model captures the observation that few
individuals remain detectable for a long time, but for those who do,
they do not necessarily have Ct values that tend towards the limit of
detection.

``` r
sim_cts <- simulate_viral_loads_example(test_ages, pars,N=200)
print(head(sim_cts))
## # A tibble: 6 x 3
##     age i        ct
##   <dbl> <chr> <dbl>
## 1     1 1      31.3
## 2     1 2      36.4
## 3     1 3      40  
## 4     1 4      40  
## 5     1 5      40  
## 6     1 6      34.7
p_sim_cts_age <- ggplot(sim_cts %>% filter(ct < 40)) +
  geom_point(aes(x=age,y=ct),alpha=0.25) +
  scale_y_continuous(trans="reverse",limits=c(40,10)) +
  theme_bw() +
  ylab("Ct value") +
  xlab("Days since infection") +
  ggtitle("Simulated detectable Ct values on each day post infection")
p_sim_cts_age
```

![](man/figures/unnamed-chunk-7-1.png)<!-- -->

## Inference procedure

Now that we have our Ct data and understand the assumed viral kinetics
underpinning the model, we can get into the inference framework. We use
a Markov chain Monte Carlo framework to estimate the posterior
distributions of the free model parameters, conditional on the observed
Ct data (the likelihood) and any priors we wish to place on the Ct model
and incidence curve parameters. We’ll step through the general
principles of using the MCMC package, `lazymcmc`, define priors for key
model parameters, and then demonstrate how the model works using either
a single cross section of data or multiple cross sections.

### Parameter control table

`lazymcmc` uses a data frame (usually called `parTab`) to track model
parameters. The table allows users to fix/estimate certain parameters
and also to specific upper and lower bounds. See `?example_gp_partab`
for a bit more documentation, and refer to the the
[`lazymcmc`](https://github.com/jameshay218/lazymcmc) vignettes for more
detail.

``` r
data(example_gp_partab)
head(example_gp_partab)
##       values        names fixed lower_bound upper_bound steps lower_start
## 1  0.5000000 overall_prob     0           0           1   0.1         0.0
## 2  0.0000000       tshift     1           0           3   0.1         0.0
## 3  5.0000000 desired_mode     1           0           7   0.1         0.0
## 4 19.7359875   viral_peak     0           0          40   0.1        15.0
## 5  5.0000000       obs_sd     0           0          25   0.1         1.0
## 6  0.7888288       sd_mod     1           0           1   0.1         0.4
##   upper_start
## 1         1.0
## 2        10.0
## 3        10.0
## 4        25.0
## 5        10.0
## 6         0.6
## Illustration -- set the `viral_peak` parameter to be estimated during the procedure, and the `intercept` parameter to be fixed
example_gp_partab <- example_gp_partab %>% filter(names == "viral_peak") %>% mutate(fixed=0)
example_gp_partab <- example_gp_partab %>% filter(names == "intercept") %>% mutate(fixed=1)
```

### Priors

We need to specify priors on all estimated model parameters. We use
informative priors for the Ct model, as we need to constrain its shape
somewhat to ensure identifiability of the incidence curve. We can use
less informative priors for the epidemic model, as that’s what we’re
most interested in estimating. Here, we use the example parameter table
to find the prior *means* for each model parameter, and then we set the
prior *standard deviations*. We then define a function to be used later,
which takes a vector of parameters (with names corresponding to entries
in the parameter table), which returns a single log prior probability.

``` r
## Read in the SEIR model parameter control table
data(example_seir_partab)
## Pull out the current values for each parameter, and set these as the prior means
means <- example_seir_partab$values
names(means) <- example_seir_partab$names
## Set standard deviations of prior distribution
sds_seir <- c("R0"=0.6,"obs_sd"=0.5,"viral_peak"=2,
         "wane_rate2"=1,"t_switch"=3,"level_switch"=1,
         "prob_detect"=0.03,
         "incubation"=0.25, "infectious"=0.5)

## Define a function that returns the log prior probability for a given vector of parameter
## values in `pars`, given the prior means and standard deviations set above.
prior_func_seir <- function(pars,...){
  ## Ct model priors
  obs_sd_prior <- dnorm(pars["obs_sd"], means[which(names(means) == "obs_sd")], sds_seir["obs_sd"],log=TRUE)
  r0_prior <- dlnorm(pars["R0"],log(2),sds_seir["R0"],log=TRUE)
  viral_peak_prior <- dnorm(pars["viral_peak"], means[which(names(means) == "viral_peak")], sds_seir["viral_peak"],log=TRUE)
  wane_2_prior <- dnorm(pars["wane_rate2"],means[which(names(means) == "wane_rate2")],sds_seir["wane_rate2"],log=TRUE)
  tswitch_prior <- dnorm(pars["t_switch"],means[which(names(means) == "t_switch")],sds_seir["t_switch"],log=TRUE)
  level_prior <- dnorm(pars["level_switch"],means[which(names(means) == "level_switch")],sds_seir["level_switch"],log=TRUE)
  ## Beta prior on the prob_detect parameter to ensure between 0 and 1
  beta1_mean <- means[which(names(means) == "prob_detect")]
  beta1_sd <- sds_seir["prob_detect"]
  beta_alpha <- ((1-beta1_mean)/beta1_sd^2 - 1/beta1_mean)*beta1_mean^2
  beta_beta <- beta_alpha*(1/beta1_mean - 1)
  beta_prior <- dbeta(pars["prob_detect"],beta_alpha,beta_beta,log=TRUE)
  
  ## SEIR model priors
  incu_prior <- dlnorm(pars["incubation"],log(means[which(names(means) == "incubation")]), sds_seir["incubation"], TRUE)
  infectious_prior <- dlnorm(pars["infectious"],log(means[which(names(means) == "infectious")]),sds_seir["infectious"],TRUE)
  
  ## Sum up
  obs_sd_prior + viral_peak_prior + 
    wane_2_prior + tswitch_prior + level_prior + beta_prior +
    incu_prior + infectious_prior + r0_prior
}
```

### Likelihood and posterior function

We need a function to give likelihood of observing a particular set of
Ct values conditional on the underlying model parameters. Part of this
is to determine the function for the probability of infection over time
(the incidence curve). `virosolver` has some prebuilt functions for
this, and here we’ll be using the SEIR model to describe incidence. We
can then use the `create_posterior_func` function to create a new
function which calculates the posterior probability of a set of
parameter values conditional on the Ct data. This function expects a
vector of model parameters with names corresponding to the parameter
control table and returns a single log posterior probability. Note that
the parameter vector needs entries for each parameter needed to solve
both the Ct model and the incidence model.

``` r
incidence_function <- solveSEIRModel_rlsoda_wrapper
data(example_seir_partab)
posterior_func <- create_posterior_func(parTab=example_seir_partab,
                                        data=example_ct_data,
                                        PRIOR_FUNC = prior_func_seir,
                                        INCIDENCE_FUNC = incidence_function,
                                        use_pos=FALSE)
posterior_func(example_seir_partab$values)
##    obs_sd 
## -12992.16
```
