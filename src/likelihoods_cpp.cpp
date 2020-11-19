#include <Rcpp.h>
using namespace Rcpp;

double viral_load_func_single_cpp(double tshift,
                                  double desired_mode,
                                  double t_switch,
                                  double viral_peak,
                                  double obs_sd,
                                  double level_switch,
                                  double true_0,
                                  double yintercept,
                                  double lod,
                                  double wane_rate,
                                  double wane_rate2,
                                  double growth_rate,
                                  double obs_t,
                                  bool convert_ct);

// [[Rcpp::export]]
double dgumbel_jh(double x, double mu, double sigma){
  double a = (x-mu)/sigma;
  double prob = (1/sigma) * exp(-(a + exp(-a)));
  return prob;
}
// [[Rcpp::export]]
double pgumbel_jh(double x, double mu, double sigma){
  double a = (x-mu)/sigma;
  double prob = exp(-exp(-a));
  return prob;
}


// [[Rcpp::export]]
double prop_detectable_cpp(double a,
                           double viral_load,
                           double obs_sd,
                           double lod,
                           double t_switch,
                           double prob_detect){
  double days_loss = a - t_switch;
  double additional_prob = 1;
  if(days_loss >= 0){
    additional_prob = pow((1-prob_detect),days_loss);\
  }
  double main_prob = pgumbel_jh(lod, viral_load, obs_sd)*additional_prob;
  return main_prob;
}

// [[Rcpp::export]]
NumericVector likelihood_cpp(NumericVector obs,
                             double obs_time,
                             NumericVector ages,
                             NumericVector pars,
                             NumericVector prob_infection){
  // Parameters for viral kinetics model
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double prob_detect = pars["prob_detect"];

  // Parameters for observation process
  double obs_sd = pars["obs_sd"]; // Gumbel distributed observation variation
  double yintercept = pars["intercept"];
  double lod = pars["LOD"];

  // Transformed parameters
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - yintercept)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double t_switch1 = t_switch + desired_mode + tshift;

  // Store viral loads over time since infection
  NumericVector vl_ages(ages.size());
  // Re-normalize probabilities as given detectable, must sum to 1
  NumericVector renormalizes(ages.size());

  // Probability detectable over time
  NumericVector prob_detectable_dat(ages.size());
  double prob_undetectable=0;

  // There may be an age at which the probability of retaining detectable vl is 0 because the
  // viral load curve has waned so much. If we keep solving the likelihood up to this age,
  // we get NaN because the renormalization sum is 0 (and we divdide by this). So track that
  // oldest an infection can be and still be detectable
  int max_age = 0;

  // Likelihoods for data points
  NumericVector liks_tj(obs.size());

  // For each age since infection, calculated the expected viral load and normalization factor
  for(int i = 0; i < vl_ages.size(); ++i){
    vl_ages[i] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd, level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_jh(yintercept, vl_ages[i],obs_sd);
    // If still detectable from gumbel dist, then increase max age
    if(renormalizes[i] > 0) max_age++;
  }

  // For each day since infection, calculate the probability that you are still detectable
  // Conditional on the probability of infection and viral kinetics model
  // Also calculate the overall probability of remaining undetectable the whole time
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd, yintercept,
                                                 t_switch1, prob_detect);
    prob_undetectable += prob_detectable_dat[i]*prob_infection[obs_time-ages[i]-1];
  }
  prob_undetectable = 1 - prob_undetectable;


  // For each observation, find log likelihood
  for(int i = 0; i < obs.size(); ++i){

    // If undetectable, has the same probability
    if(obs[i] >= yintercept){
      liks_tj[i] = prob_undetectable;
    } else {
      // If detectable, across all past days:
      // i) Probability of infection that day
      // ii) Probability still detectable today
      // iii) Probability of observing Ct value, given time since infection and being detectable+infected
      //for(int j = 0; j < vl_ages.size(); ++j){
      for(int j = 0; j < max_age; ++j){
        liks_tj[i] += (dgumbel_jh(obs[i], vl_ages[ages[j]-1], obs_sd)*
          prob_infection[obs_time-ages[j]-1]*
          prob_detectable_dat[ages[j]-1])/
            renormalizes[ages[j]-1]  ;
      }
    }
    liks_tj[i] = log(liks_tj[i]);
  }
  return liks_tj;
}


// [[Rcpp::export]]
NumericVector likelihood_pos_only_cpp(NumericVector obs,
                             double obs_time,
                             NumericVector ages,
                             NumericVector pars,
                             NumericVector prob_infection){
  // Parameters for viral kinetics model
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double prob_detect = pars["prob_detect"];

  // Parameters for observation process
  double obs_sd = pars["obs_sd"]; // Gumbel distributed observation variation
  double yintercept = pars["intercept"];
  double lod = pars["LOD"];

  // Transformed parameters
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - yintercept)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double t_switch1 = t_switch + desired_mode + tshift;

  // Store viral loads over time since infection
  NumericVector vl_ages(ages.size());
  // Re-normalize probabilities as given detectable, must sum to 1
  NumericVector renormalizes(ages.size());

    // Probability detectable over time
  NumericVector prob_detectable_dat(ages.size());
  double prob_detectable_and_infected=0;

  // There may be an age at which the probability of retaining detectable vl is 0 because the
  // viral load curve has waned so much. If we keep solving the likelihood up to this age,
  // we get NaN because the renormalization sum is 0 (and we divdide by this). So track that
  // oldest an infection can be and still be detectable
  int max_age = 0;

  // Likelihoods for data points
  NumericVector liks_tj(obs.size());

  // For each age since infection, calculated the expected viral load and normalization factor
  for(int i = 0; i < vl_ages.size(); ++i){
    vl_ages[i] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd, level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_jh(yintercept, vl_ages[i],obs_sd);
    // If still detectable from gumbel dist, then increase max age
    if(renormalizes[i] > 0) max_age++;
  }

  // For each day since infection, calculate the probability that you are still detectable
  // Conditional on the probability of infection and viral kinetics model
  // Also calculate the overall probability of being infected on each day in the past and
  //   still being detectable today
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd, yintercept,
                                                 t_switch1, prob_detect);
    prob_detectable_and_infected += prob_detectable_dat[i]*prob_infection[obs_time-ages[i]-1];
  }

  // For each observation, find log likelihood
  // Note, all Cts are detectable here
  for(int i = 0; i < obs.size(); ++i){
      // Across all past days:
      // i) Probability of infection that day
      // ii) Probability still detectable today
      // iii) Probability of observing Ct value, given time since infection and being detectable+infected
      //for(int j = 0; j < vl_ages.size(); ++j){
        for(int j = 0; j < max_age; ++j){
        liks_tj[i] += (dgumbel_jh(obs[i], vl_ages[ages[j]-1], obs_sd)*
          prob_infection[obs_time-ages[j]-1]*
          prob_detectable_dat[ages[j]-1])/
            renormalizes[ages[j]-1]  ;
      }
    liks_tj[i] = log(liks_tj[i]/prob_detectable_and_infected);
  }
  return liks_tj;
}



// [[Rcpp::export]]
NumericVector pred_dist_cpp(NumericVector test_cts,
                            NumericVector ages,
                            double obs_time,
                            NumericVector pars,
                            NumericVector prob_infection){

  double lod = pars["LOD"];
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double obs_sd = pars["obs_sd"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double yintercept = pars["intercept"];
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - lod)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double prob_detect = pars["prob_detect"];

  double t_switch1 = t_switch + desired_mode + tshift;
  NumericVector vl_ages(ages.size());
  NumericVector renormalizes(ages.size());
  double prob_undetectable=0;
  for(int i = 0; i < vl_ages.size(); ++i){
    vl_ages[i] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd, level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_jh(yintercept, vl_ages[i],obs_sd);
  }
  NumericVector prob_detectable_dat(ages.size());
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd, yintercept,
                                                 t_switch1, prob_detect);
    prob_undetectable += prob_detectable_dat[i]*prob_infection[obs_time-ages[i]-1];

    // If no one is detectable of a certain age since infection, then we will be dividing
    // 0 by renormalizes[i], so set to 1 to prevent error
    if(prob_detectable_dat[i] == 0) {
      renormalizes[i] = 1;
    }
  }
  prob_undetectable = 1 - prob_undetectable;

  NumericVector density(test_cts.size());
  // For each observation
  for(int i = 0; i < test_cts.size(); ++i){
    for(int j = 0; j < vl_ages.size(); ++j){
      if(test_cts[i] >= yintercept){
        density[i] = prob_undetectable;
      } else {
        density[i] += ((
          pgumbel_jh(test_cts[i]+1, vl_ages[ages[j]-1], obs_sd)-
          pgumbel_jh(test_cts[i], vl_ages[ages[j]-1], obs_sd)
                         )*
          prob_infection[obs_time-ages[j]-1]*
          prob_detectable_dat[ages[j]-1])/
            renormalizes[ages[j]-1]  ;
      }
    }
  }
  return density;
}

