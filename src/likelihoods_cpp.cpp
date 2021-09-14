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
                                  bool convert_vl);

//' @export
// [[Rcpp::export]]
double dgumbel_jh(double x, double mu, double sigma){
  double a = (x-mu)/sigma;
  double prob = (1/sigma) * exp(-(a + exp(-a)));
  return prob;
}
// [[Rcpp::export]]
double pgumbel_scale(double x, double mu, double sigma){
  // Get probability between 0 and the limit of detection, x
  double a = (x-mu)/sigma;
  double prob = exp(-exp(-a));
  return prob;
}
// [[Rcpp::export]]
double pgumbel_jh(double x, double mu, double sigma){
  // Get probability between 0 and the limit of detection, x
  double a = (x-mu)/sigma;
  double prob = exp(-exp(-a));

  double a1 = (0-mu)/sigma;
  double prob1 = exp(-exp(-a1));
  return prob - prob1;
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
    additional_prob = pow((1-prob_detect),days_loss);
  }
  double main_prob = pgumbel_scale(lod, viral_load, obs_sd)*additional_prob;
  return main_prob;
}

// [[Rcpp::export]]
NumericVector likelihood_cpp(NumericVector obs,
                             double obs_time,
                             NumericVector ages,
                             NumericVector pars,
                             NumericVector prob_infection,
                             NumericVector sd_mod_vec){
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
                                            obs_sd*sd_mod_vec[i],level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_scale(yintercept, vl_ages[i],obs_sd*sd_mod_vec[i]);
    //renormalizes[i] = pgumbel_scale(yintercept, vl_ages[i],obs_sd);
    // If still detectable from gumbel dist, then increase max age
    if(renormalizes[i] > 0) max_age++;
  }

  // For each day since infection, calculate the probability that you are still detectable
  // Conditional on the probability of infection and viral kinetics model
  // Also calculate the overall probability of remaining undetectable the whole time
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd*sd_mod_vec[i],
                                                 yintercept,t_switch1, prob_detect);
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
        liks_tj[i] += (dgumbel_jh(obs[i], vl_ages[ages[j]-1], obs_sd*sd_mod_vec[j])*
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
                             NumericVector prob_infection,
                             NumericVector sd_mod_vec){
  
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
                                            obs_sd*sd_mod_vec[i],level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_scale(yintercept, vl_ages[i],obs_sd*sd_mod_vec[i]);
    // If still detectable from gumbel dist, then increase max age
    if(renormalizes[i] > 0) max_age++;
  }

  // For each day since infection, calculate the probability that you are still detectable
  // Conditional on the probability of infection and viral kinetics model
  // Also calculate the overall probability of being infected on each day in the past and
  //   still being detectable today
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd*sd_mod_vec[i],
                                                 yintercept,t_switch1, prob_detect);
    prob_detectable_and_infected += prob_detectable_dat[i]*prob_infection[obs_time-ages[i]-1];
  }

  // For each observation, find log likelihood
  // Note, all Cts are detectable here
  for(int i = 0; i < obs.size(); ++i){
      // Across all past days:
      // i) Probability of infection that day
      // ii) Probability still detectable today
      // iii) Probability of observing Ct value, given time since infection and being detectable+infected
        for(int j = 0; j < max_age; ++j){
        liks_tj[i] += (dgumbel_jh(obs[i], vl_ages[ages[j]-1], obs_sd*sd_mod_vec[j])*
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
                            NumericVector prob_infection,
                            NumericVector sd_mod_vec){
  
  double ct_step = test_cts[1]-test_cts[0];
  
  double lod = pars["LOD"];
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double obs_sd = pars["obs_sd"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double yintercept = pars["intercept"];
  // Transformed parameters
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - yintercept)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double t_switch1 = t_switch + desired_mode + tshift;

  double prob_detect = pars["prob_detect"];

  NumericVector vl_ages(ages.size());
  NumericVector renormalizes(ages.size());
  double prob_undetectable=0;
  for(int i = 0; i < vl_ages.size(); ++i){
    vl_ages[i] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd*sd_mod_vec[i],level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],false);
    renormalizes[i] = pgumbel_scale(yintercept, vl_ages[i],obs_sd*sd_mod_vec[i]);
  }
  NumericVector prob_detectable_dat(ages.size());
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd*sd_mod_vec[i],
                                                 yintercept,t_switch1, prob_detect);
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
          pgumbel_jh(test_cts[i]+ct_step, vl_ages[ages[j]-1], obs_sd*sd_mod_vec[j])-
          pgumbel_jh(test_cts[i], vl_ages[ages[j]-1], obs_sd*sd_mod_vec[j])
                         )*
          prob_infection[obs_time-ages[j]-1]*
          prob_detectable_dat[ages[j]-1])/
            renormalizes[ages[j]-1]  ;
      }
    }
  }
  
  return density;
}


// [[Rcpp::export]]
NumericVector pred_dist_cpp_symptoms(NumericVector test_cts,
                            int max_incu_period,
                            int max_sampling_delay,
                            double obs_time,
                            NumericVector pars,
                            NumericVector prob_infection,
                            NumericVector sd_mod_vec){
  
  double ct_step = test_cts[1]-test_cts[0];
  
  // Viral kinetics parameters
  double lod = pars["LOD"];
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double obs_sd = pars["obs_sd"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double yintercept = pars["intercept"];
  double prob_detect = pars["prob_detect"];
  
  // Transformed viral kinetics parameters
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - yintercept)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double t_switch1 = t_switch + desired_mode + tshift;
  
  // Incubation period parameters
  double incu_par1 = pars["incu_par1"];
  double incu_par2 = pars["incu_par2"];
  
  // Sampling delay parameters
  double sampling_par1 = pars["sampling_par1"];
  double sampling_par2 = pars["sampling_par2"];
  
  int age_max = max_incu_period + max_sampling_delay;
  
  NumericVector vl_ages(age_max+1);
  NumericVector renormalizes(age_max+1);
  double prob_undetectable=0;
  
  // Solve model for expected viral load on each day of infection
  for(int a = 0; a < vl_ages.size(); ++a){
    vl_ages[a] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd*sd_mod_vec[a],level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            a,false);
    renormalizes[a] = pgumbel_scale(yintercept, vl_ages[a],obs_sd*sd_mod_vec[a]);
  }
  
  // Solve model for proportion detectable (including uninfected) on each day of infection
  NumericVector prob_detectable_dat(age_max+1);
  for(int a = 0; a < prob_detectable_dat.size(); ++a){
    prob_detectable_dat[a] = prop_detectable_cpp(a, vl_ages[a],obs_sd*sd_mod_vec[a],
                                                 yintercept,t_switch1, prob_detect);
    
    prob_undetectable += prob_detectable_dat[a]*prob_infection[obs_time-a];
    
    // If no one is detectable of a certain age since infection, then we will be dividing
    // 0 by renormalizes[a], so set to 1 to prevent error
    if(prob_detectable_dat[a] == 0) {
      renormalizes[a] = 1;
    }
  }
  prob_undetectable = 1 - prob_undetectable;
  
  // Find probability of each discrete incubation period
  NumericVector prob_incu_period(max_incu_period+1);
  for(int o = 0; o < prob_incu_period.size(); ++o){
    prob_incu_period[o] = (R::plnorm(o+1, incu_par1, incu_par2, true, false) - R::plnorm(o, incu_par1, incu_par2, true, false))/
      R::plnorm(max_incu_period+1, incu_par1, incu_par2, true, false);
  }
  
  // Find probability of each discrete sampling delay
  NumericVector prob_sampling_delay(max_sampling_delay+1);
  for(int d = 0; d < prob_sampling_delay.size(); ++d){
    prob_sampling_delay[d] = (R::pgamma(d+1, sampling_par1, sampling_par2, true, false) - R::pgamma(d, sampling_par1, sampling_par2, true, false))/
      R::pgamma(max_sampling_delay+1, sampling_par1, sampling_par2, true, false);
  }
  
  NumericVector density(test_cts.size());
  // For each observation
  int a;
  for(int i = 0; i < test_cts.size(); ++i){
    for(int d = 0; d < prob_sampling_delay.size(); ++d){
      for(int o = 0; o < prob_incu_period.size(); ++o){
        a = d + o; // Time since infection is days since onset + sampling delay
        if(test_cts[i] >= yintercept){
          density[i] = prob_undetectable;
        } else {
          density[i] += ((
            pgumbel_jh(test_cts[i]+ct_step, vl_ages[a], obs_sd*sd_mod_vec[a])-
              pgumbel_jh(test_cts[i], vl_ages[a], obs_sd*sd_mod_vec[a])
          )*
            prob_sampling_delay[d]* // Probability of having this sampling delay
            prob_incu_period[o]* // Probability of having this incubation period
            prob_infection[obs_time-a]*
            prob_detectable_dat[a])/
              renormalizes[a]  ;
        }
      }
    }
  }
  double renormalize_age = 0;
  for(int d = 0; d < prob_sampling_delay.size(); ++d){
    for(int o = 0; o < prob_incu_period.size(); ++o){
      a = d + o; // Time since infection is days since onset + sampling delay
      renormalize_age += prob_sampling_delay[d]* // Probability of having this sampling delay
        prob_incu_period[o]* // Probability of having this incubation period
        prob_infection[obs_time-a]*
        prob_detectable_dat[a];
      }
    }
  density = density/renormalize_age;
  
  return density;
}


// [[Rcpp::export]]
NumericVector pred_age_since_inf_symptomatic(int max_incu_period,
                                             int max_sampling_delay,
                                             double obs_time,
                                             NumericVector pars,
                                             NumericVector prob_infection,
                                             NumericVector sd_mod_vec){
  
  
  // Viral kinetics parameters
  double lod = pars["LOD"];
  double tshift = pars["tshift"];
  double desired_mode = pars["desired_mode"];
  double t_switch = pars["t_switch"];
  double viral_peak = pars["viral_peak"];
  double obs_sd = pars["obs_sd"];
  double level_switch = pars["level_switch"];
  double true_0 = pars["true_0"];
  double yintercept = pars["intercept"];
  double prob_detect = pars["prob_detect"];
  
  // Transformed viral kinetics parameters
  double wane_rate = (viral_peak - level_switch)/t_switch;
  double wane_rate2 = (level_switch - yintercept)/pars["wane_rate2"];
  double growth_rate = (viral_peak - true_0)/desired_mode;
  double t_switch1 = t_switch + desired_mode + tshift;
  
  // Incubation period parameters
  double incu_par1 = pars["incu_par1"];
  double incu_par2 = pars["incu_par2"];
  
  // Sampling delay parameters
  double sampling_par1 = pars["sampling_par1"];
  double sampling_par2 = pars["sampling_par2"];
  
  int age_max = max_incu_period + max_sampling_delay + 1;
  
  NumericVector vl_ages(age_max);
  NumericVector renormalizes(age_max);
  double prob_undetectable=0;
  
  // Solve model for expected viral load on each day of infection
  for(int a = 0; a < vl_ages.size(); ++a){
    vl_ages[a] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd*sd_mod_vec[a],level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            a,false);
    renormalizes[a] = pgumbel_scale(yintercept, vl_ages[a],obs_sd*sd_mod_vec[a]);
  }
  
  // Solve model for proportion detectable (including uninfected) on each day of infection
  NumericVector prob_detectable_dat(age_max);
  for(int a = 0; a < prob_detectable_dat.size(); ++a){
    prob_detectable_dat[a] = prop_detectable_cpp(a, vl_ages[a],obs_sd*sd_mod_vec[a],
                                                 yintercept,t_switch1, prob_detect);
    
    prob_undetectable += prob_detectable_dat[a]*prob_infection[obs_time-a];
    
    // If no one is detectable of a certain age since infection, then we will be dividing
    // 0 by renormalizes[a], so set to 1 to prevent error
    if(prob_detectable_dat[a] == 0) {
      renormalizes[a] = 1;
    }
  }
  prob_undetectable = 1 - prob_undetectable;
  
  // Find probability of each discrete incubation period
  NumericVector prob_incu_period(max_incu_period);
  for(int o = 0; o < prob_incu_period.size(); ++o){
    prob_incu_period[o] = (R::plnorm(o+1, incu_par1, incu_par2, true, false) - R::plnorm(o, incu_par1, incu_par2, true, false))/
      R::plnorm(max_incu_period+1, incu_par1, incu_par2, true, false);
  }
  
  // Find probability of each discrete sampling delay
  NumericVector prob_sampling_delay(max_sampling_delay);
  for(int d = 0; d < prob_sampling_delay.size(); ++d){
    prob_sampling_delay[d] = (R::pgamma(d+1, sampling_par1, sampling_par2, true, false) - R::pgamma(d, sampling_par1, sampling_par2, true, false))/
      R::pgamma(max_sampling_delay+1, sampling_par1, sampling_par2, true, false);
  }
  
  NumericVector density(age_max);
  int a;
  for(int d = 0; d < prob_sampling_delay.size(); ++d){
    for(int o = 0; o < prob_incu_period.size(); ++o){
      a = d + o; // Time since infection is days since onset + sampling delay
      density[a] += prob_sampling_delay[d]* // Probability of having this sampling delay
        prob_incu_period[o]* // Probability of having this incubation period
        prob_infection[obs_time-a]*
        prob_detectable_dat[a];
    }
  }
  density = density/sum(density);
  return density;
}
