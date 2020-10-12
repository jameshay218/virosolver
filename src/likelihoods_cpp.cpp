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
  for(int i = 0; i < vl_ages.size(); ++i){
    vl_ages[i] = viral_load_func_single_cpp(tshift,desired_mode, t_switch,viral_peak,
                                            obs_sd, level_switch,true_0, yintercept,
                                            lod,wane_rate, wane_rate2,growth_rate,
                                            ages[i],true);
    renormalizes[i] = pgumbel_jh(yintercept, vl_ages[i],obs_sd);
  }
  //return vl_ages;
  NumericVector prob_detectable_dat(ages.size());
  double prob_undetectable=0;
  for(int i = 0; i < prob_detectable_dat.size(); ++i){
    prob_detectable_dat[i] = prop_detectable_cpp(ages[i], vl_ages[i],obs_sd, yintercept,
                                                 t_switch1, prob_detect);
    prob_undetectable += prob_detectable_dat[i]*prob_infection[obs_time-ages[i]-1];
  }
  prob_undetectable = 1 - prob_undetectable;

  NumericVector liks_tj(obs.size());

  // For each observation
  for(int i = 0; i < obs.size(); ++i){
    if(obs[i] >= yintercept){
      liks_tj[i] = prob_undetectable;
    } else {
      for(int j = 0; j < vl_ages.size(); ++j){
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
                                            ages[i],true);
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

