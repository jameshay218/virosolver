#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
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
                                  bool convert_ct) {
  double y;

  if(obs_t <= tshift){
    y = true_0;
  } else if(obs_t > tshift & obs_t <= (desired_mode + tshift)) {
    y = growth_rate * (obs_t - tshift) + true_0;
  } else if(obs_t > (desired_mode + tshift) & obs_t <= (desired_mode + tshift + t_switch)) {
    y = viral_peak - wane_rate * (obs_t - (desired_mode+tshift));
  } else if(obs_t > (desired_mode + tshift + t_switch)) {
    y = level_switch - wane_rate2*(obs_t - (desired_mode + tshift + t_switch));
  } else {
    y = true_0;
  }
  if(convert_ct){
    y = yintercept - log2(10) * (y-lod);
  }
  return y;
}
