#include <R.h>
#include <Rinternals.h>

void SEIR_switch_rlsoda(double t, double *y, double *ydot, void *data)
{
  double *parms = (double*)data;
 // 3 betas
  double beta1 = parms[0], beta2=parms[1], beta3=parms[2], sigma = parms[3], gamma = parms[4];
  
  // 2 switch points between the betas
  double t_switch1=parms[5], t_switch2=parms[6];
  
  // Find the beta for this time period
  double beta;
  if(t < t_switch1){
    beta = beta1;
  } else if(t > t_switch1 && t < t_switch2) {
    beta = beta2;
  } else {
    beta = beta3;
  }
  
  double S = y[0]; // Susceptible humans
  double E = y[1]; // Exposed humans
  double I = y[2]; // Infected humans
  double R = y[3]; // Recovered humans
  
  double N = S + E + I + R;
  
  ydot[0]= -beta*I*S/N;
  ydot[1]= beta*I*S/N - sigma*E;
  ydot[2]= sigma*E - gamma*I;
  ydot[3]= gamma*I;
  ydot[4] = beta*I*S/N;
}
