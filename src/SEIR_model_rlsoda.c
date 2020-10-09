#include <R.h>
#include <Rinternals.h>

void SEIR_model_rlsoda(double t, double *y, double *ydot, void *data)
{
  double *parms = (double*)data;
  double beta = parms[0], sigma = parms[1], gamma = parms[2];

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
