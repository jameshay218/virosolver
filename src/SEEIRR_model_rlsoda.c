#include <R.h>
#include <Rinternals.h>

void SEEIRR_model_rlsoda(double t, double *y, double *ydot, void *data)
{
  double *parms = (double*)data;
  double beta = parms[0], sigma = parms[1], alpha = parms[2], gamma = parms[3], omega = parms[4];

  double S = y[0]; // Susceptible
  double E1 = y[1]; // Infected, not detectable, not infectious
  double E2 = y[2]; // Infected, detectable, not infectious
  double I = y[3]; // Infected, detectable, infectious. Assume that if infectious, you are detectable by PCR
  double R1 = y[4]; // Recovered, detectable, not infectious
  double R2 = y[5]; // Recovered, not detectable, not infectious, immune

  double N = S + E1 + E2 + I + R1 + R2;

  ydot[0]= -beta*I*S/N; // Rate of getting infected
  ydot[1]= beta*I*S/N - sigma*E1; // Rate of becoming detectable
  ydot[2] = sigma*E1 - alpha*E2; // Rate of becoming infectious
  ydot[3] = alpha*E2 - gamma*I;// Rate of becoming not infectious
  ydot[4]= gamma*I - omega*R1; // Rate of becoming undetectable
  ydot[5] = omega*R1; // Rate of full recovery
  ydot[6] = beta*I*S/N;// Incidence
}
