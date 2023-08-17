#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static double parms[2];
#define beta parms[0]
#define gamma parms[1]


/* initializer */
void initmodSIR(void (* odeparms)(int *, double *))
{
  int N=2;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void SIR_model_lsoda (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
  double S = y[0]; // Susceptible
  double I = y[1]; // Infected
  double R = y[2]; // Recovered
  
  double N = S + I + R;
  
  ydot[0]= -beta*I*S/N;
  ydot[1]= beta*I*S/N - gamma*I;
  ydot[2]= gamma*I;
  ydot[3] = beta*I*S/N;
  
}