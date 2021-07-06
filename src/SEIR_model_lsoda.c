#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>


static double parms[3];
#define beta parms[0]
#define sigma parms[1]
#define gamma parms[2]


/* initializer */
void initmodSEIR(void (* odeparms)(int *, double *))
{
  int N=3;
  odeparms(&N, parms);
}

/* Derivatives and 1 output variable */
void SEIR_model_lsoda (int *neq, double *t, double *y, double *ydot, double *yout, int *ip)
{
  double S = y[0]; // Susceptible
  double E = y[1]; // Exposed
  double I = y[2]; // Infected
  double R = y[3]; // Recovered
  
  double N = S + E + I + R;
  
  ydot[0]= -beta*I*S/N;
  ydot[1]= beta*I*S/N - sigma*E;
  ydot[2]= sigma*E - gamma*I;
  ydot[3]= gamma*I;
  ydot[4] = beta*I*S/N;
  
}