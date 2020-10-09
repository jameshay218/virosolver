#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <Rversion.h>

static const R_CMethodDef cMethods[] = {
  // Disabled until the rlsoda registration issue is fixed
  //{"C_SEIR_model_lsoda", (DL_FUNC) &SEIR_model_lsoda, 4, NULL},
  //{"C_SEIR_model_rlsoda", (DL_FUNC) &SEIR_model_rlsoda, 4, NULL},
  {NULL, NULL, 0, NULL}
};

void R_init_virosolver(DllInfo *dll) {
  R_registerRoutines(dll, cMethods, NULL, NULL, NULL);
  // Rich needs to work out how to get rlsoda to behave with registered routines, so this is disabled for now.
  // R_useDynamicSymbols(dll, FALSE);
  // R_useDynamicSymbols(dll, TRUE);
  //R_useDynamicSymbols(dll, TRUE);
  // R_forceSymbols(dll, FALSE);
}
