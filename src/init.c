#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>

// Declare the function
extern SEXP yame_summary_cfunc(SEXP, SEXP);

// Register the function
static const R_CallMethodDef CallEntries[] = {
    {"yame_summary_cfunc", (DL_FUNC) &yame_summary_cfunc, 2},
    {NULL, NULL, 0}
};

void R_init_MethScope(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
