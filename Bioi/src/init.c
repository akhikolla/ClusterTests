#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _Bioi_euclidean_linker_cpp(SEXP, SEXP, SEXP);
extern SEXP _Bioi_find_min_dists_cpp(SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
    {"_Bioi_euclidean_linker_cpp", (DL_FUNC) &_Bioi_euclidean_linker_cpp, 3},
    {"_Bioi_find_min_dists_cpp",   (DL_FUNC) &_Bioi_find_min_dists_cpp,   2},
    {NULL, NULL, 0}
};

void R_init_Bioi(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
