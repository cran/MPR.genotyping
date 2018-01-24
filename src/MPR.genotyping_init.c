#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void core_localMPR(void *, void *, void *, void *, void *);
extern void core_NumRecomEvents(void *, void *, void *, void *);

static const R_CMethodDef CEntries[] = {
    {"core_localMPR",       (DL_FUNC) &core_localMPR,       5},
    {"core_NumRecomEvents", (DL_FUNC) &core_NumRecomEvents, 4},
    {NULL, NULL, 0}
};

void R_init_MPR_genotyping(DllInfo *dll)
{
    R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
