
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .C calls */
extern void adajoint_chr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void artp2_chr(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void artp2_main(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void artp2_select_genes(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void check_nthread(void *);
extern void check_os(void *);
extern void ReadBED(void *, void *, void *, void *, void *, void *);
extern void test_c();

static const R_CMethodDef CEntries[] = {
  {"adajoint_chr",       (DL_FUNC) &adajoint_chr,       20},
  {"artp2_chr",          (DL_FUNC) &artp2_chr,          20},
  {"artp2_main",         (DL_FUNC) &artp2_main,         11},
  {"artp2_select_genes", (DL_FUNC) &artp2_select_genes, 12},
  {"check_nthread",      (DL_FUNC) &check_nthread,       1},
  {"check_os",           (DL_FUNC) &check_os,            1},
  {"ReadBED",            (DL_FUNC) &ReadBED,             6},
  {"test_c",             (DL_FUNC) &test_c,              0},
  {NULL, NULL, 0}
};

void R_init_ARTP2(DllInfo *dll)
{
  R_registerRoutines(dll, CEntries, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
