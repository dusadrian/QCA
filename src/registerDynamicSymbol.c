#include "qca_rinternals.h"
#include "qca_rinternals.h"
#include <R_ext/Rdynload.h>
#include "solvePIchart_gurobi.h"

extern SEXP C_findminLpSolveInternal(SEXP chart);
extern SEXP C_findminScpInternal(SEXP chart);
extern SEXP C_findminLagrangianInfo(SEXP chart);
extern SEXP C_getScpProfile(void);
extern SEXP C_resetScpProfile(void);
extern SEXP C_getSA(SEXP solution_list, SEXP expressions, SEXP noflevels, SEXP mbaseexpr,
                    SEXP inputt, SEXP mbaseplus, SEXP mbase);

static const R_CallMethodDef CallEntries[] = {
  {"C_findminLpSolveInternal", (DL_FUNC) &C_findminLpSolveInternal, 1},
  {"C_findminScpInternal", (DL_FUNC) &C_findminScpInternal, 1},
  {"C_findminLagrangianInfo", (DL_FUNC) &C_findminLagrangianInfo, 1},
  {"C_getScpProfile", (DL_FUNC) &C_getScpProfile, 0},
  {"C_resetScpProfile", (DL_FUNC) &C_resetScpProfile, 0},
  {"C_getSA", (DL_FUNC) &C_getSA, 7},
  {NULL, NULL, 0}
};

void R_init_QCA(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, TRUE);
}

void R_unload_QCA(DllInfo* info) {
  (void) info;
  gurobi_release_env();
}
