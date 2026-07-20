#include "qca_rinternals.h"
#include <R_ext/Rdynload.h>
#include "findmin_gurobi.h"

extern SEXP C_findminLpSolveInternal(SEXP chart);
extern SEXP C_findminHybridInternal(SEXP chart);
extern SEXP C_findminLagrangianInfo(SEXP chart);
extern SEXP C_findminLagrangian(SEXP chart);
extern SEXP C_gurobiRuntimeAvailable(void);
extern SEXP C_findminExact(SEXP chart);
extern SEXP C_getScpProfile(void);
extern SEXP C_resetScpProfile(void);
extern SEXP C_getSA(SEXP solution_list, SEXP expressions, SEXP noflevels, SEXP mbaseexpr,
                    SEXP inputt, SEXP mbaseplus, SEXP mbase);
extern SEXP C_solveChart(SEXP pichart, SEXP allsol, SEXP vdepth, SEXP k, SEXP maxcomb, SEXP firstmin);
extern SEXP C_getRow(SEXP input);
extern SEXP C_createMatrix(SEXP input);
extern SEXP C_superSubset(SEXP x, SEXP noflevels, SEXP fuz, SEXP vo,
                          SEXP nec, SEXP inclcut, SEXP covcut, SEXP depth);
extern SEXP C_QMC(SEXP tt, SEXP noflevels);
extern SEXP C_removeRedundants(SEXP rowno, SEXP noflevels, SEXP mbase);
extern SEXP C_findSubsets(SEXP rowno, SEXP noflevels, SEXP mbase, SEXP max);
extern SEXP C_pof(SEXP x, SEXP y, SEXP nec);
extern SEXP C_omplexity(SEXP list);
extern SEXP C_expand(SEXP mat, SEXP noflevels, SEXP partial);
extern SEXP C_simplify(SEXP mat, SEXP noflevels, SEXP partial);
extern SEXP C_Cubes(SEXP list);
extern SEXP C_getEC(SEXP dem, SEXP cexpr, SEXP csolm, SEXP pexpr, SEXP psolm, SEXP SA, SEXP noflevels);
extern SEXP C_truthTable(SEXP x, SEXP vo, SEXP tt, SEXP fuz);

static const R_CallMethodDef CallEntries[] = {
  {"C_findminLpSolveInternal", (DL_FUNC) &C_findminLpSolveInternal, 1},
  {"C_findminHybridInternal", (DL_FUNC) &C_findminHybridInternal, 1},
  {"C_findminScpInternal", (DL_FUNC) &C_findminHybridInternal, 1},
  {"C_findminLagrangianInfo", (DL_FUNC) &C_findminLagrangianInfo, 1},
  {"C_findminLagrangian", (DL_FUNC) &C_findminLagrangian, 1},
  {"C_gurobiRuntimeAvailable", (DL_FUNC) &C_gurobiRuntimeAvailable, 0},
  {"C_findminExact", (DL_FUNC) &C_findminExact, 1},
  {"C_getScpProfile", (DL_FUNC) &C_getScpProfile, 0},
  {"C_resetScpProfile", (DL_FUNC) &C_resetScpProfile, 0},
  {"C_getSA", (DL_FUNC) &C_getSA, 7},
  {"C_solveChart", (DL_FUNC) &C_solveChart, 6},
  {"C_getRow", (DL_FUNC) &C_getRow, 1},
  {"C_createMatrix", (DL_FUNC) &C_createMatrix, 1},
  {"C_superSubset", (DL_FUNC) &C_superSubset, 8},
  {"C_QMC", (DL_FUNC) &C_QMC, 2},
  {"C_removeRedundants", (DL_FUNC) &C_removeRedundants, 3},
  {"C_findSubsets", (DL_FUNC) &C_findSubsets, 4},
  {"C_pof", (DL_FUNC) &C_pof, 3},
  {"C_omplexity", (DL_FUNC) &C_omplexity, 1},
  {"C_expand", (DL_FUNC) &C_expand, 3},
  {"C_simplify", (DL_FUNC) &C_simplify, 3},
  {"C_Cubes", (DL_FUNC) &C_Cubes, 1},
  {"C_getEC", (DL_FUNC) &C_getEC, 7},
  {"C_truthTable", (DL_FUNC) &C_truthTable, 4},
  {NULL, NULL, 0}
};

void R_init_QCA(DllInfo* info) {
  R_registerRoutines(info, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(info, FALSE);
}

void R_unload_QCA(DllInfo* info) {
  (void) info;
  gurobi_release_env();
}
