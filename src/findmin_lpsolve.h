#ifndef QCA_FINDMIN_LPSOLVE_H
#define QCA_FINDMIN_LPSOLVE_H

#include <R_ext/Boolean.h>
#include <Rinternals.h>

SEXP C_findminLpSolveInternal(SEXP chart);
Rboolean solvePIchart_lpsolve(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
);

#endif
