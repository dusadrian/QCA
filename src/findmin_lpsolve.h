#ifndef QCA_FINDMIN_LPSOLVE_H
#define QCA_FINDMIN_LPSOLVE_H

#include "qca_rinternals.h"
#include "qca_rinternals.h"

SEXP C_findminLpSolveInternal(SEXP chart);
Rboolean solvePIchart_lpsolve(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
);

#endif
