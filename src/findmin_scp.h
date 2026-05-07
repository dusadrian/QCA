#ifndef QCA_FINDMIN_SCP_H
#define QCA_FINDMIN_SCP_H

#include <R_ext/Boolean.h>
#include <Rinternals.h>

SEXP C_findminScpInternal(SEXP chart);
SEXP C_getScpProfile(void);
SEXP C_resetScpProfile(void);
Rboolean solvePIchart_scp(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
);

#endif
