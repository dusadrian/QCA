#ifndef QCA_FINDMIN_SCP_H
#define QCA_FINDMIN_SCP_H

#include "qca_rinternals.h"
#include "qca_rinternals.h"

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
