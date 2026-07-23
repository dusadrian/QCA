#ifndef QCA_FINDMIN_HYBRID_H
#define QCA_FINDMIN_HYBRID_H

#include "qca_rinternals.h"

SEXP C_findminHybridInternal(SEXP chart);
SEXP C_getScpProfile(void);
SEXP C_resetScpProfile(void);
Rboolean solvePIchart_hybrid(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
);
Rboolean solvePIchart_hybrid_active(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    int *indices,
    int *solmin
);
Rboolean solvePIchart_hybrid_active_with_incumbent(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    const int *initial_indices,
    int initial_solmin,
    int *indices,
    int *solmin
);

#endif
