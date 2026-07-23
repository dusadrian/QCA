#ifndef QCA_COVER_VALIDATION_H
#define QCA_COVER_VALIDATION_H

#include "qca_rinternals.h"

Rboolean qca_cover_is_feasible(
    const int *chart,
    int nrows,
    int ncols,
    const int *indices,
    int solmin
);

#endif
