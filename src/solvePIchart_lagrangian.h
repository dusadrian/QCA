#ifndef QCA_SOLVE_PICHART_LAGRANGIAN_H
#define QCA_SOLVE_PICHART_LAGRANGIAN_H

#include <Rinternals.h>

void solvePIchart_lagrangian(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin
);

SEXP C_findminLagrangian(SEXP chart);

#endif
