#ifndef QCA_FINDMIN_LAGRANGIAN_H
#define QCA_FINDMIN_LAGRANGIAN_H

#include "qca_rinternals.h"

void solvePIchart_lagrangian(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out
);

void solvePIchart_lagrangian_prepare(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out,
    unsigned char *improving_core_out,
    int *improving_core_size_out
);

SEXP C_findminLagrangian(SEXP chart);
SEXP C_findminLagrangianInfo(SEXP chart);

#endif
