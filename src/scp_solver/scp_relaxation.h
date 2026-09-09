#ifndef QCA_SCP_RELAXATION_H
#define QCA_SCP_RELAXATION_H

#include "scp_solver.h"
#include <stdint.h>

#define QCA_SCP_DUAL_SCALE INT64_C(65536)

/* Fixed-point arithmetic makes every returned bound and fixing exact. */
int qca_scp_lagrangian_lb(
    const qca_scp_problem *problem, qca_scp_state *state,
    int64_t *dual, int target, int iterations,
    int64_t *costs, int *gradient, int64_t *best_dual,
    int *fixed_columns, int *iterations_used
);

int qca_scp_relaxation_dual_lb(
    const qca_scp_problem *problem,
    const qca_scp_state *state
);

double qca_scp_relaxation_dual_info(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    double *slack_out
);

#endif
