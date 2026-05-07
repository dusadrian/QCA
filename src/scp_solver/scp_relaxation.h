#ifndef QCA_SCP_RELAXATION_H
#define QCA_SCP_RELAXATION_H

#include "scp_solver.h"

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
