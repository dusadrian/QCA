#ifndef QCA_SCP_SOLVER_H
#define QCA_SCP_SOLVER_H

#include <stddef.h>

typedef struct {
    int nr;
    int nc;
    int nwords_rows;
    int *row_starts;
    int *row_cols;
    unsigned long long *col_masks;
    const double *branch_priority;
} qca_scp_problem;

typedef struct {
    unsigned long long *uncovered;
    unsigned char *active;
    unsigned char *chosen;
    int chosen_count;
} qca_scp_state;

typedef struct {
    double total_seconds;
    double reductions_seconds;
    double lower_bound_seconds;
    double greedy_seconds;
    double branching_seconds;
    unsigned long long nodes;
    unsigned long long leaves;
    unsigned long long reduction_calls;
    unsigned long long lower_bound_calls;
    unsigned long long greedy_calls;
} qca_scp_profile;

int qca_scp_solve_exact(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size
);

int qca_scp_solve_exact_with_incumbent(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size
);

void qca_scp_profile_reset(void);
qca_scp_profile qca_scp_profile_get(void);

#endif
