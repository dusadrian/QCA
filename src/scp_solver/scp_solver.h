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
    unsigned long long root_branches_total;
    unsigned long long root_branches_completed;
    unsigned long long adaptive_extensions;
} qca_scp_profile;

typedef enum {
    QCA_SCP_ERROR = -1,
    QCA_SCP_NO_SOLUTION = 0,
    QCA_SCP_SOLUTION = 1,
    QCA_SCP_LIMIT = 2
} qca_scp_result;

qca_scp_result qca_scp_solve_exact(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size
);

qca_scp_result qca_scp_solve_exact_with_incumbent(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size
);

qca_scp_result qca_scp_solve_exact_with_incumbent_bounded(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size,
    unsigned long long node_limit,
    double time_limit_seconds
);

qca_scp_result qca_scp_solve_exact_with_incumbent_adaptive(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size,
    unsigned long long probe_node_limit,
    double probe_time_limit_seconds,
    unsigned long long hard_node_limit,
    double hard_time_limit_seconds
);

void qca_scp_profile_reset(void);
qca_scp_profile qca_scp_profile_get(void);

#endif
