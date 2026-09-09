#include "scp_relaxation.h"

#include <math.h>
#include <limits.h>
#include <R_ext/RS.h>
#include <string.h>

static int row_word_relax(int row) {
    return row >> 6;
}

static unsigned long long row_mask_relax(int row) {
    return 1ULL << (row & 63);
}

static int row_is_uncovered_relax(const qca_scp_state *state, int row) {
    return (state->uncovered[row_word_relax(row)] & row_mask_relax(row)) != 0ULL;
}

static int active_row_support_relax(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int row
) {
    int support = 0;
    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
        int col = problem->row_cols[i];
        if (state->active[col]) {
            support++;
        }
    }
    return support;
}

static double row_priority_relax(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int row
) {
    double score = 0.0;
    int count = 0;

    if (problem->branch_priority == NULL) {
        return 0.0;
    }

    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
        int col = problem->row_cols[i];
        if (!state->active[col]) {
            continue;
        }
        if (problem->branch_priority[col] < 0.0) {
            score += -problem->branch_priority[col];
        }
        count++;
    }

    if (count == 0) {
        return 0.0;
    }
    return score / (double) count;
}

static void sort_rows_by_support(
    int *rows,
    int nrows,
    const int *support,
    const double *priority,
    int descending
) {
    for (int i = 1; i < nrows; i++) {
        int row = rows[i];
        int j = i - 1;
        while (j >= 0) {
            int prev = rows[j];
            int move = 0;

            if (!descending) {
                if (support[prev] > support[row]) {
                    move = 1;
                }
                else if (support[prev] == support[row] && priority[prev] < priority[row]) {
                    move = 1;
                }
            }
            else {
                if (support[prev] < support[row]) {
                    move = 1;
                }
                else if (support[prev] == support[row] && priority[prev] < priority[row]) {
                    move = 1;
                }
            }

            if (!move) {
                break;
            }
            rows[j + 1] = prev;
            j--;
        }
        rows[j + 1] = row;
    }
}

static double dual_bound_for_order(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    const int *order,
    int nrows,
    double *slack_out
) {
    double *slack = (double *) R_Calloc((size_t) problem->nc, double);
    double bound = 0.0;

    if (slack == NULL) {
        return 0.0;
    }

    for (int col = 0; col < problem->nc; col++) {
        slack[col] = state->active[col] ? 1.0 : 0.0;
    }

    for (int idx = 0; idx < nrows; idx++) {
        int row = order[idx];
        double delta = HUGE_VAL;
        int seen = 0;

        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (!state->active[col]) {
                continue;
            }
            if (slack[col] < delta) {
                delta = slack[col];
            }
            seen = 1;
        }

        if (!seen || delta <= 1e-12) {
            continue;
        }

        bound += delta;
        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (state->active[col]) {
                slack[col] -= delta;
                if (slack[col] < 0.0) {
                    slack[col] = 0.0;
                }
            }
        }
    }

    if (slack_out != NULL) {
        memcpy(slack_out, slack, (size_t) problem->nc * sizeof(double));
    }

    R_Free(slack);
    return bound;
}

double qca_scp_relaxation_dual_info(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    double *slack_out
) {
    int *rows = NULL;
    int *support = NULL;
    double *priority = NULL;
    double *slack_tmp = NULL;
    int nrows = 0;
    double best = 0.0;

    rows = (int *) R_Calloc((size_t) problem->nr, int);
    support = (int *) R_Calloc((size_t) problem->nr, int);
    priority = (double *) R_Calloc((size_t) problem->nr, double);
    slack_tmp = (double *) R_Calloc((size_t) problem->nc, double);

    if (rows == NULL || support == NULL || priority == NULL || slack_tmp == NULL) {
        if (rows) R_Free(rows);
        if (support) R_Free(support);
        if (priority) R_Free(priority);
        if (slack_tmp) R_Free(slack_tmp);
        return 0;
    }

    for (int row = 0; row < problem->nr; row++) {
        if (!row_is_uncovered_relax(state, row)) {
            continue;
        }
        rows[nrows++] = row;
        support[row] = active_row_support_relax(problem, state, row);
        priority[row] = row_priority_relax(problem, state, row);
    }

    if (nrows == 0) {
        R_Free(rows);
        R_Free(support);
        R_Free(priority);
        R_Free(slack_tmp);
        if (slack_out != NULL) {
            memset(slack_out, 0, (size_t) problem->nc * sizeof(double));
        }
        return 0;
    }

    sort_rows_by_support(rows, nrows, support, priority, 0);
    best = dual_bound_for_order(problem, state, rows, nrows, slack_out);

    sort_rows_by_support(rows, nrows, support, priority, 1);
    {
        double alt = dual_bound_for_order(problem, state, rows, nrows, slack_tmp);
        if (alt > best) {
            best = alt;
            if (slack_out != NULL) {
                memcpy(slack_out, slack_tmp, (size_t) problem->nc * sizeof(double));
            }
        }
    }

    R_Free(rows);
    R_Free(support);
    R_Free(priority);
    R_Free(slack_tmp);

    return best;
}

int qca_scp_relaxation_dual_lb(
    const qca_scp_problem *problem,
    const qca_scp_state *state
) {
    return (int) ceil(qca_scp_relaxation_dual_info(problem, state, NULL) - 1e-12);
}

/*
For the residual binary cover, any y >= 0 gives
  L(y) = sum_r y[r] + sum_j min(0, 1 - sum_r A[r,j] y[r]).
Selecting j increases this lower bound by max(0, reduced_cost[j]).
All proof arithmetic below uses integers scaled by 65536; floating-point
subgradient steps only choose the next nonnegative integer multipliers.
Covered/deleted rows contribute zero, and excluded columns never contribute.
*/
int qca_scp_lagrangian_lb(
    const qca_scp_problem *problem, qca_scp_state *state,
    int64_t *dual, int target, int iterations,
    int64_t *costs, int *gradient, int64_t *best_dual,
    int *fixed_columns, int *iterations_used
) {
    const int64_t scale = QCA_SCP_DUAL_SCALE;
    const int64_t cutoff = (int64_t)target * scale;
    int64_t best = 0;
    *fixed_columns = 0;
    *iterations_used = 0;
    memcpy(best_dual, dual, (size_t)problem->nr * sizeof(int64_t));

    for (int iteration = 0; iteration < iterations; ++iteration) {
        ++*iterations_used;
        int64_t bound = 0;
        for (int col = 0; col < problem->nc; ++col) costs[col] = scale;
        for (int row = 0; row < problem->nr; ++row) {
            if (!row_is_uncovered_relax(state, row)) {
                dual[row] = 0;
                gradient[row] = 0;
                continue;
            }
            gradient[row] = 1;
            bound += dual[row];
            for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; ++i) {
                int col = problem->row_cols[i];
                if (state->active[col]) costs[col] -= dual[row];
            }
        }
        for (int col = 0; col < problem->nc; ++col) {
            if (state->active[col] && costs[col] < 0) {
                /* Impractically large charts still cannot overflow a proof. */
                if (bound < INT64_MIN - costs[col]) return 0;
                bound += costs[col];
            }
        }
        if (bound > best) {
            best = bound;
            memcpy(best_dual, dual, (size_t)problem->nr * sizeof(int64_t));
        }
        if (bound > cutoff) break;

        for (int col = 0; col < problem->nc; ++col) {
            if (state->active[col] && costs[col] > 0 &&
                bound > cutoff - costs[col]) {
                state->active[col] = 0;
                ++*fixed_columns;
            }
        }

        double norm = 0.0;
        for (int row = 0; row < problem->nr; ++row) {
            if (!row_is_uncovered_relax(state, row)) continue;
            for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; ++i) {
                int col = problem->row_cols[i];
                if (state->active[col] && costs[col] < 0) --gradient[row];
            }
            norm += (double)gradient[row] * gradient[row];
        }
        if (norm == 0.0 || iteration + 1 == iterations) break;
        double alpha = iteration < 4 ? 1.0 : 0.5;
        double step = alpha * ((double)cutoff + scale - (double)bound) / norm;
        for (int row = 0; row < problem->nr; ++row) {
            if (!row_is_uncovered_relax(state, row)) continue;
            double next = (double)dual[row] + step * gradient[row];
            dual[row] = next <= 0.0 ? 0 : next >= (double)scale
                ? scale : (int64_t)next;
        }
    }
    memcpy(dual, best_dual, (size_t)problem->nr * sizeof(int64_t));
    return (int)((best + scale - 1) / scale);
}
