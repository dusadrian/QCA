/*
    Copyright (c) 2016–2026, Adrian Dusa
    All rights reserved.

    License: Academic Non-Commercial License (see LICENSE file for details).
    SPDX-License-Identifier: LicenseRef-ANCL-AdrianDusa
*/

#include "findmin_lagrangian.h"
#include "qca_threads.h"
#include "scp_solver.h"

#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

typedef enum {
    LAGR_STOP_NOT_RUN = 0,
    LAGR_STOP_OPTIMAL,
    LAGR_STOP_MAX_ITER,
    LAGR_STOP_NO_UPDATE,
    LAGR_STOP_STEP_MIN,
    LAGR_STOP_INFEASIBLE,
    LAGR_STOP_OOM,
    LAGR_STOP_NO_CANDIDATE,
    LAGR_STOP_EXACT_COMPLETION
} LagrangianStopReason;

typedef struct {
    int rows;
    int cols;
    int iterations;
    int best_ub;
    int best_lb;
    int gap;
    int pool_mode;
    double best_zlb;
    double last_zlb;
    double step_coef;
    LagrangianStopReason stop_reason;
} LagrangianStats;

#define EPS 1e-12

static LagrangianStats lagr_last_stats = {
    .rows = 0,
    .cols = 0,
    .iterations = 0,
    .best_ub = -1,
    .best_lb = INT_MIN,
    .gap = -1,
    .pool_mode = 0,
    .best_zlb = -DBL_MAX,
    .last_zlb = -DBL_MAX,
    .step_coef = 0.0,
    .stop_reason = LAGR_STOP_NOT_RUN
};

void lagrangian_reset_stats(void) {
    lagr_last_stats.rows = 0;
    lagr_last_stats.cols = 0;
    lagr_last_stats.iterations = 0;
    lagr_last_stats.best_ub = -1;
    lagr_last_stats.best_lb = INT_MIN;
    lagr_last_stats.gap = -1;
    lagr_last_stats.pool_mode = 0;
    lagr_last_stats.best_zlb = -DBL_MAX;
    lagr_last_stats.last_zlb = -DBL_MAX;
    lagr_last_stats.step_coef = 0.0;
    lagr_last_stats.stop_reason = LAGR_STOP_NOT_RUN;
}

const LagrangianStats *lagrangian_last_stats(void) {
    return &lagr_last_stats;
}

const char *lagrangian_stop_reason_name(LagrangianStopReason reason) {
    switch (reason) {
        case LAGR_STOP_OPTIMAL: return "optimal";
        case LAGR_STOP_MAX_ITER: return "max_iter";
        case LAGR_STOP_NO_UPDATE: return "no_update";
        case LAGR_STOP_STEP_MIN: return "step_min";
        case LAGR_STOP_INFEASIBLE: return "infeasible";
        case LAGR_STOP_OOM: return "oom";
        case LAGR_STOP_NO_CANDIDATE: return "no_candidate";
        case LAGR_STOP_EXACT_COMPLETION: return "exact_completion";
        case LAGR_STOP_NOT_RUN:
        default: return "not_run";
    }
}

static void lagr_stats_begin(int rows, int cols, int pool_mode) {
    lagrangian_reset_stats();
    lagr_last_stats.rows = rows;
    lagr_last_stats.cols = cols;
    lagr_last_stats.pool_mode = pool_mode;
}

static void lagr_stats_finish(
    int best_ub,
    double best_lb,
    double best_zlb,
    double last_zlb,
    double step_coef,
    int iterations,
    LagrangianStopReason stop_reason
) {
    lagr_last_stats.best_ub = best_ub;
    if (best_lb > (double)INT_MIN && best_lb < (double)INT_MAX) {
        lagr_last_stats.best_lb = (int)best_lb;
    } else {
        lagr_last_stats.best_lb = INT_MIN;
    }
    if (lagr_last_stats.best_ub >= 0 && lagr_last_stats.best_lb != INT_MIN) {
        lagr_last_stats.gap = lagr_last_stats.best_ub - lagr_last_stats.best_lb;
        if (lagr_last_stats.gap < 0) lagr_last_stats.gap = 0;
    } else {
        lagr_last_stats.gap = -1;
    }
    lagr_last_stats.best_zlb = best_zlb;
    lagr_last_stats.last_zlb = last_zlb;
    lagr_last_stats.step_coef = step_coef;
    lagr_last_stats.iterations = iterations;
    lagr_last_stats.stop_reason = stop_reason;
}

typedef struct {
    const int *pichart;
    int rows;
    int cols;
    int *counts;
} lagr_count_context;

static void lagr_count_col_rows_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    (void)worker_id;
    lagr_count_context *ctx = (lagr_count_context *)data;
    for (unsigned long long c0 = start; c0 < end; ++c0) {
        int c = (int)c0;
        int cnt = 0;
        for (int r = 0; r < ctx->rows; ++r) {
            cnt += ctx->pichart[c * ctx->rows + r];
        }
        ctx->counts[c] = cnt;
    }
}

static void lagr_count_row_cols_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    (void)worker_id;
    lagr_count_context *ctx = (lagr_count_context *)data;
    for (unsigned long long r0 = start; r0 < end; ++r0) {
        int r = (int)r0;
        int cnt = 0;
        for (int c = 0; c < ctx->cols; ++c) {
            cnt += ctx->pichart[c * ctx->rows + r];
        }
        ctx->counts[r] = cnt;
    }
}

typedef struct {
    const int *pichart;
    int rows;
    int cols;
    int **out;
} lagr_fill_context;

static void lagr_fill_rows_covered_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    (void)worker_id;
    lagr_fill_context *ctx = (lagr_fill_context *)data;
    for (unsigned long long c0 = start; c0 < end; ++c0) {
        int c = (int)c0;
        int k = 0;
        for (int r = 0; r < ctx->rows; ++r) {
            if (ctx->pichart[c * ctx->rows + r] != 0) {
                ctx->out[c][k++] = r;
            }
        }
    }
}

static void lagr_fill_cols_covering_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    (void)worker_id;
    lagr_fill_context *ctx = (lagr_fill_context *)data;
    for (unsigned long long r0 = start; r0 < end; ++r0) {
        int r = (int)r0;
        int k = 0;
        for (int c = 0; c < ctx->cols; ++c) {
            if (ctx->pichart[c * ctx->rows + r] != 0) {
                ctx->out[r][k++] = c;
            }
        }
    }
}

typedef struct {
    int **rowsCovered;
    int *rowsCoveredCount;
    const double *t;
    double *ls;
    double *partials;
    const double *col_costs;
    unsigned char *x;
} lagr_reduced_context;

static void lagr_reduced_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    lagr_reduced_context *ctx = (lagr_reduced_context *)data;
    double neg = 0.0;
    for (unsigned long long c0 = start; c0 < end; ++c0) {
        int c = (int)c0;
        double sum = 0.0;
        for (int k = 0; k < ctx->rowsCoveredCount[c]; ++k) {
            int r = ctx->rowsCovered[c][k];
            sum += ctx->t[r];
        }
        double cost = ctx->col_costs ? ctx->col_costs[c] : 1.0;
        ctx->ls[c] = cost - sum;
        if (ctx->x) ctx->x[c] = (ctx->ls[c] < 0.0) ? 1 : 0;
        if (ctx->ls[c] < 0.0) neg += ctx->ls[c];
    }
    ctx->partials[worker_id] = neg;
}

typedef struct {
    int **colsCovering;
    int *colsCoveringCount;
    const unsigned char *x;
    double *slack;
    double *partials;
} lagr_slack_context;

static void lagr_slack_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    lagr_slack_context *ctx = (lagr_slack_context *)data;
    double sum_s2 = 0.0;
    for (unsigned long long i0 = start; i0 < end; ++i0) {
        int i = (int)i0;
        int covered = 0;
        for (int k = 0; k < ctx->colsCoveringCount[i]; ++k) {
            int c = ctx->colsCovering[i][k];
            covered += ctx->x[c];
        }
        double s = 1.0 - (double)covered;
        ctx->slack[i] = s;
        sum_s2 += s * s;
    }
    ctx->partials[worker_id] = sum_s2;
}

/*
Build adjacency:
rowsCovered[c] = list of rows that column c covers
colsCovering[r] = list of columns that cover row r
Returns:
0 on success
-2 if any row has no covering column (infeasible)
-1 on out of memory.
*/

static int build_adjacency(
    const int *pichart,
    int rows, int cols,
    int ***rowsCovered,
    int **rowsCoveredCount,
    int ***colsCovering,
    int **colsCoveringCount
) {
    int **rc = NULL, *rcc = NULL;
    int **cr = NULL, *crc = NULL;

    rcc = (int*)calloc((size_t)cols, sizeof(int));
    crc = (int*)calloc((size_t)rows, sizeof(int));
    if (!rcc || !crc) goto oom;

    lagr_count_context count_cols = {
        .pichart = pichart,
        .rows = rows,
        .cols = cols,
        .counts = rcc
    };
    if (!qca_parallel_for((unsigned long long) cols, qca_default_thread_count(), lagr_count_col_rows_worker, &count_cols)) {
        goto oom;
    }

    lagr_count_context count_rows = {
        .pichart = pichart,
        .rows = rows,
        .cols = cols,
        .counts = crc
    };
    if (!qca_parallel_for((unsigned long long) rows, qca_default_thread_count(), lagr_count_row_cols_worker, &count_rows)) {
        goto oom;
    }

    for (int r = 0; r < rows; ++r) {
        if (crc[r] == 0) {
            free(rcc);
            free(crc);
            return -2;
        }
    }

    rc = (int**)calloc((size_t)cols, sizeof(int*));
    cr = (int**)calloc((size_t)rows, sizeof(int*));
    if (!rc || !cr) goto oom;

    for (int c = 0; c < cols; ++c) {
        rc[c] = (rcc[c] > 0) ? (int*)malloc((size_t)rcc[c] * sizeof(int)) : NULL;
        if (rcc[c] > 0 && !rc[c]) goto oom;
    }
    for (int r = 0; r < rows; ++r) {
        cr[r] = (crc[r] > 0) ? (int*)malloc((size_t)crc[r] * sizeof(int)) : NULL;
        if (crc[r] > 0 && !cr[r]) goto oom;
    }

    lagr_fill_context fill_rows = {
        .pichart = pichart,
        .rows = rows,
        .cols = cols,
        .out = rc
    };
    if (!qca_parallel_for((unsigned long long) cols, qca_default_thread_count(), lagr_fill_rows_covered_worker, &fill_rows)) {
        goto oom;
    }

    lagr_fill_context fill_cols = {
        .pichart = pichart,
        .rows = rows,
        .cols = cols,
        .out = cr
    };
    if (!qca_parallel_for((unsigned long long) rows, qca_default_thread_count(), lagr_fill_cols_covering_worker, &fill_cols)) {
        goto oom;
    }

    *rowsCovered = rc;
    *rowsCoveredCount = rcc;
    *colsCovering = cr;
    *colsCoveringCount = crc;
    return 0;

oom: // out of memory
    if (rc) {
        for (int c = 0; c < cols; ++c) free(rc[c]);
        free(rc);
    }
    if (cr) {
        for (int r = 0; r < rows; ++r) free(cr[r]);
        free(cr);
    }
    free(rcc);
    free(crc);
    return -1;
}

static void free_adjacency(
    int **rowsCovered,
    int *rowsCoveredCount,
    int cols,
    int **colsCovering,
    int *colsCoveringCount,
    int rows
) {
    if (rowsCovered) {
        for (int c = 0; c < cols; ++c) free(rowsCovered[c]);
        free(rowsCovered);
    }

    free(rowsCoveredCount);

    if (colsCovering) {
        for (int r = 0; r < rows; ++r) free(colsCovering[r]);
        free(colsCovering);
    }

    free(colsCoveringCount);
}

/*
Column presolve (single-solution path only): deactivate any column whose
coverage is a subset of another column's coverage, provided the dominating
column has at least equal weight. This preserves both the minimum-size
objective and the max-weight tie-break, while shrinking the instance for
the subgradient loop, the heuristics and the exact finish.
Deactivated columns get rowsCoveredCount[c] = 0 and are removed from the
colsCovering lists; column indices remain stable, so solutions still refer
to the original pichart columns.
Returns the number of deactivated columns (0 on OOM: presolve is skipped).
*/
static int presolve_dominated_columns(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *weights
) {
    if (rows <= 0 || cols <= 1) return 0;

    int words = (rows + 63) / 64;
    uint64_t *bits = (uint64_t*)calloc((size_t)cols * (size_t)words, sizeof(uint64_t));
    unsigned char *active = (unsigned char*)malloc((size_t)cols);
    if (!bits || !active) {
        free(bits);
        free(active);
        return 0;
    }

    for (int c = 0; c < cols; ++c) {
        active[c] = 1;
        uint64_t *bc = bits + (size_t)c * words;
        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
            int r = rowsCovered[c][k];
            bc[r >> 6] |= 1ULL << (r & 63);
        }
    }

    int removed = 0;
    for (int c = 0; c < cols; ++c) {
        if (!active[c] || rowsCoveredCount[c] <= 0) continue;

        /* rarest row covered by c: smallest candidate list to scan */
        int rmin = -1;
        int best = INT_MAX;
        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
            int r = rowsCovered[c][k];
            if (colsCoveringCount[r] < best) {
                best = colsCoveringCount[r];
                rmin = r;
            }
        }
        if (rmin < 0) continue;

        double wc = weights ? weights[c] : 0.0;
        const uint64_t *bc = bits + (size_t)c * words;

        for (int k = 0; k < colsCoveringCount[rmin]; ++k) {
            int d = colsCovering[rmin][k];
            if (d == c || !active[d]) continue;
            if (rowsCoveredCount[d] < rowsCoveredCount[c]) continue;

            double wd = weights ? weights[d] : 0.0;
            if (wd + EPS < wc) continue; /* keep the heavier column for the tie-break */

            if (rowsCoveredCount[d] == rowsCoveredCount[c] && fabs(wd - wc) <= EPS && d > c) {
                continue; /* identical coverage and weight: keep the lower index */
            }

            const uint64_t *bd = bits + (size_t)d * words;
            int subset = 1;
            for (int w = 0; w < words; ++w) {
                if (bc[w] & ~bd[w]) {
                    subset = 0;
                    break;
                }
            }

            if (subset) {
                active[c] = 0;
                ++removed;
                break;
            }
        }
    }

    if (removed > 0) {
        for (int c = 0; c < cols; ++c) {
            if (!active[c]) {
                free(rowsCovered[c]);
                rowsCovered[c] = NULL;
                rowsCoveredCount[c] = 0;
            }
        }
        for (int r = 0; r < rows; ++r) {
            int w = 0;
            for (int k = 0; k < colsCoveringCount[r]; ++k) {
                int cc = colsCovering[r][k];
                if (active[cc]) colsCovering[r][w++] = cc;
            }
            colsCoveringCount[r] = w;
        }
    }

    free(bits);
    free(active);
    return removed;
}

/* Redundancy pruning: remove columns that are not needed to keep all rows covered (reverse order) */
static void prune_redundancy(
    int rows,
    int **rowsCovered,
    int *rowsCoveredCount,
    int *sol, int *sol_len
) {
    int n = *sol_len;
    if (n <= 0) return;

    int *coverCount = (int*)calloc((size_t)rows, sizeof(int));
    if (!coverCount) return;

    for (int i = 0; i < n; ++i) {
        int c = sol[i];
        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
            int r = rowsCovered[c][k];
            coverCount[r]++;
        }
    }

    int write = n;
    for (int i = n - 1; i >= 0; --i) {
        int c = sol[i];
        bool removable = true;

        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
            int r = rowsCovered[c][k];
            if (coverCount[r] <= 1) {
                removable = false;
                break;
            }
        }

        if (removable) {
            for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                int r = rowsCovered[c][k];
                coverCount[r]--;
            }

            for (int j = i + 1; j < write; ++j) {
                sol[j - 1] = sol[j];
            }
            write--;
        }
    }

    *sol_len = write;
    free(coverCount);
}

/*
Greedy construction guided by Lagrangian scores
Primary: maximize newCover (minimize number of columns)
Tie 1: smaller lagr_score (more negative is better)
Tie 2: higher weight
Tie 3: smaller index
*/
static void greedy_from_lagr_scores(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights, /* weights can be NULL */
    int *sol,
    int *sol_len
) {
    (void)colsCovering;
    (void)colsCoveringCount; /* not used here, but kept for symmetry */

    bool *rowCovered = (bool*)calloc((size_t)rows, sizeof(bool));
    bool *colSelected = (bool*)calloc((size_t)cols, sizeof(bool));
    int covered = 0, out = 0;

    if (!rowCovered || !colSelected) {
        free(rowCovered);
        free(colSelected);
        *sol_len = -1;
        return;
    }

    while (covered < rows) {
        int best = -1;
        int bestNew = -1;
        double bestLS = DBL_MAX;
        double bestW = -DBL_MAX;

        for (int c = 0; c < cols; ++c) {
            if (colSelected[c]) continue;

            int newCover = 0;
            for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                int r = rowsCovered[c][k];
                if (!rowCovered[r]) newCover++;
            }
            if (newCover <= 0) continue;

            double ls = lagr_score ? lagr_score[c] : 0.0;
            double w = weights ? weights[c] : 0.0;

            bool better = false;
            if (newCover > bestNew) better = true;
            else if (newCover == bestNew) {
                if (ls < bestLS) {
                    better = true;
                } else if (fabs(ls - bestLS) <= EPS) {
                    if (w > bestW) {
                        better = true;
                    } else if (fabs(w - bestW) <= EPS) {
                        if (best == -1 || c < best) better = true;
                    }
                }
            }
            if (better) {
                best = c;
                bestNew = newCover;
                bestLS = ls;
                bestW = w;
            }
        }

        if (best == -1) { /* cannot cover remaining rows -> infeasible construct */
            free(rowCovered);
            free(colSelected);
            *sol_len = -1;
            return;
        }

        colSelected[best] = true;
        sol[out++] = best;

        for (int k = 0; k < rowsCoveredCount[best]; ++k) {
            int r = rowsCovered[best][k];
            if (!rowCovered[r]) {
                rowCovered[r] = true; covered++;
            }
        }
    }

    *sol_len = out;
    prune_redundancy(
        rows,
        rowsCovered,
        rowsCoveredCount,
        sol,
        sol_len
    );

    free(rowCovered);
    free(colSelected);
}

/*
Compute Lagrangian scores:
rc_c = 1.0 - sum_{r in rowsCovered_c} t[r]
ZLB = sum_i t_i + sum_c min(0, rc_c)
*/
static double compute_reduced_and_lb(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    const double *t,
    const double *col_costs, // can be NULL
    double *ls, /* out */
    unsigned char *x /* optional sign flags for subgradient */
) {
    double ZLB = 0.0;
    double ZLB_rows = 0.0;
    double ZLB_neg = 0.0;

    int threads = qca_default_thread_count();
    double *partials = (double*)calloc((size_t)threads, sizeof(double));
    if (!partials) return 0.0;

    for (int i = 0; i < rows; ++i) ZLB_rows += t[i];

    lagr_reduced_context reduced_ctx = {
        .rowsCovered = rowsCovered,
        .rowsCoveredCount = rowsCoveredCount,
        .t = t,
        .ls = ls,
        .partials = partials,
        .col_costs = col_costs,
        .x = x
    };
    if (!qca_parallel_for((unsigned long long) cols, threads, lagr_reduced_worker, &reduced_ctx)) {
        free(partials);
        return 0.0;
    }
    for (int i = 0; i < threads; ++i) {
        ZLB_neg += partials[i];
    }
    free(partials);

    ZLB = ZLB_rows + ZLB_neg;
    return ZLB;
}

/*
A "slack" measures how far a constraint is from being tight.
Positive slack means the constraint is not yet satisfied, zero slack means tight, negative means surplus.
In this code, a slack is the per-row violation of the covering constraint under the current Lagrangian/primal guess
- Constraint per row i: sum over columns covering i of x_c >= 1
- Slack definition: slack_i = 1 − Σ_{c covers i} x_c
- slack_i > 0: row i is under-covered (violation)
- slack_i = 0: row i is exactly satisfied
- slack_i < 0: row i is over-covered

How it’s used:
- Build a tentative primal x by taking columns with negative lagr_score (x_c = 1 if rc_c < 0, else 0).
- The vector slack is exactly the subgradient of the dual Lagrangian at the current multipliers t.
- Step size uses (UB − ZLB) / Σ slack_i^2, and multipliers update as t_i := max(0, t_i + step · slack_i).
- If all slack_i are zero, there’s no subgradient direction to update (no progress possible).

Subgradient update:
x_c = 1 if rc_c < 0 else 0 (LP solution)
slack_i = 1 - sum_{c covering i} x_c
step = step_coef * (UB - ZLB) / sum(slack^2), t := max(0, t + step*slack)
Returns 1 if step applied, 0 if slacks zero or no progress possible.
*/

static int subgradient_update(
    int rows,
    int **colsCovering,
    int *colsCoveringCount,
    const unsigned char *x,
    double UB, double ZLB,
    double *t,
    double *step_coef,
    double step_min,
    int *stuck_iter,
    int stuck_halve_period,
    double *prev_direction,
    double deflection_alpha,
    double step_contract,
    double stabilization_beta
) {
    int threads = qca_default_thread_count();
    if (!x) return 0;

    double *slack = (double*)malloc((size_t)rows * sizeof(double));
    if (!slack) return 0;

    double sum_s2 = 0.0;

    double *partials = (double*)calloc((size_t)threads, sizeof(double));
    if (!partials) {
        free(slack);
        return 0;
    }

    lagr_slack_context slack_ctx = {
        .colsCovering = colsCovering,
        .colsCoveringCount = colsCoveringCount,
        .x = x,
        .slack = slack,
        .partials = partials
    };
    if (!qca_parallel_for((unsigned long long) rows, threads, lagr_slack_worker, &slack_ctx)) {
        free(slack);
        free(partials);
        return 0;
    }
    for (int i = 0; i < threads; ++i) {
        sum_s2 += partials[i];
    }
    free(partials);

    if (deflection_alpha > 0.0 && prev_direction) {
        sum_s2 = 0.0;
        for (int i = 0; i < rows; ++i) {
            slack[i] += deflection_alpha * prev_direction[i];
            sum_s2 += slack[i] * slack[i];
        }
    }

    if (sum_s2 <= 1e-15) {
        free(slack);
        return 0;
    }

    double gap = UB - ZLB;
    if (gap <= 0.0) {
        free(slack);
        return 0;
    }

    double step = (*step_coef) * (gap / sum_s2);
    if (step <= 0.0) {
        free(slack);
        return 0;
    }

    for (int i = 0; i < rows; ++i) {
        double trial = t[i] + step * slack[i];
        if (trial < 0.0) trial = 0.0;
        double val = (stabilization_beta >= 1.0)
            ? trial
            : ((1.0 - stabilization_beta) * t[i] + stabilization_beta * trial);
        t[i] = (val > 0.0) ? val : 0.0;
    }

    if (deflection_alpha > 0.0 && prev_direction) {
        memcpy(prev_direction, slack, (size_t)rows * sizeof(double));
    }

    (*stuck_iter)++;
    if (*stuck_iter >= stuck_halve_period) {
        *stuck_iter = 0;
        *step_coef *= step_contract;
        if (*step_coef < step_min) *step_coef = step_min;
    }

    free(slack);
    return 1;
}

typedef struct {
    int rows;
    int memory;
    int count;
    int head;
    double prox_mu;
    double min_mu;
    double max_mu;
    double null_expand;
    double serious_shrink;
    double serious_fraction;
    double serious_tol;
    double momentum;
    double center_value;
    double predicted_gain;
    double *center_t;
    double *prev_center_t;
    double *cut_alpha;
    double *cut_slacks;
    int null_streak;
    int serious_streak;
    int has_trial;
} LagrangianBundleState;

static void bundle_project_simplex(double *v, int n) {
    if (!v || n <= 0) return;

    double *u = (double*)malloc((size_t)n * sizeof(double));
    if (!u) {
        double uniform = 1.0 / (double)n;
        for (int i = 0; i < n; ++i) v[i] = uniform;
        return;
    }

    memcpy(u, v, (size_t)n * sizeof(double));
    for (int i = 1; i < n; ++i) {
        double key = u[i];
        int j = i - 1;
        while (j >= 0 && u[j] < key) {
            u[j + 1] = u[j];
            --j;
        }
        u[j + 1] = key;
    }

    double cssv = 0.0;
    int rho = 0;
    for (int i = 0; i < n; ++i) {
        cssv += u[i];
        double theta = (cssv - 1.0) / (double)(i + 1);
        if (u[i] - theta > 0.0) rho = i + 1;
    }

    cssv = 0.0;
    for (int i = 0; i < rho; ++i) cssv += u[i];
    double theta = rho > 0 ? (cssv - 1.0) / (double)rho : 0.0;

    for (int i = 0; i < n; ++i) {
        double projected = v[i] - theta;
        v[i] = projected > 0.0 ? projected : 0.0;
    }

    free(u);
}

static int bundle_update(
    int rows,
    int **colsCovering,
    int *colsCoveringCount,
    const unsigned char *x,
    double UB,
    double ZLB,
    double *t,
    LagrangianBundleState *bs
) {
    if (!x || !bs || bs->memory <= 0 || !bs->center_t || !bs->cut_alpha || !bs->cut_slacks) {
        return 0;
    }

    int threads = qca_default_thread_count();

    double *slack = (double*)malloc((size_t)rows * sizeof(double));
    if (!slack) return 0;

    double *partials = (double*)calloc((size_t)threads, sizeof(double));
    if (!partials) {
        free(slack);
        return 0;
    }

    lagr_slack_context slack_ctx = {
        .colsCovering = colsCovering,
        .colsCoveringCount = colsCoveringCount,
        .x = x,
        .slack = slack,
        .partials = partials
    };
    if (!qca_parallel_for((unsigned long long) rows, threads, lagr_slack_worker, &slack_ctx)) {
        free(slack);
        free(partials);
        return 0;
    }

    double sum_s2 = 0.0;
    for (int i = 0; i < threads; ++i) sum_s2 += partials[i];
    free(partials);

    if (sum_s2 <= 1e-15 || UB - ZLB <= 0.0) {
        free(slack);
        return 0;
    }

    if (bs->center_value <= -DBL_MAX / 2.0) {
        memcpy(bs->center_t, t, (size_t)rows * sizeof(double));
        if (bs->prev_center_t) {
            memcpy(bs->prev_center_t, t, (size_t)rows * sizeof(double));
        }
        bs->center_value = ZLB;
        bs->has_trial = 0;
    } else if (bs->has_trial) {
        double actual_gain = ZLB - bs->center_value;
        if (actual_gain > bs->serious_tol) {
            if (bs->prev_center_t) {
                memcpy(bs->prev_center_t, bs->center_t, (size_t)rows * sizeof(double));
            }
            memcpy(bs->center_t, t, (size_t)rows * sizeof(double));
            bs->center_value = ZLB;
            bs->prox_mu = fmax(bs->min_mu, bs->prox_mu * bs->serious_shrink);
            bs->serious_streak++;
            bs->null_streak = 0;
        } else {
            bs->prox_mu = fmin(bs->max_mu, bs->prox_mu * bs->null_expand);
            bs->null_streak++;
            bs->serious_streak = 0;
        }
        bs->has_trial = 0;
    } else if (ZLB > bs->center_value + bs->serious_tol) {
        if (bs->prev_center_t) {
            memcpy(bs->prev_center_t, bs->center_t, (size_t)rows * sizeof(double));
        }
        memcpy(bs->center_t, t, (size_t)rows * sizeof(double));
        bs->center_value = ZLB;
        bs->serious_streak++;
        bs->null_streak = 0;
    } else {
        bs->prox_mu = fmin(bs->max_mu, bs->prox_mu * bs->null_expand);
        bs->null_streak++;
        bs->serious_streak = 0;
    }

    double dot_gt = 0.0;
    for (int r = 0; r < rows; ++r) dot_gt += slack[r] * t[r];
    double alpha = -ZLB + dot_gt;

    if (bs->count < bs->memory) {
        int slot = (bs->head + bs->count) % bs->memory;
        memcpy(bs->cut_slacks + (size_t)slot * (size_t)rows, slack, (size_t)rows * sizeof(double));
        bs->cut_alpha[slot] = alpha;
        bs->count++;
    } else {
        int slot = bs->head;
        memcpy(bs->cut_slacks + (size_t)slot * (size_t)rows, slack, (size_t)rows * sizeof(double));
        bs->cut_alpha[slot] = alpha;
        bs->head = (bs->head + 1) % bs->memory;
    }

    int m = bs->count;
    if (m <= 0 || (bs->prox_mu >= bs->max_mu - 1e-12 && bs->null_streak >= bs->memory * 8)) {
        free(slack);
        return 0;
    }

    double *base_t = (double*)malloc((size_t)rows * sizeof(double));
    double *linear = (double*)malloc((size_t)m * sizeof(double));
    double *gram = (double*)calloc((size_t)m * (size_t)m, sizeof(double));
    double *lambda = (double*)malloc((size_t)m * sizeof(double));
    double *gradient = (double*)malloc((size_t)m * sizeof(double));
    double *gbar = (double*)calloc((size_t)rows, sizeof(double));
    if (!base_t || !linear || !gram || !lambda || !gradient || !gbar) {
        free(base_t);
        free(linear);
        free(gram);
        free(lambda);
        free(gradient);
        free(gbar);
        free(slack);
        return 0;
    }

    for (int r = 0; r < rows; ++r) {
        double base = bs->center_t[r];
        if (bs->momentum > 0.0 && bs->prev_center_t && bs->serious_streak > 0) {
            base += bs->momentum * (bs->center_t[r] - bs->prev_center_t[r]);
        }
        base_t[r] = base > 0.0 ? base : 0.0;
    }

    int start = bs->head;
    for (int i = 0; i < m; ++i) {
        int slot_i = (start + i) % bs->memory;
        const double *g_i = bs->cut_slacks + (size_t)slot_i * (size_t)rows;
        double dot_center = 0.0;
        for (int r = 0; r < rows; ++r) dot_center += g_i[r] * base_t[r];
        linear[i] = bs->cut_alpha[slot_i] - dot_center;

        for (int j = i; j < m; ++j) {
            int slot_j = (start + j) % bs->memory;
            const double *g_j = bs->cut_slacks + (size_t)slot_j * (size_t)rows;
            double dot = 0.0;
            for (int r = 0; r < rows; ++r) dot += g_i[r] * g_j[r];
            gram[i * m + j] = dot;
            gram[j * m + i] = dot;
        }
    }

    double max_row_sum = 0.0;
    for (int i = 0; i < m; ++i) {
        double row_sum = 0.0;
        for (int j = 0; j < m; ++j) row_sum += fabs(gram[i * m + j]);
        if (row_sum > max_row_sum) max_row_sum = row_sum;
    }

    double pg_step = 1.0 / ((max_row_sum / bs->prox_mu) + 1e-9);
    for (int i = 0; i < m; ++i) lambda[i] = 1.0 / (double)m;

    for (int iter = 0; iter < 80; ++iter) {
        for (int i = 0; i < m; ++i) {
            double gram_lambda = 0.0;
            for (int j = 0; j < m; ++j) gram_lambda += gram[i * m + j] * lambda[j];
            gradient[i] = linear[i] - gram_lambda / bs->prox_mu;
        }

        for (int i = 0; i < m; ++i) lambda[i] += pg_step * gradient[i];
        bundle_project_simplex(lambda, m);
    }

    for (int i = 0; i < m; ++i) {
        int slot = (start + i) % bs->memory;
        const double *g_i = bs->cut_slacks + (size_t)slot * (size_t)rows;
        for (int r = 0; r < rows; ++r) gbar[r] += lambda[i] * g_i[r];
    }

    for (int r = 0; r < rows; ++r) {
        double next = base_t[r] + gbar[r] / bs->prox_mu;
        t[r] = next > 0.0 ? next : 0.0;
    }

    double model_f = -DBL_MAX;
    for (int i = 0; i < m; ++i) {
        int slot = (start + i) % bs->memory;
        const double *g_i = bs->cut_slacks + (size_t)slot * (size_t)rows;
        double cut = bs->cut_alpha[slot];
        for (int r = 0; r < rows; ++r) cut -= g_i[r] * t[r];
        if (cut > model_f) model_f = cut;
    }

    double predicted_gain = -bs->center_value - model_f;
    if (predicted_gain <= bs->serious_tol) {
        free(base_t);
        free(linear);
        free(gram);
        free(lambda);
        free(gradient);
        free(gbar);
        free(slack);
        return 0;
    }

    bs->predicted_gain = predicted_gain;
    bs->has_trial = 1;

    free(base_t);
    free(linear);
    free(gram);
    free(lambda);
    free(gradient);
    free(gbar);
    free(slack);
    return 1;
}

/*
Row-targeted heuristic: for each uncovered row, pick a covering column with minimal Lagrangian score.
Tie order:
(1) more negative Lagrangian score
(2) higher weight
(3) covers more rows
(4) lower index.
*/
static void heuristic_row_min_rc(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights, /* tie-break: higher is better */
    int *sol, int *sol_len
) {
    (void)colsCoveringCount;

    bool *rowCovered = (bool*)calloc((size_t)rows, sizeof(bool));
    unsigned char *x = (unsigned char*)calloc((size_t)cols, 1);
    int out = 0, covered = 0;

    if (!rowCovered || !x) {
        free(rowCovered);
        free(x);
        *sol_len = -1;
        return;
    }

    for (int i = 0; i < rows; ++i) {
        if (rowCovered[i]) continue;

        /* choose candidate among columns covering row i */
        int best = -1;
        double bestLS = DBL_MAX;
        int bestOnes = -1;
        double bestWeight = -DBL_MAX;

        for (int k = 0; k < colsCoveringCount[i]; ++k) {
            int c = colsCovering[i][k];
            double ls = lagr_score ? lagr_score[c] : 0.0;
            double w  = weights ? weights[c] : 0.0;
            int ones = rowsCoveredCount[c];

            bool better = false;
            if (ls < bestLS) {
                better = true;
            } else if (fabs(ls - bestLS) <= EPS) {
                if (w > bestWeight) {
                    better = true;
                } else if (fabs(w - bestWeight) <= EPS) {
                    if (ones > bestOnes) {
                        better = true;
                    } else if (ones == bestOnes) {
                        if (best == -1 || c < best) {
                            better = true;
                        }
                    }
                }
            }
            if (better) {
                best = c;
                bestLS = ls;
                bestOnes = ones;
                bestWeight = w;
            }
        }

        if (best == -1) { /* should not happen if adjacency is valid */
            free(rowCovered);
            free(x);
            *sol_len = -1;
            return;
        }

        if (!x[best]) {
            x[best] = 1;
            sol[out++] = best;
            for (int c = 0; c < rowsCoveredCount[best]; ++c) {
                int r = rowsCovered[best][c];
                if (!rowCovered[r]) {
                    rowCovered[r] = true;
                    covered++;
                }
            }
        }

        if (covered >= rows) break;
    }

    *sol_len = out;

    /* Removal: try to drop columns in order of worst Lagrangian score first */
    if (out > 0) {
        /* Build cover counts */
        int *coverCount = (int*)calloc((size_t)rows, sizeof(int));
        int *order = (int*)malloc((size_t)out * sizeof(int));
        double *ls_sel = (double*)malloc((size_t)out * sizeof(double));

        if (coverCount && order && ls_sel) {
            for (int i = 0; i < out; ++i) {
                int s = sol[i];
                order[i] = i; /* index into sol */
                ls_sel[i] = lagr_score ? lagr_score[s] : 0.0;
                for (int c = 0; c < rowsCoveredCount[s]; ++c) {
                    int r = rowsCovered[s][c];
                    coverCount[r]++;
                }
            }
            /* sort order by score descending (worst first) */
            for (int i = 0; i < out; ++i) {
                for (int j = i + 1; j < out; ++j) {
                    if (ls_sel[order[j]] > ls_sel[order[i]]) {
                        int tmp = order[i];
                        order[i] = order[j];
                        order[j] = tmp;
                    }
                }
            }
            /* attempt removals */
            for (int o = 0; o < out; ++o) {
                int idx = order[o];
                if (idx < 0) continue;
                int s = sol[idx];
                bool can_remove = true;

                for (int c = 0; c < rowsCoveredCount[s]; ++c) {
                    int r = rowsCovered[s][c];
                    if (coverCount[r] <= 1) {
                        can_remove = false;
                        break;
                    }
                }

                if (can_remove) {
                    for (int c = 0; c < rowsCoveredCount[s]; ++c) {
                        int r = rowsCovered[s][c];
                        coverCount[r]--;
                    }

                    /* remove from sol: shift left */
                    for (int j = idx + 1; j < out; ++j) {
                        sol[j - 1] = sol[j];
                    }
                    out--;

                    /* adjust later order indices */
                    for (int j = o + 1; j < *sol_len; ++j) {
                        if (order[j] > idx) order[j]--;
                    }
                }
            }
            *sol_len = out;
        }

        free(coverCount);
        free(order);
        free(ls_sel);
    }

    prune_redundancy(rows, rowsCovered, rowsCoveredCount, sol, sol_len);

    free(rowCovered);
    free(x);
}

static double solution_total_weight(const int *sol, int sol_len, const double *weights) {
    if (!sol || sol_len <= 0 || !weights) return 0.0;
    double s = 0.0;
    for (int i = 0; i < sol_len; ++i) s += weights[sol[i]];
    return s;
}

static void lagr_initialize_multipliers(
    int rows,
    const int *colsCoveringCount,
    int rarity_init,
    double *t
) {
    if (!t || rows <= 0) return;

    if (!rarity_init || !colsCoveringCount) {
        for (int i = 0; i < rows; ++i) {
            t[i] = 1.0;
        }
        return;
    }

    double mean_support = 0.0;
    for (int i = 0; i < rows; ++i) {
        mean_support += (double)colsCoveringCount[i];
    }
    mean_support /= (double)rows;
    if (mean_support <= 0.0) mean_support = 1.0;

    for (int i = 0; i < rows; ++i) {
        double support = (double)colsCoveringCount[i];
        if (support <= 0.0) support = 1.0;
        t[i] = mean_support / support;
        if (t[i] < 0.1) t[i] = 0.1;
    }
}

/* Proxy for acceptance: if weights exist, use their sum; else use coverage volume as a cheap heuristic */
static double solution_proxy_value(
    const int *sol,
    int sol_len,
    const double *weights,
    const int *rowsCoveredCount
) {
    if (!sol || sol_len <= 0) return 0.0;
    double v = 0.0;
    if (weights) {
        for (int i = 0; i < sol_len; ++i) v += weights[sol[i]];
    } else if (rowsCoveredCount) {
        for (int i = 0; i < sol_len; ++i) v += (double)rowsCoveredCount[sol[i]];
    }
    return v;
}

/* Local search: try drop-1 + greedy repair + prune; accept only if shorter */
static void local_search_drop1_repair(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int *sol,
    int *sol_len,
    int max_passes
) {
    if (!sol || !sol_len || *sol_len <= 1) return;
    if (max_passes <= 0) return;

    const double PROXY_EPS = 1e-12;
    double curr_proxy = solution_proxy_value(sol, *sol_len, weights, rowsCoveredCount);

    bool *rowCovered = (bool*)malloc((size_t)rows * sizeof(bool));
    unsigned char *colSelected = (unsigned char*)calloc((size_t)cols, 1);
    int *cand = (int*)malloc((size_t)cols * sizeof(int));
    if (!rowCovered || !colSelected || !cand) {
        free(rowCovered);
        free(colSelected);
        free(cand);
        return;
    }

    for (int pass = 0; pass < max_passes; ++pass) {
        int improved = 0;

        for (int drop_idx = 0; drop_idx < *sol_len; ++drop_idx) {
            int cand_len = 0;
            memset(colSelected, 0, (size_t)cols);

            for (int i = 0; i < *sol_len; ++i) {
                if (i == drop_idx) continue;
                int c = sol[i];
                cand[cand_len++] = c;
                colSelected[c] = 1;
            }

            memset(rowCovered, 0, (size_t)rows * sizeof(bool));
            int covered = 0;
            for (int i = 0; i < cand_len; ++i) {
                int c = cand[i];
                for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                    int r = rowsCovered[c][k];
                    if (!rowCovered[r]) { rowCovered[r] = true; covered++; }
                }
            }

            while (covered < rows) {
                int r0 = -1;
                for (int r = 0; r < rows; ++r) {
                    if (!rowCovered[r]) { r0 = r; break; }
                }
                if (r0 < 0) break;

                int best = -1, bestNew = -1;
                double bestLS = DBL_MAX, bestW = -DBL_MAX;

                for (int kk = 0; kk < colsCoveringCount[r0]; ++kk) {
                    int c = colsCovering[r0][kk];
                    if (colSelected[c]) continue;

                    int newCover = 0;
                    for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                        int rr = rowsCovered[c][k];
                        if (!rowCovered[rr]) newCover++;
                    }
                    if (newCover <= 0) continue;

                    double ls = lagr_score ? lagr_score[c] : 0.0;
                    double w  = weights ? weights[c] : 0.0;

                    bool better = false;
                    if (newCover > bestNew) better = true;
                    else if (newCover == bestNew) {
                        if (ls < bestLS) better = true;
                        else if (fabs(ls - bestLS) <= EPS) {
                            if (w > bestW) better = true;
                            else if (fabs(w - bestW) <= EPS) {
                                if (best == -1 || c < best) better = true;
                            }
                        }
                    }

                    if (better) {
                        best = c; bestNew = newCover; bestLS = ls; bestW = w;
                    }
                }

                if (best < 0) { cand_len = -1; break; }

                colSelected[best] = 1;
                cand[cand_len++] = best;
                for (int k = 0; k < rowsCoveredCount[best]; ++k) {
                    int rr = rowsCovered[best][k];
                    if (!rowCovered[rr]) { rowCovered[rr] = true; covered++; }
                }
            }

            if (cand_len <= 0 || cand_len > *sol_len) continue;

            prune_redundancy(rows, rowsCovered, rowsCoveredCount, cand, &cand_len);

            if (cand_len > 0 && cand_len <= *sol_len) {
                double cand_proxy = solution_proxy_value(cand, cand_len, weights, rowsCoveredCount);
                /* Accept if proxy strictly improves, even for equal length (swap), or for shorter */
                if (cand_proxy > curr_proxy + PROXY_EPS) {
                    memcpy(sol, cand, (size_t)cand_len * sizeof(int));
                    *sol_len = cand_len;
                    curr_proxy = cand_proxy;
                    improved = 1;
                    break; /* restart after improvement */
                }
            }
        }

        if (!improved) break;
    }

    free(rowCovered);
    free(colSelected);
    free(cand);
}

/* Local search: drop-2 + greedy repair + prune; accept only if shorter */
static void local_search_drop2_repair(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int *sol,
    int *sol_len,
    int max_passes
) {
    if (!sol || !sol_len || *sol_len <= 2) return;
    if (max_passes <= 0) return;

    bool *rowCovered = (bool*)malloc((size_t)rows * sizeof(bool));
    unsigned char *colSelected = (unsigned char*)calloc((size_t)cols, 1);
    int *cand = (int*)malloc((size_t)cols * sizeof(int));
    if (!rowCovered || !colSelected || !cand) {
        free(rowCovered);
        free(colSelected);
        free(cand);
        return;
    }

    for (int pass = 0; pass < max_passes; ++pass) {
        int improved = 0;

        for (int drop_i = 0; drop_i < *sol_len; ++drop_i) {
            for (int drop_j = drop_i + 1; drop_j < *sol_len; ++drop_j) {
                int cand_len = 0;
                memset(colSelected, 0, (size_t)cols);

                for (int i = 0; i < *sol_len; ++i) {
                    if (i == drop_i || i == drop_j) continue;
                    int c = sol[i];
                    cand[cand_len++] = c;
                    colSelected[c] = 1;
                }

                memset(rowCovered, 0, (size_t)rows * sizeof(bool));
                int covered = 0;
                for (int i = 0; i < cand_len; ++i) {
                    int c = cand[i];
                    for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                        int r = rowsCovered[c][k];
                        if (!rowCovered[r]) { rowCovered[r] = true; covered++; }
                    }
                }

                while (covered < rows) {
                    int r0 = -1;
                    for (int r = 0; r < rows; ++r) {
                        if (!rowCovered[r]) { r0 = r; break; }
                    }
                    if (r0 < 0) break;

                    int best = -1, bestNew = -1;
                    double bestLS = DBL_MAX, bestW = -DBL_MAX;

                    for (int kk = 0; kk < colsCoveringCount[r0]; ++kk) {
                        int c = colsCovering[r0][kk];
                        if (colSelected[c]) continue;

                        int newCover = 0;
                        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                            int rr = rowsCovered[c][k];
                            if (!rowCovered[rr]) newCover++;
                        }
                        if (newCover <= 0) continue;

                        double ls = lagr_score ? lagr_score[c] : 0.0;
                        double w  = weights ? weights[c] : 0.0;

                        bool better = false;
                        if (newCover > bestNew) better = true;
                        else if (newCover == bestNew) {
                            if (ls < bestLS) better = true;
                            else if (fabs(ls - bestLS) <= EPS) {
                                if (w > bestW) better = true;
                                else if (fabs(w - bestW) <= EPS) {
                                    if (best == -1 || c < best) better = true;
                                }
                            }
                        }

                        if (better) {
                            best = c; bestNew = newCover; bestLS = ls; bestW = w;
                        }
                    }

                    if (best < 0) { cand_len = -1; break; }

                    colSelected[best] = 1;
                    cand[cand_len++] = best;
                    for (int k = 0; k < rowsCoveredCount[best]; ++k) {
                        int rr = rowsCovered[best][k];
                        if (!rowCovered[rr]) { rowCovered[rr] = true; covered++; }
                    }
                }

                if (cand_len <= 0 || cand_len >= *sol_len) continue;

                prune_redundancy(rows, rowsCovered, rowsCoveredCount, cand, &cand_len);

                if (cand_len > 0 && cand_len < *sol_len) {
                    memcpy(sol, cand, (size_t)cand_len * sizeof(int));
                    *sol_len = cand_len;
                    improved = 1;
                    break; /* restart after improvement */
                }
            }

            if (improved) break;
        }

        if (!improved) break;
    }

    free(rowCovered);
    free(colSelected);
    free(cand);
}

typedef struct {
    int col;
    int new_cover;
    int cover_count;
    double lagr_score;
    double weight;
} PolishCandidate;

typedef struct {
    int rows;
    int cols;
    int **rowsCovered;
    int *rowsCoveredCount;
    int **colsCovering;
    int *colsCoveringCount;
    const double *lagr_score;
    const double *weights;
    const unsigned char *allowed; /* NULL = every column usable */
    const int *forced;            /* columns pre-selected before the search */
    int forced_len;
    int target;                   /* free slots available to the search (excludes forced) */
    int out_stride;               /* forced_len + target */
    long node_limit;
    long nodes;
    int aborted;                  /* node budget exhausted */
    int want;                     /* number of distinct covers to collect */
    int found;                    /* covers collected so far */
    int *coverCount;
    unsigned char *colSelected;
    int *chosen;
    int *out;                     /* want * out_stride ints */
    int *out_lens;                /* want ints */
    int *scratch;                 /* 2 * out_stride ints, for dedup */
    PolishCandidate *cand_buf;    /* target * max_row_degree entries */
    int max_row_degree;
    int max_cover_static;         /* max rows covered by any usable column */
    /* bitset fast path (rows <= 64) */
    const uint64_t *col_mask;     /* per-column row coverage mask, or NULL */
    uint64_t full_mask;           /* all rows covered */
} PolishSearch;

static int lagr_popcount64(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int n = 0;
    while (x) { x &= x - 1; ++n; }
    return n;
#endif
}

static bool polish_candidate_better(const PolishCandidate *a, const PolishCandidate *b) {
    if (a->new_cover != b->new_cover) return a->new_cover > b->new_cover;
    if (fabs(a->lagr_score - b->lagr_score) > EPS) return a->lagr_score < b->lagr_score;
    if (a->cover_count != b->cover_count) return a->cover_count > b->cover_count;
    if (fabs(a->weight - b->weight) > EPS) return a->weight > b->weight;
    return a->col < b->col;
}

static void polish_sort_candidates(PolishCandidate *candidates, int count) {
    if (count > 32) {
        /*
        Wide nodes: a full sort is wasted work — on failure every branch is
        explored regardless of order, and a successful dive follows the front
        candidates. Selection-sort just the first few positions in O(8n);
        the tail keeps its deterministic discovery order.
        */
        int front = count < 8 ? count : 8;
        for (int i = 0; i < front; ++i) {
            int best = i;
            for (int j = i + 1; j < count; ++j) {
                if (polish_candidate_better(&candidates[j], &candidates[best])) best = j;
            }
            if (best != i) {
                PolishCandidate tmp = candidates[i];
                candidates[i] = candidates[best];
                candidates[best] = tmp;
            }
        }
        return;
    }
    for (int i = 1; i < count; ++i) {
        PolishCandidate value = candidates[i];
        int j = i - 1;
        while (j >= 0 && polish_candidate_better(&value, &candidates[j])) {
            candidates[j + 1] = candidates[j];
            --j;
        }
        candidates[j + 1] = value;
    }
}

static int polish_new_cover(
    const PolishSearch *search,
    int col
) {
    int new_cover = 0;
    for (int i = 0; i < search->rowsCoveredCount[col]; ++i) {
        int row = search->rowsCovered[col][i];
        if (search->coverCount[row] == 0) ++new_cover;
    }
    return new_cover;
}

static int polish_choose_row(
    const PolishSearch *search,
    int uncovered,
    int slots_left
) {
    /* bound: every column adds at most max_cover_static newly covered rows */
    if (uncovered > slots_left * search->max_cover_static) return -1;

    int best_row = -1;
    int best_count = INT_MAX;

    for (int row = 0; row < search->rows; ++row) {
        if (search->coverCount[row] > 0) continue;

        /* any unselected usable column covering this (uncovered) row
           necessarily has new_cover >= 1, so a plain count suffices */
        int viable = 0;
        for (int i = 0; i < search->colsCoveringCount[row]; ++i) {
            int col = search->colsCovering[row][i];
            if (search->colSelected[col]) continue;
            if (search->allowed && !search->allowed[col]) continue;
            ++viable;
            if (viable >= best_count) break; /* cannot become the minimum */
        }

        if (viable == 0) return -1;
        if (viable < best_count) {
            best_count = viable;
            best_row = row;
        }
    }

    return best_row;
}

static int pool_cmp_int_asc(const void *a, const void *b);
static void pool_clear(int **pool_solutions, int *pool_count);
static int pool_add_if_new(
    int **pool_solutions,
    int *pool_count,
    int max_pool,
    const int *sol_sorted,
    int len
);

static bool polish_search_rec(
    PolishSearch *search,
    int depth,
    int covered
) {
    if (covered >= search->rows) {
        int len = search->forced_len + depth;
        int *norm = search->scratch;

        for (int i = 0; i < search->forced_len; ++i) norm[i] = search->forced[i];
        memcpy(norm + search->forced_len, search->chosen, (size_t)depth * sizeof(int));
        qsort(norm, (size_t)len, sizeof(int), pool_cmp_int_asc);

        for (int s = 0; s < search->found; ++s) {
            if (
                search->out_lens[s] == len &&
                memcmp(search->out + (size_t)s * search->out_stride, norm, (size_t)len * sizeof(int)) == 0
            ) {
                return false; /* duplicate cover: keep searching */
            }
        }

        memcpy(search->out + (size_t)search->found * search->out_stride, norm, (size_t)len * sizeof(int));
        search->out_lens[search->found] = len;
        search->found++;
        return search->found >= search->want;
    }

    if (depth >= search->target) return false;
    if (search->nodes++ >= search->node_limit) {
        search->aborted = 1;
        return true; /* unwind immediately */
    }

    int slots_left = search->target - depth;
    int uncovered = search->rows - covered;
    int row = polish_choose_row(search, uncovered, slots_left);
    if (row < 0) return false;

    int raw_count = search->colsCoveringCount[row];
    PolishCandidate *candidates = search->cand_buf + (size_t)depth * (size_t)search->max_row_degree;

    int candidate_count = 0;
    for (int i = 0; i < raw_count; ++i) {
        int col = search->colsCovering[row][i];
        if (search->colSelected[col]) continue;
        if (search->allowed && !search->allowed[col]) continue;

        int new_cover = polish_new_cover(search, col);
        if (new_cover <= 0) continue;

        candidates[candidate_count++] = (PolishCandidate){
            .col = col,
            .new_cover = new_cover,
            .cover_count = search->rowsCoveredCount[col],
            .lagr_score = search->lagr_score ? search->lagr_score[col] : 0.0,
            .weight = search->weights ? search->weights[col] : 0.0
        };
    }

    polish_sort_candidates(candidates, candidate_count);

    for (int i = 0; i < candidate_count; ++i) {
        int col = candidates[i].col;
        int added = 0;

        search->colSelected[col] = 1;
        search->chosen[depth] = col;

        for (int j = 0; j < search->rowsCoveredCount[col]; ++j) {
            int covered_row = search->rowsCovered[col][j];
            if (search->coverCount[covered_row] == 0) ++added;
            ++search->coverCount[covered_row];
        }

        bool stop = polish_search_rec(search, depth + 1, covered + added);

        for (int j = 0; j < search->rowsCoveredCount[col]; ++j) {
            int covered_row = search->rowsCovered[col][j];
            --search->coverCount[covered_row];
        }

        search->colSelected[col] = 0;

        if (stop) return true;
    }

    return false;
}

/* Bitset variant of polish_search_rec for rows <= 64: coverage is a single
   uint64_t, new-cover checks are popcounts, undo restores the saved mask. */
static bool polish_search_rec_bits(
    PolishSearch *search,
    int depth,
    uint64_t covered_mask
) {
    if (covered_mask == search->full_mask) {
        int len = search->forced_len + depth;
        int *norm = search->scratch;

        for (int i = 0; i < search->forced_len; ++i) norm[i] = search->forced[i];
        memcpy(norm + search->forced_len, search->chosen, (size_t)depth * sizeof(int));
        qsort(norm, (size_t)len, sizeof(int), pool_cmp_int_asc);

        for (int s = 0; s < search->found; ++s) {
            if (
                search->out_lens[s] == len &&
                memcmp(search->out + (size_t)s * search->out_stride, norm, (size_t)len * sizeof(int)) == 0
            ) {
                return false;
            }
        }

        memcpy(search->out + (size_t)search->found * search->out_stride, norm, (size_t)len * sizeof(int));
        search->out_lens[search->found] = len;
        search->found++;
        return search->found >= search->want;
    }

    if (depth >= search->target) return false;
    if (search->nodes++ >= search->node_limit) {
        search->aborted = 1;
        return true;
    }

    uint64_t uncovered_mask = search->full_mask & ~covered_mask;
    int uncovered = lagr_popcount64(uncovered_mask);
    int slots_left = search->target - depth;

    if (uncovered > slots_left * search->max_cover_static) return false;

    /* MRV row choice: fewest unselected candidates */
    int best_row = -1;
    int best_count = INT_MAX;
    for (uint64_t m = uncovered_mask; m; m &= m - 1) {
        int row = __builtin_ctzll(m);
        int viable = 0;
        for (int i = 0; i < search->colsCoveringCount[row]; ++i) {
            int col = search->colsCovering[row][i];
            if (search->colSelected[col]) continue;
            ++viable;
            if (viable >= best_count) break;
        }
        if (viable == 0) return false;
        if (viable < best_count) {
            best_count = viable;
            best_row = row;
        }
    }
    if (best_row < 0) return false;

    PolishCandidate *candidates = search->cand_buf + (size_t)depth * (size_t)search->max_row_degree;
    int candidate_count = 0;
    for (int i = 0; i < search->colsCoveringCount[best_row]; ++i) {
        int col = search->colsCovering[best_row][i];
        if (search->colSelected[col]) continue;

        int new_cover = lagr_popcount64(search->col_mask[col] & uncovered_mask);
        if (new_cover <= 0) continue;

        candidates[candidate_count++] = (PolishCandidate){
            .col = col,
            .new_cover = new_cover,
            .cover_count = search->rowsCoveredCount[col],
            .lagr_score = search->lagr_score ? search->lagr_score[col] : 0.0,
            .weight = search->weights ? search->weights[col] : 0.0
        };
    }

    polish_sort_candidates(candidates, candidate_count);

    for (int i = 0; i < candidate_count; ++i) {
        int col = candidates[i].col;

        search->colSelected[col] = 1;
        search->chosen[depth] = col;

        bool stop = polish_search_rec_bits(
            search,
            depth + 1,
            covered_mask | search->col_mask[col]
        );

        search->colSelected[col] = 0;

        if (stop) return true;
    }

    return false;
}

/*
Bounded exact search for covers of size <= target_total.
- allowed: optional column mask (NULL = all active columns)
- forced/forced_len: columns fixed into every cover before the search starts
- want: how many distinct covers to collect (enumeration for pool mode)
- out: want * target_total ints; out_lens: length of each collected cover
Returns the number of covers found; *aborted_out is set when the node budget
was exhausted (in which case a "not found" answer is inconclusive).
*/
static int lagr_exact_cover_search(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    const unsigned char *allowed,
    const int *forced,
    int forced_len,
    int target_total,
    long node_limit,
    int want,
    int *out,
    int *out_lens,
    int *aborted_out
) {
    if (aborted_out) *aborted_out = 0;
    if (target_total <= 0 || node_limit <= 0 || want <= 0) return 0;
    if (forced_len > target_total) return 0;

    int free_slots = target_total - forced_len;

    int max_row_degree = 1;
    int max_cover_static = 1;
    for (int r = 0; r < rows; ++r) {
        if (colsCoveringCount[r] > max_row_degree) max_row_degree = colsCoveringCount[r];
        for (int i = 0; i < colsCoveringCount[r]; ++i) {
            int c = colsCovering[r][i];
            if (rowsCoveredCount[c] > max_cover_static) max_cover_static = rowsCoveredCount[c];
        }
    }

    int use_bits = (rows <= 64 && allowed == NULL);

    int *coverCount = (int*)calloc((size_t)rows, sizeof(int));
    unsigned char *colSelected = (unsigned char*)calloc((size_t)cols, 1);
    int *chosen = (int*)malloc((size_t)(free_slots > 0 ? free_slots : 1) * sizeof(int));
    int *scratch = (int*)malloc((size_t)target_total * 2 * sizeof(int));
    PolishCandidate *cand_buf = (PolishCandidate*)malloc(
        (size_t)(free_slots > 0 ? free_slots : 1) * (size_t)max_row_degree * sizeof(PolishCandidate)
    );
    uint64_t *col_mask = use_bits
        ? (uint64_t*)calloc((size_t)cols, sizeof(uint64_t))
        : NULL;

    if (!coverCount || !colSelected || !chosen || !scratch || !cand_buf || (use_bits && !col_mask)) {
        free(coverCount);
        free(colSelected);
        free(chosen);
        free(scratch);
        free(cand_buf);
        free(col_mask);
        return 0;
    }

    if (use_bits) {
        /* per-column coverage masks, filled once per core column */
        for (int r = 0; r < rows; ++r) {
            for (int i = 0; i < colsCoveringCount[r]; ++i) {
                int c = colsCovering[r][i];
                if (col_mask[c] == 0) {
                    uint64_t m = 0;
                    for (int j = 0; j < rowsCoveredCount[c]; ++j) {
                        m |= 1ULL << rowsCovered[c][j];
                    }
                    col_mask[c] = m;
                }
            }
        }
    }

    int covered = 0;
    for (int i = 0; i < forced_len; ++i) {
        int c = forced[i];
        colSelected[c] = 1;
        for (int j = 0; j < rowsCoveredCount[c]; ++j) {
            int r = rowsCovered[c][j];
            if (coverCount[r] == 0) ++covered;
            ++coverCount[r];
        }
    }

    PolishSearch search = {
        .rows = rows,
        .cols = cols,
        .rowsCovered = rowsCovered,
        .rowsCoveredCount = rowsCoveredCount,
        .colsCovering = colsCovering,
        .colsCoveringCount = colsCoveringCount,
        .lagr_score = lagr_score,
        .weights = weights,
        .allowed = allowed,
        .forced = forced,
        .forced_len = forced_len,
        .target = free_slots,
        .out_stride = target_total,
        .node_limit = node_limit,
        .nodes = 0,
        .aborted = 0,
        .want = want,
        .found = 0,
        .coverCount = coverCount,
        .colSelected = colSelected,
        .chosen = chosen,
        .out = out,
        .out_lens = out_lens,
        .scratch = scratch,
        .cand_buf = cand_buf,
        .max_row_degree = max_row_degree,
        .max_cover_static = max_cover_static,
        .col_mask = col_mask,
        .full_mask = (rows >= 64) ? ~0ULL : ((1ULL << rows) - 1ULL)
    };

    if (use_bits) {
        uint64_t covered_mask = 0;
        for (int i = 0; i < forced_len; ++i) {
            covered_mask |= col_mask[forced[i]];
        }
        polish_search_rec_bits(&search, 0, covered_mask);
    } else {
        polish_search_rec(&search, 0, covered);
    }

    if (aborted_out) *aborted_out = search.aborted;

    free(coverCount);
    free(colSelected);
    free(chosen);
    free(scratch);
    free(cand_buf);
    free(col_mask);

    return search.found;
}

/*
Lagrangian exact finish:
1. Recompute reduced costs at the best multipliers found (caller passes ls_best/zlb_best).
2. Reduced-cost variable fixing for a target size T:
   - exclude column c when zlb + max(0, rc_c) > T (cannot be in any cover of cost <= T)
   - force column c when rc_c < 0 and zlb - rc_c > T (must be in every such cover)
3. Bounded exact search at T = UB-1; iterate downward while covers are found.
   An exhaustive (non-aborted) failure at T proves UB optimal, which also lifts
   the integer lower bound to UB.
4. Optional pool enumeration at the final UB.
*/
/*
Build a core-restricted copy of the row -> columns adjacency containing only
columns that can appear in a cover of cost <= target (reduced-cost fixing at
the best multipliers). Also collects columns that MUST be in every such cover.
Returns the core size, or -1 on OOM. *core_covering / *core_counts are
allocated as one block each; caller frees both.
*/
static int lagr_build_core(
    int rows,
    int cols,
    const int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *ls_best,
    double zlb_best,
    int target,
    unsigned char *allowed,      /* cols scratch, filled */
    int *forced,                 /* cols scratch, filled */
    int *forced_len_out,
    int ***core_covering_out,
    int **core_counts_out
) {
    const double FIXTOL = 1e-7;

    int forced_len = 0;
    int core_count = 0;
    for (int c = 0; c < cols; ++c) {
        double rc_c = ls_best[c];
        if (rowsCoveredCount[c] <= 0 || zlb_best + (rc_c > 0.0 ? rc_c : 0.0) > (double)target + FIXTOL) {
            allowed[c] = 0;
        } else {
            allowed[c] = 1;
            ++core_count;
            if (rc_c < 0.0 && zlb_best - rc_c > (double)target + FIXTOL) {
                forced[forced_len++] = c;
            }
        }
    }
    *forced_len_out = forced_len;

    if (!core_covering_out || !core_counts_out) return core_count;

    int **core_covering = (int**)malloc((size_t)rows * sizeof(int*));
    int *core_counts = (int*)malloc((size_t)rows * sizeof(int));
    if (!core_covering || !core_counts) {
        free(core_covering);
        free(core_counts);
        return -1;
    }

    size_t nnz = 0;
    for (int r = 0; r < rows; ++r) {
        int cnt = 0;
        for (int k = 0; k < colsCoveringCount[r]; ++k) {
            if (allowed[colsCovering[r][k]]) ++cnt;
        }
        core_counts[r] = cnt;
        nnz += (size_t)cnt;
    }

    int *pool = (int*)malloc((nnz > 0 ? nnz : 1) * sizeof(int));
    if (!pool) {
        free(core_covering);
        free(core_counts);
        return -1;
    }

    size_t off = 0;
    for (int r = 0; r < rows; ++r) {
        core_covering[r] = pool + off;
        for (int k = 0; k < colsCoveringCount[r]; ++k) {
            int c = colsCovering[r][k];
            if (allowed[c]) pool[off++] = c;
        }
    }

    *core_covering_out = core_covering;   /* core_covering[0] owns the pool */
    *core_counts_out = core_counts;
    return core_count;
}

static void lagr_free_core(int **core_covering, int *core_counts, int rows) {
    if (core_covering) {
        if (rows > 0 && core_covering[0]) free(core_covering[0]); /* the pool */
        free(core_covering);
    }
    free(core_counts);
}

/*
Chvatal-Gomory k-cover bound over a core-restricted chart (all rows
uncovered): a set S of rows such that every core column covers at most k
rows of S forces any cover to use at least ceil(|S| / k) columns. Greedy
separation, rows with the fewest core candidates first. Used as a root
test before launching a bounded search: a bound above the target proves
the target unreachable without expanding a single node.
*/
static int lagr_kcover_bound(
    int rows,
    int cols,
    int **core_covering,
    int *core_counts,
    int k,
    int *hits,       /* cols scratch */
    int *rows_order  /* rows scratch */
) {
    if (k <= 1 || rows <= 0 || rows > 256) return 0;

    for (int r = 0; r < rows; ++r) rows_order[r] = r;

    /* ascending number of core candidates, ties by row id */
    for (int i = 1; i < rows; ++i) {
        int r = rows_order[i];
        int j = i - 1;
        while (j >= 0 && (core_counts[rows_order[j]] > core_counts[r] ||
               (core_counts[rows_order[j]] == core_counts[r] && rows_order[j] > r))) {
            rows_order[j + 1] = rows_order[j];
            --j;
        }
        rows_order[j + 1] = r;
    }

    memset(hits, 0, (size_t)cols * sizeof(int));

    int in_s = 0;
    for (int idx = 0; idx < rows; ++idx) {
        int r = rows_order[idx];
        int fits = 1;

        for (int i = 0; i < core_counts[r]; ++i) {
            if (hits[core_covering[r][i]] >= k) {
                fits = 0;
                break;
            }
        }
        if (!fits) continue;

        for (int i = 0; i < core_counts[r]; ++i) {
            hits[core_covering[r][i]]++;
        }
        ++in_s;
    }

    return (in_s + k - 1) / k;
}

static void lagr_exact_finish(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *weights,
    const double *ls_best,
    double zlb_best,
    long node_limit,
    int core_cap,
    int *solution,
    int *solmin,
    double *bestLB_io,
    int *proved_optimal_io,
    int max_pool,
    int **pool_solutions,
    int *pool_count,
    int *pool_min_len,
    int *norm_buf
) {
    if (!solution || !solmin || *solmin <= 0) return;

    unsigned char *allowed = (unsigned char*)malloc((size_t)cols);
    int *forced = (int*)malloc((size_t)cols * sizeof(int));
    int *one = (int*)malloc((size_t)(*solmin) * sizeof(int));
    int *found_buf = NULL;
    int *found_lens = NULL;
    if (!allowed || !forced || !one) {
        free(allowed);
        free(forced);
        free(one);
        return;
    }

    int UB = *solmin;
    int proved = proved_optimal_io ? *proved_optimal_io : 0;

    /*
    Improvement loop. Searching size <= target only makes sense while
    target >= LB. When the core at a target is too large (weak fixing),
    step down to a tighter target instead of burning the node budget:
    the mask sharpens as the target approaches ZLB.
    An exhaustive failure at target T proves no cover of size <= T exists,
    lifting the integer LB to T+1; when T == UB-1 that proves UB optimal.
    */
    while (!proved && UB > 1) {
        int target = UB - 1;
        int improved = 0;

        while (target >= 1 && (!bestLB_io || (double)target >= *bestLB_io - EPS)) {
            int forced_len = 0;
            int **core_covering = NULL;
            int *core_counts = NULL;

            int core_count = lagr_build_core(
                rows, cols,
                rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls_best, zlb_best,
                target,
                allowed, forced, &forced_len,
                NULL, NULL
            );

            if (forced_len > target) {
                /* more must-have columns than slots: size <= target impossible */
                if (bestLB_io && (double)(target + 1) > *bestLB_io) *bestLB_io = (double)(target + 1);
                if (target == UB - 1) proved = 1;
                break;
            }

            if (core_count > core_cap) {
                --target; /* tighter target, smaller core */
                continue;
            }

            core_count = lagr_build_core(
                rows, cols,
                rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls_best, zlb_best,
                target,
                allowed, forced, &forced_len,
                &core_covering, &core_counts
            );
            if (core_count < 0) break; /* OOM: skip quietly */

            /* root k-cover cut: may prove target unreachable for free */
            {
                int *hits = (int*)malloc((size_t)cols * sizeof(int));
                int *rows_order = (int*)malloc((size_t)rows * sizeof(int));
                int cut_lb = 0;
                if (hits && rows_order) {
                    for (int k = 2; k <= 3; ++k) {
                        int b = lagr_kcover_bound(
                            rows, cols, core_covering, core_counts, k, hits, rows_order
                        );
                        if (b > cut_lb) cut_lb = b;
                    }
                }
                free(hits);
                free(rows_order);

                if (cut_lb > target) {
                    lagr_free_core(core_covering, core_counts, rows);
                    if (bestLB_io && (double)(target + 1) > *bestLB_io) {
                        *bestLB_io = (double)(target + 1);
                    }
                    if (target == UB - 1) proved = 1;
                    break;
                }
            }

            int found_len_one = 0;
            int aborted = 0;
            int n = lagr_exact_cover_search(
                rows, cols,
                rowsCovered, rowsCoveredCount,
                core_covering, core_counts,
                ls_best, weights,
                NULL, /* adjacency already core-restricted */
                forced, forced_len,
                target,
                node_limit,
                1,
                one, &found_len_one,
                &aborted
            );

            lagr_free_core(core_covering, core_counts, rows);

            if (n > 0) {
                memcpy(solution, one, (size_t)found_len_one * sizeof(int));
                *solmin = found_len_one;
                UB = found_len_one;
                improved = 1;
                break;
            }

            if (aborted) break; /* inconclusive: give up on this and lower targets */

            /* exhaustive failure: no cover of size <= target exists */
            if (bestLB_io && (double)(target + 1) > *bestLB_io) *bestLB_io = (double)(target + 1);
            if (target == UB - 1) proved = 1;
            break;
        }

        if (!improved) break;
    }

    if (proved) {
        if (bestLB_io && (double)UB > *bestLB_io) *bestLB_io = (double)UB;
        if (proved_optimal_io) *proved_optimal_io = 1;
    }

    /* ---- pool enumeration at the final size ---- */
    if (max_pool > 1 && pool_solutions && pool_count && pool_min_len && norm_buf) {
        if (UB < *pool_min_len) {
            /* the improvement loop beat every pooled solution: resync the pool */
            pool_clear(pool_solutions, pool_count);
            *pool_min_len = UB;
            memcpy(norm_buf, solution, (size_t)UB * sizeof(int));
            qsort(norm_buf, (size_t)UB, sizeof(int), pool_cmp_int_asc);
            pool_add_if_new(pool_solutions, pool_count, max_pool, norm_buf, UB);
        }

        int restart = 1;
        int rounds = 0;

        while (restart && rounds < 4) {
            restart = 0;
            ++rounds;

            int target = UB;
            int forced_len = 0;
            int **core_covering = NULL;
            int *core_counts = NULL;

            int core_count = lagr_build_core(
                rows, cols,
                rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls_best, zlb_best,
                target,
                allowed, forced, &forced_len,
                NULL, NULL
            );
            if (forced_len > target) break;
            if (core_count > core_cap) break; /* enumeration too expensive here */

            core_count = lagr_build_core(
                rows, cols,
                rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls_best, zlb_best,
                target,
                allowed, forced, &forced_len,
                &core_covering, &core_counts
            );
            if (core_count < 0) break;

            free(found_buf);
            free(found_lens);
            found_buf = (int*)malloc((size_t)max_pool * (size_t)target * sizeof(int));
            found_lens = (int*)malloc((size_t)max_pool * sizeof(int));
            if (!found_buf || !found_lens) {
                lagr_free_core(core_covering, core_counts, rows);
                break;
            }

            int aborted = 0;
            int n = lagr_exact_cover_search(
                rows, cols,
                rowsCovered, rowsCoveredCount,
                core_covering, core_counts,
                ls_best, weights,
                NULL,
                forced, forced_len,
                target,
                node_limit,
                max_pool,
                found_buf, found_lens,
                &aborted
            );

            lagr_free_core(core_covering, core_counts, rows);

            for (int s = 0; s < n; ++s) {
                int len = found_lens[s];
                const int *cover = found_buf + (size_t)s * (size_t)target;

                if (len < *pool_min_len || len < UB) {
                    /* better cover surfaced during enumeration */
                    memcpy(solution, cover, (size_t)len * sizeof(int));
                    *solmin = len;
                    UB = len;
                    pool_clear(pool_solutions, pool_count);
                    *pool_min_len = len;
                    memcpy(norm_buf, cover, (size_t)len * sizeof(int));
                    qsort(norm_buf, (size_t)len, sizeof(int), pool_cmp_int_asc);
                    pool_add_if_new(pool_solutions, pool_count, max_pool, norm_buf, len);
                    restart = 1;
                    break;
                }

                if (len == *pool_min_len) {
                    memcpy(norm_buf, cover, (size_t)len * sizeof(int));
                    qsort(norm_buf, (size_t)len, sizeof(int), pool_cmp_int_asc);
                    pool_add_if_new(pool_solutions, pool_count, max_pool, norm_buf, len);
                }
            }
        }
    }

    free(found_buf);
    free(found_lens);
    free(allowed);
    free(forced);
    free(one);
}


/* ---------- Solution pool helpers ---------- */
static int pool_cmp_int_asc(const void *a, const void *b) {
    int ia = *(const int*)a, ib = *(const int*)b;
    return (ia > ib) - (ia < ib);
}

static int pool_equal_sets(const int *a, int la, const int *b, int lb) {
    if (la != lb) return 0;
    for (int i = 0; i < la; ++i) {
        if (a[i] != b[i]) return 0;
    }
    return 1;
}

static void pool_clear(int **pool_solutions, int *pool_count) {
    if (!pool_solutions || !pool_count) return;
    for (int i = 0; i < *pool_count; ++i) {
        free(pool_solutions[i]);
        pool_solutions[i] = NULL;
    }
    *pool_count = 0;
}

static int pool_add_if_new(
    int **pool_solutions,
    int *pool_count,
    int max_pool,
    const int *sol_sorted,
    int len
) {
    if (!pool_solutions || !pool_count || max_pool <= 0) return 0;
    for (int i = 0; i < *pool_count; ++i) {
        if (pool_equal_sets(sol_sorted, len, pool_solutions[i], len)) return 0; /* duplicate */
    }

    if (*pool_count >= max_pool) return 0; /* full */

    int *dst = (int*)malloc((size_t)len * sizeof(int));
    if (!dst) return 0;
    memcpy(dst, sol_sorted, (size_t)len * sizeof(int));
    pool_solutions[(*pool_count)++] = dst;
    return 1;
}

typedef struct {
    int max_iter;
    int heur_every;
    double step_coef;
    double step_min;
    double step_contract;
    double stabilization_beta;
    int halve_period;
    int local_passes;
    int drop2_enabled;
    int small_gap_threshold;
    int small_gap_extra_passes;
    int rarity_init;
    double deflection_alpha;
    int polish_enabled;
    long polish_node_limit;
    int polish_core_cap;
    int portfolio_enabled;
    int portfolio_max_profiles;
    double portfolio_deflection_alpha;
    int hybrid_bundle_portfolio;
    int bundle_enabled;
    int bundle_interval;
    int bundle_memory;
    double bundle_prox_mu;
    double bundle_min_mu;
    double bundle_max_mu;
    double bundle_null_expand;
    double bundle_serious_shrink;
    double bundle_serious_fraction;
    double bundle_serious_tol;
    double bundle_momentum;
} LagrangianConfig;

#define LAGR_PORTFOLIO_PROFILE_COUNT 6

/*
Exact completion (-l3): decide optimality with the ported scp_solver.
For target = UB-1 downward, reduced-cost fixing at the best multipliers
restricts the chart to the columns that can appear in a cover of size
<= target (a cover of size <= target has cost <= target, so any column
with zlb + max(0, rc) > target is provably excluded). The exact solver
then either finds such a cover (improvement: continue with the smaller
UB) or exhausts the space, proving that UB is optimal. No node budgets
are involved, so the outcome is always a certificate, memory permitting.
Returns 1 when optimality was proven (bestLB_io lifted to *solmin),
0 when the completion could not run (OOM).
*/
static int lagr_exact_complete(
    const int *pichart,
    int rows,
    int cols,
    const int *rowsCoveredCount,
    const double *ls_best,
    double zlb_best,
    int *solution,
    int *solmin,
    double *bestLB_io
) {
    const double FIXTOL = 1e-7;

    if (!pichart || !solution || !solmin || *solmin <= 0 || rows <= 0 || cols <= 0) {
        return 0;
    }

    const int nwords_rows = (rows + 63) / 64;

    int *col_map = (int*)malloc((size_t)cols * sizeof(int));
    int *sol01 = NULL;
    if (!col_map) return 0;

    int UB = *solmin;
    int proved = 0;

    while (UB > 1) {
        int target = UB - 1;

        if (bestLB_io && (double)target < *bestLB_io - EPS) {
            proved = 1; /* below the proven lower bound: nothing can exist */
            break;
        }

        /* core columns for this target */
        int ncore = 0;
        for (int c = 0; c < cols; ++c) {
            double rc_c = ls_best ? ls_best[c] : 0.0;
            if (rowsCoveredCount && rowsCoveredCount[c] <= 0) continue;
            if (ls_best && zlb_best + (rc_c > 0.0 ? rc_c : 0.0) > (double)target + FIXTOL) continue;
            col_map[ncore++] = c;
        }

        if (ncore <= 0) {
            proved = 1; /* no usable columns: no cover of size <= target */
            break;
        }

        /* compact problem */
        size_t nnz = 0;
        unsigned long long *col_masks = (unsigned long long*)calloc(
            (size_t)ncore * (size_t)nwords_rows, sizeof(unsigned long long));
        int *row_counts = (int*)calloc((size_t)rows, sizeof(int));
        int *row_starts = (int*)calloc((size_t)rows + 1, sizeof(int));
        double *priority = ls_best ? (double*)malloc((size_t)ncore * sizeof(double)) : NULL;

        if (!col_masks || !row_counts || !row_starts || (ls_best && !priority)) {
            free(col_masks);
            free(row_counts);
            free(row_starts);
            free(priority);
            break;
        }

        for (int i = 0; i < ncore; ++i) {
            int c = col_map[i];
            if (priority) priority[i] = ls_best[c];
            const int *colv = pichart + (size_t)c * (size_t)rows;
            for (int r = 0; r < rows; ++r) {
                if (colv[r]) {
                    col_masks[(size_t)i * nwords_rows + (r >> 6)] |= 1ULL << (r & 63);
                    row_counts[r]++;
                    ++nnz;
                }
            }
        }

        row_starts[0] = 0;
        for (int r = 0; r < rows; ++r) {
            row_starts[r + 1] = row_starts[r] + row_counts[r];
        }
        memset(row_counts, 0, (size_t)rows * sizeof(int));

        int *row_cols = (int*)malloc((nnz > 0 ? nnz : 1) * sizeof(int));
        if (!row_cols) {
            free(col_masks);
            free(row_counts);
            free(row_starts);
            free(priority);
            break;
        }

        for (int i = 0; i < ncore; ++i) {
            int c = col_map[i];
            const int *colv = pichart + (size_t)c * (size_t)rows;
            for (int r = 0; r < rows; ++r) {
                if (colv[r]) {
                    row_cols[row_starts[r] + row_counts[r]] = i;
                    row_counts[r]++;
                }
            }
        }

        free(sol01);
        sol01 = (int*)malloc((size_t)ncore * sizeof(int));
        if (!sol01) {
            free(col_masks);
            free(row_counts);
            free(row_starts);
            free(row_cols);
            free(priority);
            break;
        }

        qca_scp_problem problem = {
            .nr = rows,
            .nc = ncore,
            .nwords_rows = nwords_rows,
            .row_starts = row_starts,
            .row_cols = row_cols,
            .col_masks = col_masks,
            .branch_priority = priority
        };

        int size = 0;
        qca_scp_result exact_result = qca_scp_solve_exact_with_incumbent(
            &problem,
            sol01,
            &size,
            NULL,
            0,
            target
        );

        free(col_masks);
        free(row_counts);
        free(row_starts);
        free(row_cols);
        free(priority);

        if (exact_result == QCA_SCP_SOLUTION && size > 0 && size <= target) {
            /* improvement: map the core solution back to original columns */
            int out = 0;
            for (int i = 0; i < ncore; ++i) {
                if (sol01[i]) solution[out++] = col_map[i];
            }
            *solmin = out;
            UB = out;
            continue;
        }

        if (exact_result == QCA_SCP_ERROR) {
            break;
        }

        if (exact_result == QCA_SCP_NO_SOLUTION) {
            /* Exhausted without a cover of size <= target: UB is optimal. */
            proved = 1;
            break;
        }

        /* A successful proof-target solve must have returned above. Treat
           any other result conservatively rather than claiming a proof. */
        break;
    }

    if (UB <= 1) proved = 1; /* a single-column cover is trivially optimal */

    if (proved && bestLB_io && (double)UB > *bestLB_io) {
        *bestLB_io = (double)UB;
    }

    free(col_map);
    free(sol01);
    return proved;
}

static int lagr_max_int(int a, int b) {
    return a > b ? a : b;
}

static int lagr_normalize_effort_level(int effort_level) {
    if (effort_level < 0) return 0;
    if (effort_level > 3) return 3;
    return effort_level;
}

static void lagr_config_for_effort(LagrangianConfig *cfg, int effort_level) {
    effort_level = lagr_normalize_effort_level(effort_level);

    cfg->max_iter = 200;           /* total iterations of subgradient */
    cfg->heur_every = 3;           /* construct a feasible solution every N iterations */
    cfg->step_coef = 1.5;          /* initial step coefficient */
    cfg->step_min = 0.005;         /* minimal step coefficient */
    cfg->step_contract = 0.5;      /* current CCubes behavior: halve after stagnation */
    cfg->stabilization_beta = 1.0; /* no damping unless explicitly requested */
    cfg->halve_period = 8;         /* halve step if no progress for this many iterations */
    cfg->deflection_alpha = 0.0;   /* deflected subgradient direction, disabled by default */
    cfg->polish_enabled = 1;       /* reduced-cost fixing + bounded exact finish */
    cfg->local_passes = 0;
    cfg->drop2_enabled = 0;
    cfg->small_gap_threshold = 0;
    cfg->small_gap_extra_passes = cfg->local_passes;
    cfg->rarity_init = 0;
    cfg->polish_node_limit = 80000;
    cfg->polish_core_cap = 1500;
    cfg->portfolio_enabled = 0;
    cfg->portfolio_max_profiles = 1;
    cfg->portfolio_deflection_alpha = 0.75;
    cfg->hybrid_bundle_portfolio = 0;
    cfg->bundle_enabled = 0;
    cfg->bundle_interval = 1;
    cfg->bundle_memory = 8;
    cfg->bundle_prox_mu = 8.0;
    cfg->bundle_min_mu = 1.0;
    cfg->bundle_max_mu = 1024.0;
    cfg->bundle_null_expand = 2.0;
    cfg->bundle_serious_shrink = 0.8;
    cfg->bundle_serious_fraction = 0.05;
    cfg->bundle_serious_tol = 1e-9;
    cfg->bundle_momentum = 0.0;

    if (effort_level == 1) {
        cfg->max_iter = 1000;
        cfg->heur_every = 5;
        cfg->step_coef = 2.0;
        cfg->local_passes = 10;
        cfg->drop2_enabled = 0;
        cfg->small_gap_extra_passes = cfg->local_passes;
        cfg->step_contract = 0.95;
        cfg->polish_enabled = 1;
        cfg->polish_node_limit = 400000;
        cfg->polish_core_cap = 8000;
        cfg->portfolio_enabled = 1;
        cfg->portfolio_max_profiles = 2;
    }

    if (effort_level >= 2) {
        /* level 3 runs the same Lagrangian machinery as level 2 and then
           completes any uncertified chart with the exact scp_solver */
        cfg->max_iter = 5000;
        cfg->heur_every = 10;
        cfg->step_coef = 2.0;
        cfg->local_passes = 10;
        cfg->drop2_enabled = 1;
        cfg->small_gap_extra_passes = cfg->local_passes;
        cfg->step_contract = 0.95;
        cfg->polish_enabled = 1;
        cfg->polish_node_limit = 2000000;
        cfg->polish_core_cap = 30000;
        cfg->portfolio_enabled = 1;
        cfg->portfolio_max_profiles = 6;
        cfg->hybrid_bundle_portfolio = 1;
        cfg->bundle_interval = 8;
        cfg->bundle_momentum = 0.35;
    }

    if (effort_level == 3) {
        /* Exact completion is expensive enough that the fuller portfolio
           pays for itself when it closes the gap before branch-and-bound.
           Keep level 2 lean; level 3 is the explicit exactness choice. */
        cfg->local_passes = 10;
        cfg->small_gap_extra_passes = cfg->local_passes;
        cfg->step_contract = 0.95;
        cfg->portfolio_max_profiles = 6;
    }
}

static int lagr_build_portfolio_configs(
    const LagrangianConfig *baseline,
    double deflection_alpha,
    LagrangianConfig profiles[LAGR_PORTFOLIO_PROFILE_COUNT]
) {
    profiles[0] = *baseline;

    if (baseline->hybrid_bundle_portfolio) {
        LagrangianConfig bundle = *baseline;
        bundle.bundle_enabled = 1;
        bundle.deflection_alpha = 0.0;

        profiles[1] = *baseline;
        profiles[1].deflection_alpha = deflection_alpha;

        profiles[2] = bundle;
        profiles[2].deflection_alpha = 0.0;

        profiles[3] = bundle;
        profiles[3].deflection_alpha = deflection_alpha;

        profiles[4] = *baseline;
        profiles[4].max_iter = 500;
        profiles[4].step_contract = 0.99;
        profiles[4].deflection_alpha = 0.0;

        profiles[5] = profiles[4];
        profiles[5].deflection_alpha = deflection_alpha;

        return LAGR_PORTFOLIO_PROFILE_COUNT;
    }

    profiles[1] = *baseline;
    profiles[1].deflection_alpha = deflection_alpha;

    profiles[2] = *baseline;
    profiles[2].max_iter = lagr_max_int(baseline->max_iter, 80000);
    profiles[2].halve_period = lagr_max_int(baseline->halve_period, 16);

    profiles[3] = profiles[2];
    profiles[3].deflection_alpha = deflection_alpha;

    return 4;
}

/* ---------- Main functions ---------- */
static void solve_scp_lagrangian_config(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double weights[],
    const LagrangianConfig *cfg,
    const int *initial_solution,
    int initial_solmin,
    int *solution,
    int *solmin,
    double *ls_out, /* optional: reduced costs at the best multipliers */
    double *row_dual_out
) {
    if (solmin) *solmin = -1;
    lagr_stats_begin(rows, cols, 0);
    if (
        rows <= 0 ||
        cols <= 0 ||
        !rowsCovered ||
        !rowsCoveredCount ||
        !colsCovering ||
        !colsCoveringCount ||
        !cfg ||
        !solution ||
        !solmin
    ) {
        lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_NO_CANDIDATE);
        return;
    }

    double *t = (double*)calloc((size_t)rows, sizeof(double));
    double *t_best = (double*)calloc((size_t)rows, sizeof(double));
    double *ls = (double*)malloc((size_t)cols * sizeof(double));
    int *sol_tmp = (int*)malloc((size_t)cols * sizeof(int));
    int *sol_tmp2 = (int*)malloc((size_t)cols * sizeof(int));
    double *col_costs = (double*)malloc((size_t)cols * sizeof(double));
    unsigned char *dual_x = (unsigned char*)malloc((size_t)cols);
    int best_sol_size = -1;
    int stuck_iter = 0;

    if (!t || !t_best || !ls || !sol_tmp || !sol_tmp2 || !col_costs || !dual_x) {
        free(t);
        free(t_best);
        free(ls);
        free(sol_tmp);
        free(sol_tmp2);
        free(col_costs);
        free(dual_x);

        lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_OOM);
        return;
    }

    double max_w = 0.0;
    if (weights) {
        for (int j = 0; j < cols; ++j) {
            if (weights[j] > max_w) max_w = weights[j];
        }
    }
    double eps = 0.0;
    if (max_w > 0.0) {
        eps = 1e-6 / (max_w + 1.0);
    }
    for (int j = 0; j < cols; ++j) {
        col_costs[j] = 1.0 - eps * (weights ? weights[j] : 0.0);
    }

    int local_passes = cfg->local_passes;
    int max_iter = cfg->max_iter;
    int heur_every = cfg->heur_every;
    double step_coef = cfg->step_coef;
    double step_min = cfg->step_min;
    double step_contract = cfg->step_contract;
    double stabilization_beta = cfg->stabilization_beta;
    int halve_period = cfg->halve_period;
    double deflection_alpha = cfg->deflection_alpha;
    double *prev_direction = NULL;
    LagrangianBundleState bundle_state = {0};

    if (deflection_alpha > 0.0) {
        prev_direction = (double*)calloc((size_t)rows, sizeof(double));
        if (!prev_direction) {
            free(t);
            free(t_best);
            free(ls);
            free(sol_tmp);
            free(sol_tmp2);
            free(col_costs);
            free(dual_x);
            lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_OOM);
            return;
        }
    }

    if (cfg->bundle_enabled) {
        int memory = cfg->bundle_memory > 0 ? cfg->bundle_memory : 1;
        bundle_state.rows = rows;
        bundle_state.memory = memory;
        bundle_state.prox_mu = cfg->bundle_prox_mu;
        bundle_state.min_mu = cfg->bundle_min_mu;
        bundle_state.max_mu = cfg->bundle_max_mu;
        bundle_state.null_expand = cfg->bundle_null_expand;
        bundle_state.serious_shrink = cfg->bundle_serious_shrink;
        bundle_state.serious_fraction = cfg->bundle_serious_fraction;
        bundle_state.serious_tol = cfg->bundle_serious_tol;
        bundle_state.momentum = cfg->bundle_momentum;
        bundle_state.center_value = -DBL_MAX;
        bundle_state.center_t = (double*)calloc((size_t)rows, sizeof(double));
        bundle_state.prev_center_t = (double*)calloc((size_t)rows, sizeof(double));
        bundle_state.cut_alpha = (double*)calloc((size_t)memory, sizeof(double));
        bundle_state.cut_slacks = (double*)calloc((size_t)memory * (size_t)rows, sizeof(double));
        if (!bundle_state.center_t || !bundle_state.prev_center_t || !bundle_state.cut_alpha || !bundle_state.cut_slacks) {
            free(bundle_state.center_t);
            free(bundle_state.prev_center_t);
            free(bundle_state.cut_alpha);
            free(bundle_state.cut_slacks);
            free(prev_direction);
            free(t);
            free(t_best);
            free(ls);
            free(sol_tmp);
            free(sol_tmp2);
            free(col_costs);
            free(dual_x);
            lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_OOM);
            return;
        }
    }

    lagr_initialize_multipliers(rows, colsCoveringCount, cfg->rarity_init, t);

    double bestZLB = -DBL_MAX;
    double lastZLB = -DBL_MAX;
    int iterations = 0;
    LagrangianStopReason stop_reason = LAGR_STOP_MAX_ITER;

    /* Initial heuristics with current ls (computed from t) */
    double initialZLB = compute_reduced_and_lb(
        rows,
        cols,
        rowsCovered,
        rowsCoveredCount,
        t,
        col_costs,
        ls,
        dual_x
    );
    bestZLB = initialZLB;
    lastZLB = initialZLB;
    memcpy(t_best, t, (size_t)rows * sizeof(double));

    int sol_size = -1, sol_size2 = -1;

    heuristic_row_min_rc(
        rows,
        cols,
        rowsCovered,
        rowsCoveredCount,
        colsCovering,
        colsCoveringCount,
        ls,
        weights,
        sol_tmp,
        &sol_size
    );

    greedy_from_lagr_scores(
        rows,
        cols,
        rowsCovered,
        rowsCoveredCount,
        colsCovering,
        colsCoveringCount,
        ls,
        weights,
        sol_tmp2,
        &sol_size2
    );

    /* Optional deeper search: disabled by default */
    if (local_passes > 0) {
        if (sol_size > 0) {
            local_search_drop1_repair(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp,
                &sol_size,
                local_passes
            );

            local_search_drop2_repair(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp,
                &sol_size,
                local_passes
            );

        }

        if (sol_size2 > 0) {
            local_search_drop1_repair(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp2,
                &sol_size2,
                local_passes
            );

            local_search_drop2_repair(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp2,
                &sol_size2,
                local_passes
            );
        }
    }

    int have_initial = initial_solution != NULL && initial_solmin > 0 && initial_solmin <= cols;

    if (sol_size == -1 && sol_size2 == -1 && !have_initial) {
        free(t);
        free(t_best);
        free(ls);
        free(sol_tmp);
        free(sol_tmp2);
        free(col_costs);
        free(dual_x);
        free(prev_direction);
        free(bundle_state.center_t);
        free(bundle_state.prev_center_t);
        free(bundle_state.cut_alpha);
        free(bundle_state.cut_slacks);
        lagr_stats_finish(-1, -DBL_MAX, bestZLB, lastZLB, step_coef, 0, LAGR_STOP_NO_CANDIDATE);
        return;
    }

    double w1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
    double w2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;

    double bestUB = DBL_MAX; /* UB equals current solution length (unit-cost) */
    int bestTerms = INT_MAX;
    double bestWeight = -DBL_MAX;

    if (sol_size != -1) {
        bestUB = (double)sol_size;
        bestTerms = sol_size;
        bestWeight = w1;
        best_sol_size = sol_size;
        memcpy(solution, sol_tmp, (size_t)sol_size * sizeof(int));
    }

    if (
        sol_size2 != -1 &&
        (
            sol_size2 < bestTerms ||
            (sol_size2 == bestTerms && w2 > bestWeight + EPS)
        )
    ) {
        bestUB = (double)sol_size2;
        bestTerms = sol_size2;
        bestWeight = w2;
        best_sol_size = sol_size2;
        memcpy(solution, sol_tmp2, (size_t)sol_size2 * sizeof(int));
    }

    if (have_initial) {
        double initial_weight = solution_total_weight(initial_solution, initial_solmin, weights);
        if (
            initial_solmin < bestTerms ||
            (initial_solmin == bestTerms && initial_weight > bestWeight + EPS)
        ) {
            bestUB = (double)initial_solmin;
            bestTerms = initial_solmin;
            bestWeight = initial_weight;
            best_sol_size = initial_solmin;
            memcpy(solution, initial_solution, (size_t)initial_solmin * sizeof(int));
        }
    }

    double bestLB = -DBL_MAX;
    for (int it = 0; it < max_iter; ++it) {
        iterations = it + 1;
        double ZLB = compute_reduced_and_lb(
            rows,
            cols,
            rowsCovered,
            rowsCoveredCount,
            t,
            col_costs,
            ls,
            dual_x
        );
        lastZLB = ZLB;
        if (ZLB > bestZLB + EPS) {
            bestZLB = ZLB;
            memcpy(t_best, t, (size_t)rows * sizeof(double));
        }

        double LBint = ceil(ZLB - 1e-12);

        if (LBint > bestLB + EPS) {
            bestLB = LBint;
            stuck_iter = 0; /* reset halving on improvement */
        } /* else no improvement */

        if (it % heur_every == 0) {
            heuristic_row_min_rc(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp,
                &sol_size
            );

            greedy_from_lagr_scores(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp2,
                &sol_size2
            );

            bool found_new_best = false;
            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;

            /* Improve by primary objective: fewer terms; tie-break by higher weight */
            if (
                sol_size != -1 &&
                (
                    sol_size < bestTerms ||
                    (
                        sol_size == bestTerms && candW1 > bestWeight + EPS
                    )
                )
            ) {
                bestUB = (double)sol_size;
                bestTerms = sol_size;
                bestWeight = candW1;
                best_sol_size = sol_size;
                memcpy(solution, sol_tmp, (size_t)sol_size * sizeof(int));
                found_new_best = true;
            }

            if (
                sol_size2 != -1 &&
                (
                    sol_size2 < bestTerms ||
                    (
                        sol_size2 == bestTerms && candW2 > bestWeight + EPS
                    )
                )
            ) {
                bestUB = (double)sol_size2;
                bestTerms = sol_size2;
                bestWeight = candW2;
                best_sol_size = sol_size2;
                memcpy(solution, sol_tmp2, (size_t)sol_size2 * sizeof(int));
                found_new_best = true;
            }

            if (found_new_best && local_passes > 0) {
                int current_sol_size = best_sol_size;
                int *current_sol = (int*)malloc((size_t)cols * sizeof(int));
                if (current_sol) {
                    memcpy(current_sol, solution, (size_t)current_sol_size * sizeof(int));

                    local_search_drop1_repair(
                        rows,
                        cols,
                        rowsCovered,
                        rowsCoveredCount,
                        colsCovering,
                        colsCoveringCount,
                        ls,
                        weights,
                        current_sol,
                        &current_sol_size,
                        local_passes
                    );

                    if (cfg->drop2_enabled) {
                        local_search_drop2_repair(
                            rows,
                            cols,
                            rowsCovered,
                            rowsCoveredCount,
                            colsCovering,
                            colsCoveringCount,
                            ls,
                            weights,
                            current_sol,
                            &current_sol_size,
                            local_passes
                        );
                    }

                    double polished_weight = solution_total_weight(current_sol, current_sol_size, weights);
                    if (
                        current_sol_size < bestTerms ||
                        (
                            current_sol_size == bestTerms && polished_weight > bestWeight + EPS
                        )
                    ) {
                        bestUB = (double)current_sol_size;
                        bestTerms = current_sol_size;
                        bestWeight = polished_weight;
                        best_sol_size = current_sol_size;
                        memcpy(solution, current_sol, (size_t)current_sol_size * sizeof(int));
                    }
                    free(current_sol);
                }
            }
        }

        if (bestUB <= LBint + EPS) {
            stop_reason = LAGR_STOP_OPTIMAL;
            break; /* proved optimal (within epsilon) */
        }

        int used_subgradient_update = 0;
        int updated = 0;

        int use_bundle_update = cfg->bundle_enabled &&
            (cfg->bundle_interval <= 1 || it % cfg->bundle_interval == 0);

        if (use_bundle_update) {
            updated = bundle_update(
                rows,
                colsCovering,
                colsCoveringCount,
                dual_x,
                bestUB,
                ZLB,
                t,
                &bundle_state
            );
        }

        if (!use_bundle_update || !updated) {
            updated = subgradient_update(
                rows,
                colsCovering,
                colsCoveringCount,
                dual_x,
                bestUB,
                ZLB,
                t,
                &step_coef,
                step_min,
                &stuck_iter,
                halve_period,
                prev_direction,
                deflection_alpha,
                step_contract,
                stabilization_beta
            );
            used_subgradient_update = 1;
        }

        if (!updated || (used_subgradient_update && step_coef <= step_min + EPS)) {
            stop_reason = updated ? LAGR_STOP_STEP_MIN : LAGR_STOP_NO_UPDATE;
            /* final refresh */
            heuristic_row_min_rc(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp,
                &sol_size
            );

            greedy_from_lagr_scores(
                rows,
                cols,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                ls,
                weights,
                sol_tmp2,
                &sol_size2
            );

            bool found_new_best = false;
            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;

            if (
                sol_size != -1 &&
                (
                    sol_size < bestTerms ||
                    (
                        sol_size == bestTerms && candW1 > bestWeight + EPS
                    )
                )
            ) {
                bestUB = (double)sol_size;
                bestTerms = sol_size;
                bestWeight = candW1;
                best_sol_size = sol_size;
                memcpy(solution, sol_tmp, (size_t)sol_size * sizeof(int));
                found_new_best = true;
            }

            if (
                sol_size2 != -1 &&
                (
                    sol_size2 < bestTerms ||
                    (
                        sol_size2 == bestTerms && candW2 > bestWeight + EPS
                    )
                )
            ) {
                bestUB = (double)sol_size2;
                bestTerms = sol_size2;
                bestWeight = candW2;
                best_sol_size = sol_size2;
                memcpy(solution, sol_tmp2, (size_t)sol_size2 * sizeof(int));
                found_new_best = true;
            }

            if (found_new_best && local_passes > 0) {
                int current_sol_size = best_sol_size;
                int *current_sol = (int*)malloc((size_t)cols * sizeof(int));
                if (current_sol) {
                    memcpy(current_sol, solution, (size_t)current_sol_size * sizeof(int));

                    local_search_drop1_repair(
                        rows,
                        cols,
                        rowsCovered,
                        rowsCoveredCount,
                        colsCovering,
                        colsCoveringCount,
                        ls,
                        weights,
                        current_sol,
                        &current_sol_size,
                        local_passes
                    );

                    if (cfg->drop2_enabled) {
                        local_search_drop2_repair(
                            rows,
                            cols,
                            rowsCovered,
                            rowsCoveredCount,
                            colsCovering,
                            colsCoveringCount,
                            ls,
                            weights,
                            current_sol,
                            &current_sol_size,
                            local_passes
                        );
                    }

                    double polished_weight = solution_total_weight(current_sol, current_sol_size, weights);
                    if (
                        current_sol_size < bestTerms ||
                        (
                            current_sol_size == bestTerms && polished_weight > bestWeight + EPS
                        )
                    ) {
                        bestUB = (double)current_sol_size;
                        bestTerms = current_sol_size;
                        bestWeight = polished_weight;
                        best_sol_size = current_sol_size;
                        memcpy(solution, current_sol, (size_t)current_sol_size * sizeof(int));
                    }
                    free(current_sol);
                }
            }
            break;
        }
    }

    if (
        cfg->polish_enabled &&
        best_sol_size > 1 &&
        stop_reason != LAGR_STOP_OPTIMAL
    ) {
        /* reduced costs at the best multipliers found */
        double zlb_best = compute_reduced_and_lb(
            rows,
            cols,
            rowsCovered,
            rowsCoveredCount,
            t_best,
            col_costs,
            ls,
            dual_x
        );

        int proved = 0;
        lagr_exact_finish(
            rows,
            cols,
            rowsCovered,
            rowsCoveredCount,
            colsCovering,
            colsCoveringCount,
            weights,
            ls,
            zlb_best,
            cfg->polish_node_limit,
            cfg->polish_core_cap,
            solution,
            &best_sol_size,
            &bestLB,
            &proved,
            1, NULL, NULL, NULL, NULL
        );

        if ((double)best_sol_size < bestUB - EPS) {
            bestUB = (double)best_sol_size;
            bestTerms = best_sol_size;
            bestWeight = solution_total_weight(solution, best_sol_size, weights);
        }

        if (proved) {
            stop_reason = LAGR_STOP_OPTIMAL;
        }
    }

    if (ls_out) {
        /* reduced costs at the best multipliers, for reuse by the caller
           (branching priorities of the exact completion) */
        compute_reduced_and_lb(
            rows,
            cols,
            rowsCovered,
            rowsCoveredCount,
            t_best,
            col_costs,
            ls_out,
            dual_x
        );
    }

    if (row_dual_out) memcpy(row_dual_out, t_best, (size_t)rows * sizeof(double));
    *solmin = best_sol_size;
    lagr_stats_finish(best_sol_size, bestLB, bestZLB, lastZLB, step_coef, iterations, stop_reason);

    free(t);
    free(t_best);
    free(ls);
    free(sol_tmp);
    free(sol_tmp2);
    free(col_costs);
    free(dual_x);
    free(prev_direction);
    free(bundle_state.center_t);
    free(bundle_state.prev_center_t);
    free(bundle_state.cut_alpha);
    free(bundle_state.cut_slacks);
}

static void solvePIchart_lagrangian_level(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out,
    int effort_level,
    unsigned char *improving_core_out,
    int *improving_core_size_out,
    const int *initial_solution,
    int initial_solmin,
    double *row_dual_out
) {
    if (row_dual_out && ON_minterms > 0) {
        memset(row_dual_out, 0, (size_t)ON_minterms * sizeof(double));
    }
    if (solmin) *solmin = -1;
    if (best_lb_out) *best_lb_out = -DBL_MAX;
    if (improving_core_size_out) *improving_core_size_out = 0;
    if (improving_core_out && foundPI > 0) {
        memset(improving_core_out, 0, (size_t)foundPI);
    }
    if (lagr_score_out && foundPI > 0) {
        memset(lagr_score_out, 0, (size_t) foundPI * sizeof(double));
    }
    lagr_stats_begin(ON_minterms, foundPI, 0);
    if (!pichart || foundPI <= 0 || ON_minterms <= 0 || !solution || !solmin) {
        lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_NO_CANDIDATE);
        return;
    }
    int **rowsCovered = NULL, *rowsCoveredCount = NULL;
    int **colsCovering = NULL, *colsCoveringCount = NULL;

    int rc_ad = build_adjacency(
        pichart,
        ON_minterms,
        foundPI,
        &rowsCovered,
        &rowsCoveredCount,
        &colsCovering,
        &colsCoveringCount
    );

    if (rc_ad == -2) { /* infeasible matrix (some row uncovered) */
        lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_INFEASIBLE);
        return;
    }

    if (rc_ad == -1) {
        lagr_stats_finish(-1, -DBL_MAX, -DBL_MAX, -DBL_MAX, 0.0, 0, LAGR_STOP_OOM);
        return;
    }

    /* shrink the instance once; benefits every portfolio profile */
    presolve_dominated_columns(
        ON_minterms,
        foundPI,
        rowsCovered,
        rowsCoveredCount,
        colsCovering,
        colsCoveringCount,
        weights
    );

    LagrangianConfig baseline;
    lagr_config_for_effort(&baseline, effort_level);
    if (improving_core_out && effort_level == 2) {
        /* The hybrid wrapper has a separate exact finisher. Keep this pass
           focused on bounds, an incumbent and a useful reduced-cost core. */
        baseline.portfolio_max_profiles = 2;
        baseline.polish_node_limit = 50000;
    }
    int portfolio_enabled = baseline.portfolio_enabled;
    double portfolio_deflection_alpha = baseline.portfolio_deflection_alpha;
    int want_exact = lagr_normalize_effort_level(effort_level) >= 3;
    int need_scores = want_exact || lagr_score_out || improving_core_out;

    if (portfolio_enabled) {
        baseline.deflection_alpha = 0.0;
    }

    /* reduced costs paired with their dual value, for the exact completion */
    double *scores = need_scores ? (double*)malloc((size_t)foundPI * sizeof(double)) : NULL;
    double *candidate_scores = NULL;
    double *candidate_dual = NULL;
    double scores_zlb = -DBL_MAX;
    int *candidate = NULL;

    if (need_scores && !scores) {
        want_exact = 0;
        need_scores = 0; /* OOM: continue without a reduced-cost handoff */
    }

    solve_scp_lagrangian_config(
        ON_minterms,
        foundPI,
        rowsCovered,
        rowsCoveredCount,
        colsCovering,
        colsCoveringCount,
        weights,
        &baseline,
        initial_solution,
        initial_solmin,
        solution,
        solmin,
        scores,
        row_dual_out
    );

    LagrangianStats merged_stats = *lagrangian_last_stats();
    int best_solmin = *solmin;
    int total_iterations = merged_stats.iterations;
    scores_zlb = merged_stats.best_zlb;

    if (
        best_solmin < 0 ||
        !portfolio_enabled ||
        (
            merged_stats.gap == 0 &&
            merged_stats.stop_reason == LAGR_STOP_OPTIMAL
        )
    ) {
        goto completion;
    }

    {
        LagrangianConfig profiles[LAGR_PORTFOLIO_PROFILE_COUNT];
        int profile_count = lagr_build_portfolio_configs(&baseline, portfolio_deflection_alpha, profiles);
        int max_profiles = baseline.portfolio_max_profiles;
        if (profile_count > max_profiles) profile_count = max_profiles;

        if (profile_count <= 1) {
            goto completion;
        }

        candidate = (int*)malloc((size_t)foundPI * sizeof(int));
        if (need_scores) {
            candidate_scores = (double*)malloc((size_t)foundPI * sizeof(double));
        }
        if (row_dual_out) candidate_dual = (double*)calloc((size_t)ON_minterms, sizeof(double));
        if (!candidate || (need_scores && !candidate_scores) ||
            (row_dual_out && !candidate_dual)) {
            goto completion;
        }

        for (int i = 1; i < profile_count; ++i) {
            R_CheckUserInterrupt();
            int candidate_solmin = -1;

            solve_scp_lagrangian_config(
                ON_minterms,
                foundPI,
                rowsCovered,
                rowsCoveredCount,
                colsCovering,
                colsCoveringCount,
                weights,
                &profiles[i],
                solution,
                best_solmin,
                candidate,
                &candidate_solmin,
                candidate_scores,
                candidate_dual
            );

            const LagrangianStats *candidate_stats = lagrangian_last_stats();
            total_iterations += candidate_stats->iterations;

            if (candidate_stats->best_lb > merged_stats.best_lb) {
                merged_stats.best_lb = candidate_stats->best_lb;
                merged_stats.last_zlb = candidate_stats->last_zlb;
                merged_stats.step_coef = candidate_stats->step_coef;
                merged_stats.stop_reason = candidate_stats->stop_reason;
            }
            if (candidate_stats->best_zlb > merged_stats.best_zlb + EPS) {
                merged_stats.best_zlb = candidate_stats->best_zlb;
                merged_stats.last_zlb = candidate_stats->last_zlb;
                merged_stats.step_coef = candidate_stats->step_coef;
            }
            if (candidate_scores && candidate_stats->best_zlb > scores_zlb + EPS) {
                /* keep (zlb, reduced costs) from the same run: the pair
                   drives the reduced-cost fixing of the completion */
                memcpy(scores, candidate_scores, (size_t)foundPI * sizeof(double));
                scores_zlb = candidate_stats->best_zlb;
                if (row_dual_out) memcpy(row_dual_out, candidate_dual, (size_t)ON_minterms * sizeof(double));
            }

            if (candidate_solmin >= 0 && candidate_solmin < best_solmin) {
                best_solmin = candidate_solmin;
                *solmin = candidate_solmin;
                memcpy(solution, candidate, (size_t)candidate_solmin * sizeof(int));
            }

            if (merged_stats.best_lb != INT_MIN && best_solmin <= merged_stats.best_lb) {
                break;
            }
        }

        LagrangianStopReason final_reason = merged_stats.stop_reason;
        if (merged_stats.best_lb != INT_MIN && best_solmin <= merged_stats.best_lb) {
            final_reason = LAGR_STOP_OPTIMAL;
        }

        lagr_stats_finish(
            best_solmin,
            merged_stats.best_lb == INT_MIN ? -DBL_MAX : (double)merged_stats.best_lb,
            merged_stats.best_zlb,
            merged_stats.last_zlb,
            merged_stats.step_coef,
            total_iterations,
            final_reason
        );
    }

completion:
    /* -l3: decide optimality with the exact scp_solver when uncertified */
    if (want_exact && best_solmin > 0) {
        const LagrangianStats *st = lagrangian_last_stats();
        int certified =
            st->gap == 0 &&
            (
                st->stop_reason == LAGR_STOP_OPTIMAL ||
                st->stop_reason == LAGR_STOP_EXACT_COMPLETION
            );

        if (!certified && scores_zlb > -DBL_MAX) {
            double bestLB_dbl = (st->best_lb == INT_MIN) ? -DBL_MAX : (double)st->best_lb;

            int proved = lagr_exact_complete(
                pichart,
                ON_minterms,
                foundPI,
                rowsCoveredCount,
                scores,
                scores_zlb,
                solution,
                solmin,
                &bestLB_dbl
            );

            if (proved) {
                best_solmin = *solmin;
                lagr_stats_finish(
                    best_solmin,
                    (double)best_solmin,
                    st->best_zlb,
                    st->last_zlb,
                    st->step_coef,
                    total_iterations,
                    LAGR_STOP_EXACT_COMPLETION
                );
            }
        }
    }

    {
        const LagrangianStats *final_stats = lagrangian_last_stats();
        if (best_lb_out && final_stats->best_lb != INT_MIN) {
            *best_lb_out = (double) final_stats->best_lb;
        }
        if (lagr_score_out && scores) {
            memcpy(lagr_score_out, scores, (size_t) foundPI * sizeof(double));
        }
        if (improving_core_out) {
            int core_size = 0;
            int target = best_solmin - 1;

            for (int c = 0; c < foundPI; ++c) {
                int allowed = rowsCoveredCount[c] > 0;
                if (
                    allowed &&
                    scores &&
                    scores_zlb > -DBL_MAX &&
                    target >= 0 &&
                    scores_zlb + (scores[c] > 0.0 ? scores[c] : 0.0) >
                        (double)target + 1e-7
                ) {
                    allowed = 0;
                }
                improving_core_out[c] = (unsigned char)allowed;
                core_size += allowed;
            }
            if (improving_core_size_out) {
                *improving_core_size_out = core_size;
            }
        }
    }

    free(candidate);
    free(candidate_scores);
    free(candidate_dual);
    free(scores);
    free_adjacency(
        rowsCovered,
        rowsCoveredCount,
        foundPI,
        colsCovering,
        colsCoveringCount,
        ON_minterms
    );
}

void solvePIchart_lagrangian(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out
) {
    solvePIchart_lagrangian_level(
        pichart, foundPI, ON_minterms, weights,
        solution, solmin, best_lb_out, lagr_score_out,
        2, NULL, NULL, NULL, 0, NULL
    );
}

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
) {
    solvePIchart_lagrangian_level(
        pichart, foundPI, ON_minterms, weights,
        solution, solmin, best_lb_out, lagr_score_out,
        2, improving_core_out, improving_core_size_out, NULL, 0, NULL
    );
}

void solvePIchart_lagrangian_prepare_with_incumbent(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    const int *initial_solution,
    int initial_solmin,
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out,
    unsigned char *improving_core_out,
    int *improving_core_size_out
) {
    solvePIchart_lagrangian_level(
        pichart, foundPI, ON_minterms, weights,
        solution, solmin, best_lb_out, lagr_score_out,
        2, improving_core_out, improving_core_size_out,
        initial_solution, initial_solmin, NULL
    );
}

void solvePIchart_lagrangian_prepare_native(
    int pichart[], int foundPI, int ON_minterms, const double weights[],
    const int *initial_solution, int initial_solmin, int *solution, int *solmin,
    double *best_lb_out, double *lagr_score_out,
    unsigned char *improving_core_out, int *improving_core_size_out,
    double *row_dual_out
) {
    solvePIchart_lagrangian_level(
        pichart, foundPI, ON_minterms, weights,
        solution, solmin, best_lb_out, lagr_score_out,
        2, improving_core_out, improving_core_size_out,
        initial_solution, initial_solmin, row_dual_out
    );
}

SEXP C_findminLagrangian(SEXP chart) {
  if (!isMatrix(chart)) {
    return ScalarInteger(0);
  }

  int rows = nrows(chart);
  int cols = ncols(chart);
  if (rows <= 0 || cols <= 0) {
    return ScalarInteger(0);
  }

  int *pichart = (int *)R_Calloc((size_t)rows * (size_t)cols, int);
  int *indices = (int *)R_Calloc((size_t)cols, int);
  if (!pichart || !indices) {
    if (pichart)
      R_Free(pichart);
    if (indices)
      R_Free(indices);
    error("Memory allocation failed in C_findminLagrangian().");
  }

  switch (TYPEOF(chart)) {
  case LGLSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = LOGICAL(chart)[i];
    }
    break;
  case INTSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = INTEGER(chart)[i];
    }
    break;
  case REALSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = REAL(chart)[i] > 0;
    }
    break;
  default:
    R_Free(pichart);
    R_Free(indices);
    return ScalarInteger(0);
  }

  int solmin = -1;
  solvePIchart_lagrangian(pichart, cols, rows, NULL, indices, &solmin, NULL,
                          NULL);

  R_Free(pichart);

  if (solmin <= 0) {
    R_Free(indices);
    return ScalarInteger(0);
  }

  SEXP solution = PROTECT(allocVector(REALSXP, cols));
  for (int i = 0; i < cols; ++i) {
    REAL(solution)[i] = 0.0;
  }

  for (int i = 0; i < solmin; ++i) {
    if (indices[i] >= 0 && indices[i] < cols) {
      REAL(solution)[indices[i]] = 1.0;
    }
  }

  R_Free(indices);
  UNPROTECT(1);
  return solution;
}

SEXP C_findminLagrangianInfo(SEXP chart) {
  if (!isMatrix(chart)) {
    return ScalarInteger(0);
  }

  int rows = nrows(chart);
  int cols = ncols(chart);
  if (rows <= 0 || cols <= 0) {
    return ScalarInteger(0);
  }

  int *pichart = (int *)R_Calloc((size_t)rows * (size_t)cols, int);
  int *indices = (int *)R_Calloc((size_t)cols, int);
  if (!pichart || !indices) {
    if (pichart)
      R_Free(pichart);
    if (indices)
      R_Free(indices);
    error("Memory allocation failed in C_findminLagrangianInfo().");
  }

  switch (TYPEOF(chart)) {
  case LGLSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = LOGICAL(chart)[i];
    }
    break;
  case INTSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = INTEGER(chart)[i];
    }
    break;
  case REALSXP:
    for (int i = 0; i < rows * cols; ++i) {
      pichart[i] = REAL(chart)[i] > 0;
    }
    break;
  default:
    R_Free(pichart);
    R_Free(indices);
    return ScalarInteger(0);
  }

  int solmin = -1;
  double best_lb = -DBL_MAX;
  solvePIchart_lagrangian(pichart, cols, rows, NULL, indices, &solmin, &best_lb,
                          NULL);

  SEXP out = PROTECT(allocVector(VECSXP, 3));
  SEXP names = PROTECT(allocVector(STRSXP, 3));
  SEXP solution = PROTECT(allocVector(REALSXP, cols));

  for (int i = 0; i < cols; ++i) {
    REAL(solution)[i] = 0.0;
  }
  if (solmin > 0) {
    for (int i = 0; i < solmin; ++i) {
      if (indices[i] >= 0 && indices[i] < cols) {
        REAL(solution)[indices[i]] = 1.0;
      }
    }
  }

  SET_STRING_ELT(names, 0, mkChar("solution"));
  SET_STRING_ELT(names, 1, mkChar("upper_bound"));
  SET_STRING_ELT(names, 2, mkChar("lower_bound"));
  SET_VECTOR_ELT(out, 0, solution);
  SET_VECTOR_ELT(out, 1, ScalarInteger(solmin > 0 ? solmin : 0));
  SET_VECTOR_ELT(out, 2,
                 ScalarReal(best_lb > -DBL_MAX / 2 ? best_lb : NA_REAL));
  setAttrib(out, R_NamesSymbol, names);

  R_Free(pichart);
  R_Free(indices);
  UNPROTECT(3);
  return out;
}
