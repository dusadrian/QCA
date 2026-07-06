/*
Copyright (c) 2016 - 2026, Adrian Dusa
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, in whole or in part, are permitted provided that the
following conditions are met:
    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * The names of its contributors may NOT be used to endorse or promote
      products derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL ADRIAN DUSA BE LIABLE FOR ANY
DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

#include "scp_relaxation.h"
#include <math.h>
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
