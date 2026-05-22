#include "solvePIchart_lagrangian.h"

#include "qca_rinternals.h"
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include "qca_rinternals.h"
#include <stdbool.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>

#define EPS 1e-12
#define INTERRUPT_EVERY 1024

static void maybe_check_user_interrupt(int iteration) {
    if (iteration > 0 && iteration % INTERRUPT_EVERY == 0) {
        R_CheckUserInterrupt();
    }
}

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

    rcc = (int*) R_Calloc((size_t)cols, int);
    crc = (int*) R_Calloc((size_t)rows, int);
    if (!rcc || !crc) goto oom;

    for (int c = 0; c < cols; ++c) {
        int cnt = 0;
        for (int r = 0; r < rows; ++r) {
            cnt += pichart[c * rows + r];
        }
        rcc[c] = cnt;
    }

    for (int r = 0; r < rows; ++r) {
        int cnt = 0;
        for (int c = 0; c < cols; ++c) {
            cnt += pichart[c * rows + r];
        }
        crc[r] = cnt;
    }

    for (int r = 0; r < rows; ++r) {
        if (crc[r] == 0) {
            R_Free(rcc);
            R_Free(crc);
            return -2;
        }
    }

    rc = (int**) R_Calloc((size_t)cols, int*);
    cr = (int**) R_Calloc((size_t)rows, int*);
    if (!rc || !cr) goto oom;

    for (int c = 0; c < cols; ++c) {
        rc[c] = (rcc[c] > 0) ? (int*) R_Calloc((size_t)rcc[c], int) : NULL;
        if (rcc[c] > 0 && !rc[c]) goto oom;
    }

    for (int r = 0; r < rows; ++r) {
        cr[r] = (crc[r] > 0) ? (int*) R_Calloc((size_t)crc[r], int) : NULL;
        if (crc[r] > 0 && !cr[r]) goto oom;
    }

    for (int c = 0; c < cols; ++c) {
        int k = 0;
        for (int r = 0; r < rows; ++r) {
            if (pichart[c * rows + r] != 0) {
                rc[c][k++] = r;
            }
        }
    }

    for (int r = 0; r < rows; ++r) {
        int k = 0;
        for (int c = 0; c < cols; ++c) {
            if (pichart[c * rows + r] != 0) {
                cr[r][k++] = c;
            }
        }
    }

    *rowsCovered = rc;
    *rowsCoveredCount = rcc;
    *colsCovering = cr;
    *colsCoveringCount = crc;
    return 0;

oom:
    if (rc) {
        for (int c = 0; c < cols; ++c) if (rc[c]) R_Free(rc[c]);
        R_Free(rc);
    }
    if (cr) {
        for (int r = 0; r < rows; ++r) if (cr[r]) R_Free(cr[r]);
        R_Free(cr);
    }
    if (rcc) R_Free(rcc);
    if (crc) R_Free(crc);
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
        for (int c = 0; c < cols; ++c) if (rowsCovered[c]) R_Free(rowsCovered[c]);
        R_Free(rowsCovered);
    }
    if (rowsCoveredCount) R_Free(rowsCoveredCount);

    if (colsCovering) {
        for (int r = 0; r < rows; ++r) if (colsCovering[r]) R_Free(colsCovering[r]);
        R_Free(colsCovering);
    }
    if (colsCoveringCount) R_Free(colsCoveringCount);
}

typedef struct {
    double phi;
    double phi_min;
    double phi_contract;
    double stabilization_beta;
    double deflection_alpha;
    double best_zlb;
    double prev_zlb;
    double prev_step;
    int stagnation_iter;
    int stagnation_period;
    int have_best;
    int have_prev;
    double *best_t;
    double *prev_t;
    double *prev_direction;
} lagr_subgradient_state;

static int lagr_state_init(lagr_subgradient_state *state, int rows) {
    if (!state || rows <= 0) return 0;

    state->phi = 2.0;
    state->phi_min = 0.005;
    state->phi_contract = 0.90;
    state->stabilization_beta = 1.0;
    state->deflection_alpha = 0.0;
    state->best_zlb = -DBL_MAX;
    state->prev_zlb = -DBL_MAX;
    state->prev_step = 0.0;
    state->stagnation_iter = 0;
    state->stagnation_period = 8;
    state->have_best = 0;
    state->have_prev = 0;
    state->best_t = (double*) R_Calloc((size_t)rows, double);
    state->prev_t = (double*) R_Calloc((size_t)rows, double);
    state->prev_direction = (double*) R_Calloc((size_t)rows, double);

    return (state->best_t != NULL && state->prev_t != NULL && state->prev_direction != NULL);
}

static void lagr_state_free(lagr_subgradient_state *state) {
    if (!state) return;
    if (state->best_t) R_Free(state->best_t);
    state->best_t = NULL;
    if (state->prev_t) R_Free(state->prev_t);
    state->prev_t = NULL;
    if (state->prev_direction) R_Free(state->prev_direction);
    state->prev_direction = NULL;
}

static void lagr_state_reset_progress(lagr_subgradient_state *state) {
    if (!state) return;
    state->stagnation_iter = 0;
}

static void lagr_state_contract_phi(lagr_subgradient_state *state) {
    if (!state) return;
    state->phi *= state->phi_contract;
    if (state->phi < state->phi_min) state->phi = state->phi_min;
}

static void lagr_state_note_best_dual(
    lagr_subgradient_state *state,
    const double *t,
    int rows,
    double zlb
) {
    if (!state || !state->best_t || !t) return;
    if (!state->have_best || zlb > state->best_zlb + EPS) {
        Memcpy(state->best_t, t, rows);
        state->best_zlb = zlb;
        state->have_best = 1;
    }
}

static void prune_redundancy(
    int rows,
    int **rowsCovered,
    int *rowsCoveredCount,
    int *sol,
    int *sol_len
) {
    int n = *sol_len;
    if (n <= 0) return;

    int *coverCount = (int*) R_Calloc((size_t)rows, int);
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
    R_Free(coverCount);
}

static void greedy_from_lagr_scores(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int *sol,
    int *sol_len
) {
    (void)colsCovering;
    (void)colsCoveringCount;

    bool *rowCovered = (bool*) R_Calloc((size_t)rows, bool);
    bool *colSelected = (bool*) R_Calloc((size_t)cols, bool);
    int covered = 0, out = 0;

    if (!rowCovered || !colSelected) {
        if (rowCovered) R_Free(rowCovered);
        if (colSelected) R_Free(colSelected);
        *sol_len = -1;
        return;
    }

    while (covered < rows) {
        maybe_check_user_interrupt(covered);
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
                }
                else if (fabs(ls - bestLS) <= EPS) {
                    if (w > bestW) {
                        better = true;
                    }
                    else if (fabs(w - bestW) <= EPS) {
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

        if (best == -1) {
            R_Free(rowCovered);
            R_Free(colSelected);
            *sol_len = -1;
            return;
        }

        colSelected[best] = true;
        sol[out++] = best;

        for (int k = 0; k < rowsCoveredCount[best]; ++k) {
            int r = rowsCovered[best][k];
            if (!rowCovered[r]) {
                rowCovered[r] = true;
                covered++;
            }
        }
    }

    *sol_len = out;
    prune_redundancy(rows, rowsCovered, rowsCoveredCount, sol, sol_len);

    R_Free(rowCovered);
    R_Free(colSelected);
}

static double compute_reduced_and_lb(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    const double *t,
    double *ls
) {
    double ZLB_rows = 0.0;
    double ZLB_neg = 0.0;

    for (int i = 0; i < rows; ++i) {
        ZLB_rows += t[i];
    }

    for (int c = 0; c < cols; ++c) {
        double sum = 0.0;
        for (int k = 0; k < rowsCoveredCount[c]; ++k) {
            int r = rowsCovered[c][k];
            sum += t[r];
        }
        ls[c] = 1.0 - sum;
        if (ls[c] < 0.0) ZLB_neg += ls[c];
    }

    return ZLB_rows + ZLB_neg;
}

static int subgradient_update(
    int rows,
    int cols,
    int **colsCovering,
    int *colsCoveringCount,
    const double *ls,
    double UB,
    double ZLB,
    double *t,
    lagr_subgradient_state *state
) {
    if (!state || !state->prev_t) return 0;

    if (state->have_prev && ZLB + EPS < state->prev_zlb) {
        Memcpy(t, state->prev_t, rows);
        lagr_state_contract_phi(state);
        state->have_prev = 0;
        state->prev_step = 0.0;
        if (state->prev_direction) {
            memset(state->prev_direction, 0, (size_t)rows * sizeof(double));
        }
        return 1;
    }

    unsigned char *x = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    if (!x) return 0;

    for (int c = 0; c < cols; ++c) {
        x[c] = (ls[c] < 0.0) ? 1 : 0;
    }

    double *slack = (double*) R_Calloc((size_t)rows, double);
    if (!slack) {
        R_Free(x);
        return 0;
    }

    double sum_s2 = 0.0;

    for (int i = 0; i < rows; ++i) {
        int covered = 0;
        for (int k = 0; k < colsCoveringCount[i]; ++k) {
            int c = colsCovering[i][k];
            covered += x[c];
        }
        double s = 1.0 - (double)covered;
        slack[i] = s;
        sum_s2 += s * s;
    }

    R_Free(x);

    if (state->deflection_alpha > 0.0 && state->prev_direction) {
        sum_s2 = 0.0;
        for (int i = 0; i < rows; ++i) {
            slack[i] += state->deflection_alpha * state->prev_direction[i];
            sum_s2 += slack[i] * slack[i];
        }
    }

    if (sum_s2 <= 1e-15) {
        R_Free(slack);
        return 0;
    }

    double gap = UB - ZLB;
    if (gap <= 0.0) {
        R_Free(slack);
        return 0;
    }

    double step = state->phi * (gap / sum_s2);
    if (step <= 0.0) {
        R_Free(slack);
        return 0;
    }

    Memcpy(state->prev_t, t, rows);
    state->prev_zlb = ZLB;
    state->prev_step = step;
    state->have_prev = 1;

    for (int i = 0; i < rows; ++i) {
        double trial = t[i] + step * slack[i];
        if (trial < 0.0) trial = 0.0;

        /* Stabilize multipliers by damping the raw projected subgradient move. */
        double beta = state->stabilization_beta;
        double val = (1.0 - beta) * t[i] + beta * trial;
        t[i] = (val > 0.0) ? val : 0.0;
    }

    if (state->deflection_alpha > 0.0 && state->prev_direction) {
        Memcpy(state->prev_direction, slack, rows);
    }

    state->stagnation_iter++;
    if (state->stagnation_iter >= state->stagnation_period) {
        state->stagnation_iter = 0;
        if (state->have_best) {
            Memcpy(t, state->best_t, rows);
            state->have_prev = 0;
        }
        lagr_state_contract_phi(state);
    }

    R_Free(slack);
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

    double *u = (double*) R_Calloc((size_t)n, double);
    if (!u) {
        double uniform = 1.0 / (double)n;
        for (int i = 0; i < n; ++i) v[i] = uniform;
        return;
    }

    Memcpy(u, v, n);
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

    R_Free(u);
}

static int bundle_update(
    int rows,
    int cols,
    int **colsCovering,
    int *colsCoveringCount,
    const double *ls,
    double UB,
    double ZLB,
    double *t,
    LagrangianBundleState *bs
) {
    if (!bs || bs->memory <= 0 || !bs->center_t || !bs->cut_alpha || !bs->cut_slacks) {
        return 0;
    }

    unsigned char *x = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    if (!x) return 0;

    for (int c = 0; c < cols; ++c) {
        x[c] = (ls[c] < 0.0) ? 1 : 0;
    }

    double *slack = (double*) R_Calloc((size_t)rows, double);
    if (!slack) {
        R_Free(x);
        return 0;
    }

    double sum_s2 = 0.0;
    for (int r = 0; r < rows; ++r) {
        int covered = 0;
        for (int k = 0; k < colsCoveringCount[r]; ++k) {
            int c = colsCovering[r][k];
            covered += x[c];
        }
        slack[r] = 1.0 - (double)covered;
        sum_s2 += slack[r] * slack[r];
    }

    R_Free(x);

    if (sum_s2 <= 1e-15 || UB - ZLB <= 0.0) {
        R_Free(slack);
        return 0;
    }

    if (bs->center_value <= -DBL_MAX / 2.0) {
        Memcpy(bs->center_t, t, rows);
        if (bs->prev_center_t) {
            Memcpy(bs->prev_center_t, t, rows);
        }
        bs->center_value = ZLB;
        bs->has_trial = 0;
    } else if (bs->has_trial) {
        double actual_gain = ZLB - bs->center_value;
        if (actual_gain > bs->serious_tol) {
            if (bs->prev_center_t) {
                Memcpy(bs->prev_center_t, bs->center_t, rows);
            }
            Memcpy(bs->center_t, t, rows);
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
            Memcpy(bs->prev_center_t, bs->center_t, rows);
        }
        Memcpy(bs->center_t, t, rows);
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
        Memcpy(bs->cut_slacks + (size_t)slot * (size_t)rows, slack, rows);
        bs->cut_alpha[slot] = alpha;
        bs->count++;
    } else {
        int slot = bs->head;
        Memcpy(bs->cut_slacks + (size_t)slot * (size_t)rows, slack, rows);
        bs->cut_alpha[slot] = alpha;
        bs->head = (bs->head + 1) % bs->memory;
    }

    int m = bs->count;
    if (m <= 0 || (bs->prox_mu >= bs->max_mu - 1e-12 && bs->null_streak >= bs->memory * 8)) {
        R_Free(slack);
        return 0;
    }

    double *base_t = (double*) R_Calloc((size_t)rows, double);
    double *linear = (double*) R_Calloc((size_t)m, double);
    double *gram = (double*) R_Calloc((size_t)m * (size_t)m, double);
    double *lambda = (double*) R_Calloc((size_t)m, double);
    double *gradient = (double*) R_Calloc((size_t)m, double);
    double *gbar = (double*) R_Calloc((size_t)rows, double);
    if (!base_t || !linear || !gram || !lambda || !gradient || !gbar) {
        if (base_t) R_Free(base_t);
        if (linear) R_Free(linear);
        if (gram) R_Free(gram);
        if (lambda) R_Free(lambda);
        if (gradient) R_Free(gradient);
        if (gbar) R_Free(gbar);
        R_Free(slack);
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
        R_Free(base_t);
        R_Free(linear);
        R_Free(gram);
        R_Free(lambda);
        R_Free(gradient);
        R_Free(gbar);
        R_Free(slack);
        return 0;
    }

    bs->predicted_gain = predicted_gain;
    bs->has_trial = 1;

    R_Free(base_t);
    R_Free(linear);
    R_Free(gram);
    R_Free(lambda);
    R_Free(gradient);
    R_Free(gbar);
    R_Free(slack);
    return 1;
}

static void heuristic_row_min_rc(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int *sol,
    int *sol_len
) {
    (void)cols;

    bool *rowCovered = (bool*) R_Calloc((size_t)rows, bool);
    unsigned char *x = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    int out = 0, covered = 0;

    if (!rowCovered || !x) {
        if (rowCovered) R_Free(rowCovered);
        if (x) R_Free(x);
        *sol_len = -1;
        return;
    }

    for (int i = 0; i < rows; ++i) {
        if (rowCovered[i]) continue;

        int best = -1;
        double bestLS = DBL_MAX;
        int bestOnes = -1;
        double bestWeight = -DBL_MAX;

        for (int k = 0; k < colsCoveringCount[i]; ++k) {
            int c = colsCovering[i][k];
            double ls = lagr_score ? lagr_score[c] : 0.0;
            double w = weights ? weights[c] : 0.0;
            int ones = rowsCoveredCount[c];

            bool better = false;
            if (ls < bestLS) {
                better = true;
            }
            else if (fabs(ls - bestLS) <= EPS) {
                if (w > bestWeight) {
                    better = true;
                }
                else if (fabs(w - bestWeight) <= EPS) {
                    if (ones > bestOnes) {
                        better = true;
                    }
                    else if (ones == bestOnes) {
                        if (best == -1 || c < best) better = true;
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

        if (best == -1) {
            R_Free(rowCovered);
            R_Free(x);
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

    if (out > 0) {
        int *coverCount = (int*) R_Calloc((size_t)rows, int);
        int *order = (int*) R_Calloc((size_t)out, int);
        double *ls_sel = (double*) R_Calloc((size_t)out, double);

        if (coverCount && order && ls_sel) {
            for (int i = 0; i < out; ++i) {
                int s = sol[i];
                order[i] = i;
                ls_sel[i] = lagr_score ? lagr_score[s] : 0.0;
                for (int c = 0; c < rowsCoveredCount[s]; ++c) {
                    int r = rowsCovered[s][c];
                    coverCount[r]++;
                }
            }

            for (int i = 0; i < out; ++i) {
                for (int j = i + 1; j < out; ++j) {
                    if (ls_sel[order[j]] > ls_sel[order[i]]) {
                        int tmp = order[i];
                        order[i] = order[j];
                        order[j] = tmp;
                    }
                }
            }

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

                    for (int j = idx + 1; j < out; ++j) {
                        sol[j - 1] = sol[j];
                    }
                    out--;

                    for (int j = o + 1; j < *sol_len; ++j) {
                        if (order[j] > idx) order[j]--;
                    }
                }
            }
            *sol_len = out;
        }

        if (coverCount) R_Free(coverCount);
        if (order) R_Free(order);
        if (ls_sel) R_Free(ls_sel);
    }

    prune_redundancy(rows, rowsCovered, rowsCoveredCount, sol, sol_len);

    R_Free(rowCovered);
    R_Free(x);
}

static void heuristic_negative_rc_repair(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int *sol,
    int *sol_len
) {
    bool *rowCovered = (bool*) R_Calloc((size_t)rows, bool);
    unsigned char *selected = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    int out = 0;
    int covered = 0;

    if (!rowCovered || !selected) {
        if (rowCovered) R_Free(rowCovered);
        if (selected) R_Free(selected);
        *sol_len = -1;
        return;
    }

    for (int c = 0; c < cols; ++c) {
        if (lagr_score && lagr_score[c] < 0.0) {
            selected[c] = 1;
            sol[out++] = c;
            for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                int r = rowsCovered[c][k];
                if (!rowCovered[r]) {
                    rowCovered[r] = true;
                    covered++;
                }
            }
        }
    }

    while (covered < rows) {
        int r0 = -1;
        for (int r = 0; r < rows; ++r) {
            if (!rowCovered[r]) {
                r0 = r;
                break;
            }
        }
        if (r0 < 0) break;

        int best = -1;
        int bestNew = -1;
        double bestLS = DBL_MAX;
        double bestW = -DBL_MAX;

        for (int k = 0; k < colsCoveringCount[r0]; ++k) {
            int c = colsCovering[r0][k];
            if (selected[c]) continue;

            int newCover = 0;
            for (int j = 0; j < rowsCoveredCount[c]; ++j) {
                int rr = rowsCovered[c][j];
                if (!rowCovered[rr]) newCover++;
            }
            if (newCover <= 0) continue;

            double ls = lagr_score ? lagr_score[c] : 0.0;
            double w = weights ? weights[c] : 0.0;

            bool better = false;
            if (newCover > bestNew) {
                better = true;
            }
            else if (newCover == bestNew) {
                if (ls < bestLS) {
                    better = true;
                }
                else if (fabs(ls - bestLS) <= EPS) {
                    if (w > bestW) {
                        better = true;
                    }
                    else if (fabs(w - bestW) <= EPS) {
                        if (best == -1 || c < best) {
                            better = true;
                        }
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

        if (best == -1) {
            R_Free(rowCovered);
            R_Free(selected);
            *sol_len = -1;
            return;
        }

        selected[best] = 1;
        sol[out++] = best;
        for (int k = 0; k < rowsCoveredCount[best]; ++k) {
            int rr = rowsCovered[best][k];
            if (!rowCovered[rr]) {
                rowCovered[rr] = true;
                covered++;
            }
        }
    }

    *sol_len = out;
    prune_redundancy(rows, rowsCovered, rowsCoveredCount, sol, sol_len);

    R_Free(rowCovered);
    R_Free(selected);
}

static double solution_total_weight(const int *sol, int sol_len, const double *weights) {
    if (!sol || sol_len <= 0 || !weights) return 0.0;
    double s = 0.0;
    for (int i = 0; i < sol_len; ++i) {
        s += weights[sol[i]];
    }
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
        mean_support += (double) colsCoveringCount[i];
    }
    mean_support /= (double) rows;
    if (mean_support <= 0.0) {
        mean_support = 1.0;
    }

    for (int i = 0; i < rows; ++i) {
        double support = (double) colsCoveringCount[i];
        if (support <= 0.0) {
            support = 1.0;
        }
        t[i] = mean_support / support;
        if (t[i] < 0.1) t[i] = 0.1;
    }
}

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
);

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
);

static void polish_candidate_solution(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *ls,
    const double *weights,
    int *sol,
    int *sol_size,
    int passes
) {
    if (!sol || !sol_size || *sol_size <= 0 || passes <= 0) return;

    local_search_drop1_repair(
        rows, cols, rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        ls, weights, sol, sol_size, passes
    );
    local_search_drop2_repair(
        rows, cols, rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        ls, weights, sol, sol_size, passes
    );
}

static double solution_proxy_value(
    const int *sol,
    int sol_len,
    const double *weights,
    const int *rowsCoveredCount
) {
    if (!sol || sol_len <= 0) return 0.0;
    double v = 0.0;
    if (weights) {
        for (int i = 0; i < sol_len; ++i) {
            v += weights[sol[i]];
        }
    }
    else if (rowsCoveredCount) {
        for (int i = 0; i < sol_len; ++i) {
            v += (double)rowsCoveredCount[sol[i]];
        }
    }
    return v;
}

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
    if (!sol || !sol_len || *sol_len <= 1 || max_passes <= 0) return;

    const double PROXY_EPS = 1e-12;
    double curr_proxy = solution_proxy_value(sol, *sol_len, weights, rowsCoveredCount);

    bool *rowCovered = (bool*) R_Calloc((size_t)rows, bool);
    unsigned char *colSelected = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    int *cand = (int*) R_Calloc((size_t)cols, int);
    if (!rowCovered || !colSelected || !cand) {
        if (rowCovered) R_Free(rowCovered);
        if (colSelected) R_Free(colSelected);
        if (cand) R_Free(cand);
        return;
    }

    for (int pass = 0; pass < max_passes; ++pass) {
        int improved = 0;

        for (int drop_idx = 0; drop_idx < *sol_len; ++drop_idx) {
            int cand_len = 0;
            Memzero(colSelected, cols);

            for (int i = 0; i < *sol_len; ++i) {
                if (i == drop_idx) continue;
                int c = sol[i];
                cand[cand_len++] = c;
                colSelected[c] = 1;
            }

            Memzero(rowCovered, rows);
            int covered = 0;
            for (int i = 0; i < cand_len; ++i) {
                int c = cand[i];
                for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                    int r = rowsCovered[c][k];
                    if (!rowCovered[r]) {
                        rowCovered[r] = true;
                        covered++;
                    }
                }
            }

            while (covered < rows) {
                maybe_check_user_interrupt(covered);
                int r0 = -1;
                for (int r = 0; r < rows; ++r) {
                    if (!rowCovered[r]) {
                        r0 = r;
                        break;
                    }
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
                    double w = weights ? weights[c] : 0.0;

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
                        best = c;
                        bestNew = newCover;
                        bestLS = ls;
                        bestW = w;
                    }
                }

                if (best < 0) {
                    cand_len = -1;
                    break;
                }

                colSelected[best] = 1;
                cand[cand_len++] = best;
                for (int k = 0; k < rowsCoveredCount[best]; ++k) {
                    int rr = rowsCovered[best][k];
                    if (!rowCovered[rr]) {
                        rowCovered[rr] = true;
                        covered++;
                    }
                }
            }

            if (cand_len <= 0 || cand_len > *sol_len) continue;

            prune_redundancy(rows, rowsCovered, rowsCoveredCount, cand, &cand_len);

            if (cand_len > 0 && cand_len <= *sol_len) {
                double cand_proxy = solution_proxy_value(cand, cand_len, weights, rowsCoveredCount);
                if (cand_proxy > curr_proxy + PROXY_EPS) {
                    Memcpy(sol, cand, cand_len);
                    *sol_len = cand_len;
                    curr_proxy = cand_proxy;
                    improved = 1;
                    break;
                }
            }
        }

        if (!improved) break;
    }

    R_Free(rowCovered);
    R_Free(colSelected);
    R_Free(cand);
}

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
    if (!sol || !sol_len || *sol_len <= 2 || max_passes <= 0) return;

    bool *rowCovered = (bool*) R_Calloc((size_t)rows, bool);
    unsigned char *colSelected = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    int *cand = (int*) R_Calloc((size_t)cols, int);
    if (!rowCovered || !colSelected || !cand) {
        if (rowCovered) R_Free(rowCovered);
        if (colSelected) R_Free(colSelected);
        if (cand) R_Free(cand);
        return;
    }

    for (int pass = 0; pass < max_passes; ++pass) {
        int improved = 0;

        for (int drop_i = 0; drop_i < *sol_len; ++drop_i) {
            for (int drop_j = drop_i + 1; drop_j < *sol_len; ++drop_j) {
                int cand_len = 0;
                Memzero(colSelected, cols);

                for (int i = 0; i < *sol_len; ++i) {
                    if (i == drop_i || i == drop_j) continue;
                    int c = sol[i];
                    cand[cand_len++] = c;
                    colSelected[c] = 1;
                }

                Memzero(rowCovered, rows);
                int covered = 0;
                for (int i = 0; i < cand_len; ++i) {
                    int c = cand[i];
                    for (int k = 0; k < rowsCoveredCount[c]; ++k) {
                        int r = rowsCovered[c][k];
                        if (!rowCovered[r]) {
                            rowCovered[r] = true;
                            covered++;
                        }
                    }
                }

                while (covered < rows) {
                    maybe_check_user_interrupt(covered);
                    int r0 = -1;
                    for (int r = 0; r < rows; ++r) {
                        if (!rowCovered[r]) {
                            r0 = r;
                            break;
                        }
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
                        double w = weights ? weights[c] : 0.0;

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
                            best = c;
                            bestNew = newCover;
                            bestLS = ls;
                            bestW = w;
                        }
                    }

                    if (best < 0) {
                        cand_len = -1;
                        break;
                    }

                    colSelected[best] = 1;
                    cand[cand_len++] = best;
                    for (int k = 0; k < rowsCoveredCount[best]; ++k) {
                        int rr = rowsCovered[best][k];
                        if (!rowCovered[rr]) {
                            rowCovered[rr] = true;
                            covered++;
                        }
                    }
                }

                if (cand_len <= 0 || cand_len >= *sol_len) continue;

                prune_redundancy(rows, rowsCovered, rowsCoveredCount, cand, &cand_len);

                if (cand_len > 0 && cand_len < *sol_len) {
                    Memcpy(sol, cand, cand_len);
                    *sol_len = cand_len;
                    improved = 1;
                    break;
                }
            }

            if (improved) break;
        }

        if (!improved) break;
    }

    R_Free(rowCovered);
    R_Free(colSelected);
    R_Free(cand);
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
    int target;
    long node_limit;
    long nodes;
    int *coverCount;
    unsigned char *colSelected;
    int *chosen;
    int *out;
} PolishSearch;

static bool polish_candidate_better(const PolishCandidate *a, const PolishCandidate *b) {
    if (a->new_cover != b->new_cover) return a->new_cover > b->new_cover;
    if (fabs(a->lagr_score - b->lagr_score) > EPS) return a->lagr_score < b->lagr_score;
    if (a->cover_count != b->cover_count) return a->cover_count > b->cover_count;
    if (fabs(a->weight - b->weight) > EPS) return a->weight > b->weight;
    return a->col < b->col;
}

static void polish_sort_candidates(PolishCandidate *candidates, int count) {
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
    int best_row = -1;
    int best_count = INT_MAX;
    int max_new_cover = 0;

    for (int row = 0; row < search->rows; ++row) {
        if (search->coverCount[row] > 0) continue;

        int viable = 0;
        for (int i = 0; i < search->colsCoveringCount[row]; ++i) {
            int col = search->colsCovering[row][i];
            if (search->colSelected[col]) continue;
            int new_cover = polish_new_cover(search, col);
            if (new_cover <= 0) continue;
            ++viable;
            if (new_cover > max_new_cover) max_new_cover = new_cover;
        }

        if (viable == 0) return -1;
        if (viable < best_count) {
            best_count = viable;
            best_row = row;
        }
    }

    if (max_new_cover <= 0) return -1;
    if (uncovered > slots_left * max_new_cover) return -1;
    return best_row;
}

static bool polish_search_rec(
    PolishSearch *search,
    int depth,
    int covered
) {
    if (covered >= search->rows) {
        Memcpy(search->out, search->chosen, depth);
        return true;
    }

    if (depth >= search->target) return false;
    if (search->nodes++ >= search->node_limit) return false;

    int slots_left = search->target - depth;
    int uncovered = search->rows - covered;
    int row = polish_choose_row(search, uncovered, slots_left);
    if (row < 0) return false;

    int raw_count = search->colsCoveringCount[row];
    PolishCandidate *candidates = (PolishCandidate*) R_Calloc((size_t)raw_count, PolishCandidate);
    if (!candidates) return false;

    int candidate_count = 0;
    for (int i = 0; i < raw_count; ++i) {
        int col = search->colsCovering[row][i];
        if (search->colSelected[col]) continue;

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

        if (polish_search_rec(search, depth + 1, covered + added)) {
            R_Free(candidates);
            return true;
        }

        for (int j = 0; j < search->rowsCoveredCount[col]; ++j) {
            int covered_row = search->rowsCovered[col][j];
            --search->coverCount[covered_row];
        }

        search->colSelected[col] = 0;
    }

    R_Free(candidates);
    return false;
}

static int polish_find_smaller_cover(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double *lagr_score,
    const double *weights,
    int target,
    long node_limit,
    int *out
) {
    if (target <= 0 || node_limit <= 0) return 0;

    int *coverCount = (int*) R_Calloc((size_t)rows, int);
    unsigned char *colSelected = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    int *chosen = (int*) R_Calloc((size_t)target, int);

    if (!coverCount || !colSelected || !chosen) {
        if (coverCount) R_Free(coverCount);
        if (colSelected) R_Free(colSelected);
        if (chosen) R_Free(chosen);
        return 0;
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
        .target = target,
        .node_limit = node_limit,
        .nodes = 0,
        .coverCount = coverCount,
        .colSelected = colSelected,
        .chosen = chosen,
        .out = out
    };

    int found = polish_search_rec(&search, 0, 0) ? 1 : 0;

    R_Free(coverCount);
    R_Free(colSelected);
    R_Free(chosen);

    return found;
}

typedef struct {
    int max_iter;
    int heur_every;
    int local_passes;
    int small_gap_threshold;
    int small_gap_extra_passes;
    double step_coef;
    double step_min;
    int halve_period;
    double phi_contract;
    double stabilization_beta;
    int rarity_init;
    double deflection_alpha;
    int polish_enabled;
    long polish_node_limit;
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

static int lagr_max_int(int a, int b) {
    return a > b ? a : b;
}

static void lagr_config_for_best_bound(LagrangianConfig *cfg) {
    cfg->max_iter = 20000;
    cfg->heur_every = 1;
    cfg->local_passes = 10;
    cfg->small_gap_threshold = 0;
    cfg->small_gap_extra_passes = cfg->local_passes;
    cfg->step_coef = 2.0;
    cfg->step_min = 0.005;
    cfg->halve_period = 8;
    cfg->phi_contract = 0.95;
    cfg->stabilization_beta = 1.0;
    cfg->rarity_init = 0;
    cfg->deflection_alpha = 0.0;
    cfg->polish_enabled = 1;
    cfg->polish_node_limit = 5000;
    cfg->portfolio_enabled = 1;
    cfg->portfolio_max_profiles = 6;
    cfg->portfolio_deflection_alpha = 0.75;
    cfg->hybrid_bundle_portfolio = 1;
    cfg->bundle_enabled = 0;
    cfg->bundle_interval = 8;
    cfg->bundle_memory = 8;
    cfg->bundle_prox_mu = 8.0;
    cfg->bundle_min_mu = 1.0;
    cfg->bundle_max_mu = 1024.0;
    cfg->bundle_null_expand = 2.0;
    cfg->bundle_serious_shrink = 0.8;
    cfg->bundle_serious_fraction = 0.05;
    cfg->bundle_serious_tol = 1e-9;
    cfg->bundle_momentum = 0.35;
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
        profiles[4].phi_contract = 0.99;
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

static void solvePIchart_lagrangian_core(
    int rows,
    int cols,
    int **rowsCovered,
    int *rowsCoveredCount,
    int **colsCovering,
    int *colsCoveringCount,
    const double weights[],
    const LagrangianConfig *cfg,
    int *solution,
    int *solmin,
    double *best_lb_out,
    double *lagr_score_out
) {
    *solmin = -1;
    if (rows <= 0 || cols <= 0 || !rowsCovered || !rowsCoveredCount || !colsCovering || !colsCoveringCount || !cfg || !solution || !solmin) return;

    int local_passes = cfg->local_passes;
    int small_gap_threshold = cfg->small_gap_threshold;
    int small_gap_extra_passes = cfg->small_gap_extra_passes;

    double *t = (double*) R_Calloc((size_t)rows, double);
    double *ls = (double*) R_Calloc((size_t)cols, double);
    int *sol_tmp = (int*) R_Calloc((size_t)cols, int);
    int *sol_tmp2 = (int*) R_Calloc((size_t)cols, int);
    int *sol_tmp3 = (int*) R_Calloc((size_t)cols, int);
    int best_sol_size = -1;
    lagr_subgradient_state sg_state = (lagr_subgradient_state){0};
    LagrangianBundleState bundle_state = {0};

    if (!t || !ls || !sol_tmp || !sol_tmp2 || !sol_tmp3 || !lagr_state_init(&sg_state, rows)) {
        if (t) R_Free(t);
        if (ls) R_Free(ls);
        if (sol_tmp) R_Free(sol_tmp);
        if (sol_tmp2) R_Free(sol_tmp2);
        if (sol_tmp3) R_Free(sol_tmp3);
        lagr_state_free(&sg_state);
        *solmin = -1;
        if (best_lb_out) *best_lb_out = -DBL_MAX;
        if (lagr_score_out) {
            for (int c = 0; c < cols; ++c) lagr_score_out[c] = 0.0;
        }
        return;
    }

    int max_iter = cfg->max_iter;
    int heur_every = cfg->heur_every;
    sg_state.phi = cfg->step_coef;
    sg_state.phi_min = cfg->step_min;
    sg_state.stagnation_period = cfg->halve_period;
    sg_state.phi_contract = cfg->phi_contract;
    sg_state.stabilization_beta = cfg->stabilization_beta;
    sg_state.deflection_alpha = cfg->deflection_alpha;

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
        bundle_state.center_t = (double*) R_Calloc((size_t)rows, double);
        bundle_state.prev_center_t = (double*) R_Calloc((size_t)rows, double);
        bundle_state.cut_alpha = (double*) R_Calloc((size_t)memory, double);
        bundle_state.cut_slacks = (double*) R_Calloc((size_t)memory * (size_t)rows, double);
        if (!bundle_state.center_t || !bundle_state.prev_center_t ||
            !bundle_state.cut_alpha || !bundle_state.cut_slacks) {
            if (bundle_state.center_t) R_Free(bundle_state.center_t);
            if (bundle_state.prev_center_t) R_Free(bundle_state.prev_center_t);
            if (bundle_state.cut_alpha) R_Free(bundle_state.cut_alpha);
            if (bundle_state.cut_slacks) R_Free(bundle_state.cut_slacks);
            R_Free(t);
            R_Free(ls);
            R_Free(sol_tmp);
            R_Free(sol_tmp2);
            R_Free(sol_tmp3);
            lagr_state_free(&sg_state);
            *solmin = -1;
            if (best_lb_out) *best_lb_out = -DBL_MAX;
            if (lagr_score_out) {
                for (int c = 0; c < cols; ++c) lagr_score_out[c] = 0.0;
            }
            return;
        }
    }

    lagr_initialize_multipliers(rows, colsCoveringCount, cfg->rarity_init, t);

    compute_reduced_and_lb(rows, cols, rowsCovered, rowsCoveredCount, t, ls);

    int sol_size = -1, sol_size2 = -1, sol_size3 = -1;

    heuristic_row_min_rc(
        rows, cols, rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        ls, weights, sol_tmp, &sol_size
    );

    greedy_from_lagr_scores(
        rows, cols, rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        ls, weights, sol_tmp2, &sol_size2
    );
    heuristic_negative_rc_repair(
        rows, cols, rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        ls, weights, sol_tmp3, &sol_size3
    );

    if (local_passes > 0) {
        polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                  colsCovering, colsCoveringCount,
                                  ls, weights, sol_tmp, &sol_size, local_passes);
        polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                  colsCovering, colsCoveringCount,
                                  ls, weights, sol_tmp2, &sol_size2, local_passes);
        polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                  colsCovering, colsCoveringCount,
                                  ls, weights, sol_tmp3, &sol_size3, local_passes);
    }

    if (sol_size == -1 && sol_size2 == -1 && sol_size3 == -1) {
        R_Free(t);
        R_Free(ls);
        R_Free(sol_tmp);
        R_Free(sol_tmp2);
        R_Free(sol_tmp3);
        lagr_state_free(&sg_state);
        if (bundle_state.center_t) R_Free(bundle_state.center_t);
        if (bundle_state.prev_center_t) R_Free(bundle_state.prev_center_t);
        if (bundle_state.cut_alpha) R_Free(bundle_state.cut_alpha);
        if (bundle_state.cut_slacks) R_Free(bundle_state.cut_slacks);
        *solmin = -1;
        if (best_lb_out) *best_lb_out = -DBL_MAX;
        if (lagr_score_out) {
            for (int c = 0; c < cols; ++c) lagr_score_out[c] = 0.0;
        }
        return;
    }

    double w1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
    double w2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
    double w3 = (sol_size3 != -1) ? solution_total_weight(sol_tmp3, sol_size3, weights) : -DBL_MAX;

    double bestUB;
    int bestTerms;
    double bestWeight;

    if (
        sol_size != -1 &&
        (sol_size2 == -1 || sol_size < sol_size2 || (sol_size == sol_size2 && w1 >= w2)) &&
        (sol_size3 == -1 || sol_size < sol_size3 || (sol_size == sol_size3 && w1 >= w3))
    ) {
        bestUB = (double)sol_size;
        bestTerms = sol_size;
        bestWeight = w1;
        best_sol_size = sol_size;
        Memcpy(solution, sol_tmp, sol_size);
    }
    else if (
        sol_size2 != -1 &&
        (sol_size3 == -1 || sol_size2 < sol_size3 || (sol_size2 == sol_size3 && w2 >= w3))
    ) {
        bestUB = (double)sol_size2;
        bestTerms = sol_size2;
        bestWeight = w2;
        best_sol_size = sol_size2;
        Memcpy(solution, sol_tmp2, sol_size2);
    }
    else {
        bestUB = (double)sol_size3;
        bestTerms = sol_size3;
        bestWeight = w3;
        best_sol_size = sol_size3;
        Memcpy(solution, sol_tmp3, sol_size3);
    }

    double bestLB = -DBL_MAX;

    for (int it = 0; it < max_iter; ++it) {
        maybe_check_user_interrupt(it);
        double ZLB = compute_reduced_and_lb(rows, cols, rowsCovered, rowsCoveredCount, t, ls);
        double LBint = ceil(ZLB - 1e-12);

        if (LBint > bestLB + EPS) {
            bestLB = LBint;
            lagr_state_note_best_dual(&sg_state, t, rows, ZLB);
            lagr_state_reset_progress(&sg_state);
        }

        if (it % heur_every == 0) {
            heuristic_row_min_rc(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp, &sol_size
            );
            greedy_from_lagr_scores(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp2, &sol_size2
            );
            heuristic_negative_rc_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp3, &sol_size3
            );

            if (local_passes > 0) {
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp, &sol_size, local_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp2, &sol_size2, local_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp3, &sol_size3, local_passes);
            }

            if (small_gap_extra_passes > local_passes &&
                small_gap_threshold > 0 &&
                bestTerms < INT_MAX &&
                bestTerms - (int) LBint <= small_gap_threshold) {
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp, &sol_size, small_gap_extra_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp2, &sol_size2, small_gap_extra_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp3, &sol_size3, small_gap_extra_passes);
            }

            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
            double candW3 = (sol_size3 != -1) ? solution_total_weight(sol_tmp3, sol_size3, weights) : -DBL_MAX;

            if (
                sol_size != -1 &&
                (sol_size < bestTerms || (sol_size == bestTerms && candW1 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size;
                bestTerms = sol_size;
                bestWeight = candW1;
                best_sol_size = sol_size;
                Memcpy(solution, sol_tmp, sol_size);
                lagr_state_reset_progress(&sg_state);
            }

            if (
                sol_size2 != -1 &&
                (sol_size2 < bestTerms || (sol_size2 == bestTerms && candW2 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size2;
                bestTerms = sol_size2;
                bestWeight = candW2;
                best_sol_size = sol_size2;
                Memcpy(solution, sol_tmp2, sol_size2);
                lagr_state_reset_progress(&sg_state);
            }

            if (
                sol_size3 != -1 &&
                (sol_size3 < bestTerms || (sol_size3 == bestTerms && candW3 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size3;
                bestTerms = sol_size3;
                bestWeight = candW3;
                best_sol_size = sol_size3;
                Memcpy(solution, sol_tmp3, sol_size3);
                lagr_state_reset_progress(&sg_state);
            }
        }

        if (bestUB <= LBint + EPS) break;

        int used_subgradient_update = 0;
        int updated = 0;

        int use_bundle_update = cfg->bundle_enabled &&
            (cfg->bundle_interval <= 1 || it % cfg->bundle_interval == 0);

        if (use_bundle_update) {
            updated = bundle_update(
                rows,
                cols,
                colsCovering,
                colsCoveringCount,
                ls,
                bestUB,
                ZLB,
                t,
                &bundle_state
            );
        }

        if (!use_bundle_update || !updated) {
            updated = subgradient_update(
                rows, cols, colsCovering, colsCoveringCount,
                ls, bestUB, ZLB, t, &sg_state
            );
            used_subgradient_update = 1;
        }

        if (!updated || (used_subgradient_update && sg_state.phi <= sg_state.phi_min + EPS)) {
            heuristic_row_min_rc(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp, &sol_size
            );
            greedy_from_lagr_scores(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp2, &sol_size2
            );
            heuristic_negative_rc_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp3, &sol_size3
            );

            if (local_passes > 0) {
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp, &sol_size, local_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp2, &sol_size2, local_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp3, &sol_size3, local_passes);
            }

            if (small_gap_extra_passes > local_passes &&
                small_gap_threshold > 0 &&
                bestTerms < INT_MAX &&
                bestTerms - (int) ceil(bestLB - 1e-12) <= small_gap_threshold) {
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp, &sol_size, small_gap_extra_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp2, &sol_size2, small_gap_extra_passes);
                polish_candidate_solution(rows, cols, rowsCovered, rowsCoveredCount,
                                          colsCovering, colsCoveringCount,
                                          ls, weights, sol_tmp3, &sol_size3, small_gap_extra_passes);
            }

            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
            double candW3 = (sol_size3 != -1) ? solution_total_weight(sol_tmp3, sol_size3, weights) : -DBL_MAX;

            if (
                sol_size != -1 &&
                (sol_size < bestTerms || (sol_size == bestTerms && candW1 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size;
                bestTerms = sol_size;
                bestWeight = candW1;
                best_sol_size = sol_size;
                Memcpy(solution, sol_tmp, sol_size);
            }

            if (
                sol_size2 != -1 &&
                (sol_size2 < bestTerms || (sol_size2 == bestTerms && candW2 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size2;
                bestTerms = sol_size2;
                bestWeight = candW2;
                best_sol_size = sol_size2;
                Memcpy(solution, sol_tmp2, sol_size2);
            }

            if (
                sol_size3 != -1 &&
                (sol_size3 < bestTerms || (sol_size3 == bestTerms && candW3 > bestWeight + EPS))
            ) {
                bestUB = (double)sol_size3;
                bestTerms = sol_size3;
                bestWeight = candW3;
                best_sol_size = sol_size3;
                Memcpy(solution, sol_tmp3, sol_size3);
            }
            break;
        }
    }

    if (cfg->polish_enabled && best_sol_size > 1) {
        int polish_target = best_sol_size - 1;
        if (
            bestLB >= (double)polish_target - EPS &&
            bestLB <= (double)polish_target + EPS &&
            polish_find_smaller_cover(
                rows, cols,
                rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights,
                polish_target,
                cfg->polish_node_limit,
                sol_tmp
            )
        ) {
            best_sol_size = polish_target;
            bestUB = (double)polish_target;
            bestTerms = polish_target;
            bestWeight = solution_total_weight(sol_tmp, polish_target, weights);
            Memcpy(solution, sol_tmp, polish_target);
        }
    }

    (void)bestUB;
    (void)bestTerms;
    (void)bestWeight;

    *solmin = best_sol_size;
    if (best_lb_out) *best_lb_out = bestLB;
    if (lagr_score_out) {
        Memcpy(lagr_score_out, ls, cols);
    }

    R_Free(t);
    R_Free(ls);
    R_Free(sol_tmp);
    R_Free(sol_tmp2);
    R_Free(sol_tmp3);
    lagr_state_free(&sg_state);
    if (bundle_state.center_t) R_Free(bundle_state.center_t);
    if (bundle_state.prev_center_t) R_Free(bundle_state.prev_center_t);
    if (bundle_state.cut_alpha) R_Free(bundle_state.cut_alpha);
    if (bundle_state.cut_slacks) R_Free(bundle_state.cut_slacks);
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
    if (solmin) *solmin = -1;
    if (!pichart || foundPI <= 0 || ON_minterms <= 0 || !solution || !solmin) return;

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

    if (rc_ad == -2 || rc_ad == -1) {
        *solmin = -1;
        if (best_lb_out) *best_lb_out = -DBL_MAX;
        if (lagr_score_out) {
            for (int c = 0; c < foundPI; ++c) lagr_score_out[c] = 0.0;
        }
        return;
    }

    const int rows = ON_minterms;
    const int cols = foundPI;
    LagrangianConfig baseline;
    lagr_config_for_best_bound(&baseline);

    int portfolio_enabled = baseline.portfolio_enabled;
    double portfolio_deflection_alpha = baseline.portfolio_deflection_alpha;
    if (portfolio_enabled) {
        baseline.deflection_alpha = 0.0;
    }

    double best_report_lb = -DBL_MAX;
    solvePIchart_lagrangian_core(
        rows, cols,
        rowsCovered, rowsCoveredCount,
        colsCovering, colsCoveringCount,
        weights,
        &baseline,
        solution,
        solmin,
        &best_report_lb,
        lagr_score_out
    );

    int best_solmin = *solmin;
    if (best_lb_out) *best_lb_out = best_report_lb;

    if (portfolio_enabled) {
        LagrangianConfig profiles[LAGR_PORTFOLIO_PROFILE_COUNT];
        int profile_count = lagr_build_portfolio_configs(&baseline, portfolio_deflection_alpha, profiles);
        int max_profiles = baseline.portfolio_max_profiles;
        if (max_profiles < 0) max_profiles = 0;
        if (max_profiles > profile_count) max_profiles = profile_count;

        int *candidate_solution = (int*) R_Calloc((size_t)cols, int);
        double *candidate_scores = lagr_score_out ? (double*) R_Calloc((size_t)cols, double) : NULL;

        if (candidate_solution != NULL && (!lagr_score_out || candidate_scores != NULL)) {
            for (int profile = 1; profile < max_profiles; ++profile) {
                int candidate_solmin = -1;
                double candidate_lb = -DBL_MAX;
                memset(candidate_solution, 0, (size_t)cols * sizeof(int));
                if (candidate_scores) {
                    memset(candidate_scores, 0, (size_t)cols * sizeof(double));
                }

                solvePIchart_lagrangian_core(
                    rows, cols,
                    rowsCovered, rowsCoveredCount,
                    colsCovering, colsCoveringCount,
                    weights,
                    &profiles[profile],
                    candidate_solution,
                    &candidate_solmin,
                    &candidate_lb,
                    candidate_scores
                );

                if (candidate_solmin > 0 && (best_solmin <= 0 || candidate_solmin < best_solmin)) {
                    best_solmin = candidate_solmin;
                    *solmin = candidate_solmin;
                    Memcpy(solution, candidate_solution, candidate_solmin);
                }

                if (candidate_lb > best_report_lb + EPS) {
                    best_report_lb = candidate_lb;
                    if (best_lb_out) *best_lb_out = best_report_lb;
                    if (lagr_score_out && candidate_scores) {
                        Memcpy(lagr_score_out, candidate_scores, cols);
                    }
                }
            }
        }

        if (candidate_solution) R_Free(candidate_solution);
        if (candidate_scores) R_Free(candidate_scores);
    }

    free_adjacency(rowsCovered, rowsCoveredCount, cols, colsCovering, colsCoveringCount, rows);
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

    int *pichart = (int*) R_Calloc((size_t)rows * (size_t)cols, int);
    int *indices = (int*) R_Calloc((size_t)cols, int);
    if (!pichart || !indices) {
        if (pichart) R_Free(pichart);
        if (indices) R_Free(indices);
        error("Memory allocation failed in C_findminLagrangian().");
    }

    switch(TYPEOF(chart)) {
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
    solvePIchart_lagrangian(pichart, cols, rows, NULL, indices, &solmin, NULL, NULL);

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

    int *pichart = (int*) R_Calloc((size_t)rows * (size_t)cols, int);
    int *indices = (int*) R_Calloc((size_t)cols, int);
    if (!pichart || !indices) {
        if (pichart) R_Free(pichart);
        if (indices) R_Free(indices);
        error("Memory allocation failed in C_findminLagrangianInfo().");
    }

    switch(TYPEOF(chart)) {
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
    solvePIchart_lagrangian(pichart, cols, rows, NULL, indices, &solmin, &best_lb, NULL);

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
    SET_VECTOR_ELT(out, 2, ScalarReal(best_lb > -DBL_MAX / 2 ? best_lb : NA_REAL));
    setAttrib(out, R_NamesSymbol, names);

    R_Free(pichart);
    R_Free(indices);
    UNPROTECT(3);
    return out;
}
