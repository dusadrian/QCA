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

#include "solvePIchart_lagrangian.h"
#include <R.h>
#include <R_ext/RS.h>
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <stdbool.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#ifdef _OPENMP
    #undef match
    #include <omp.h>
#endif
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
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < cols; ++c) {
        int cnt = 0;
        for (int r = 0; r < rows; ++r) {
            cnt += pichart[c * rows + r];
        }
        rcc[c] = cnt;
    }
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
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
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < cols; ++c) {
        int k = 0;
        for (int r = 0; r < rows; ++r) {
            if (pichart[c * rows + r] != 0) {
                rc[c][k++] = r;
            }
        }
    }
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
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
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:ZLB_rows) schedule(static)
    #endif
    for (int i = 0; i < rows; ++i) {
        ZLB_rows += t[i];
    }
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:ZLB_neg) schedule(static)
    #endif
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
    double *step_coef,
    double step_min,
    int *stuck_iter,
    int stuck_halve_period
) {
    unsigned char *x = (unsigned char*) R_Calloc((size_t)cols, unsigned char);
    if (!x) return 0;
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
    for (int c = 0; c < cols; ++c) {
        x[c] = (ls[c] < 0.0) ? 1 : 0;
    }
    double *slack = (double*) R_Calloc((size_t)rows, double);
    if (!slack) {
        R_Free(x);
        return 0;
    }
    double sum_s2 = 0.0;
    #ifdef _OPENMP
        #pragma omp parallel for reduction(+:sum_s2) schedule(static)
    #endif
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
    if (sum_s2 <= 1e-15) {
        R_Free(slack);
        return 0;
    }
    double gap = UB - ZLB;
    if (gap <= 0.0) {
        R_Free(slack);
        return 0;
    }
    double step = (*step_coef) * (gap / sum_s2);
    if (step <= 0.0) {
        R_Free(slack);
        return 0;
    }
    #ifdef _OPENMP
        #pragma omp parallel for schedule(static)
    #endif
    for (int i = 0; i < rows; ++i) {
        double val = t[i] + step * slack[i];
        t[i] = (val > 0.0) ? val : 0.0;
    }
    (*stuck_iter)++;
    if (*stuck_iter >= stuck_halve_period) {
        *stuck_iter = 0;
        *step_coef *= 0.5;
        if (*step_coef < step_min) *step_coef = step_min;
    }
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
static double solution_total_weight(const int *sol, int sol_len, const double *weights) {
    if (!sol || sol_len <= 0 || !weights) return 0.0;
    double s = 0.0;
    for (int i = 0; i < sol_len; ++i) {
        s += weights[sol[i]];
    }
    return s;
}
static int lagr_local_passes_from_env(void) {
    static int cached = -1;
    if (cached >= 0) return cached;
    int def = 10;
    const char *env = getenv("CCUBES_LAGR_LOCAL_PASSES");
    if (env) {
        long v = strtol(env, NULL, 10);
        if (v >= 0 && v < 100000000) {
            cached = (int)v;
            return cached;
        }
    }
    cached = def;
    return cached;
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
void solvePIchart_lagrangian(
    int pichart[],
    const int foundPI,
    const int ON_minterms,
    const double weights[],
    int *solution,
    int *solmin
) {
    *solmin = -1;
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
        return;
    }
    int rows = ON_minterms, cols = foundPI;
    int local_passes = lagr_local_passes_from_env();
    double *t = (double*) R_Calloc((size_t)rows, double);
    double *ls = (double*) R_Calloc((size_t)cols, double);
    int *sol_tmp = (int*) R_Calloc((size_t)cols, int);
    int *sol_tmp2 = (int*) R_Calloc((size_t)cols, int);
    int best_sol_size = -1;
    int stuck_iter = 0;
    if (!t || !ls || !sol_tmp || !sol_tmp2) {
        if (t) R_Free(t);
        if (ls) R_Free(ls);
        if (sol_tmp) R_Free(sol_tmp);
        if (sol_tmp2) R_Free(sol_tmp2);
        free_adjacency(rowsCovered, rowsCoveredCount, cols, colsCovering, colsCoveringCount, rows);
        *solmin = -1;
        return;
    }
    int max_iter = 20000;
    int heur_every = 1;
    double step_coef = 2.0;
    double step_min = 0.005;
    int halve_period = 8;
    const char *s_env = NULL;
    s_env = getenv("CCUBES_LAGR_MAX_ITER");
    if (s_env) {
        long v = strtol(s_env, NULL, 10);
        if (v > 10 && v < 50000000) max_iter = (int)v;
    }
    s_env = getenv("CCUBES_LAGR_HEUR_EVERY");
    if (s_env) {
        long v = strtol(s_env, NULL, 10);
        if (v >= 1 && v < 10000) heur_every = (int)v;
    }
    s_env = getenv("CCUBES_LAGR_STEP_COEF");
    if (s_env) {
        double v = strtod(s_env, NULL);
        if (v > 0.0 && v < 1000.0) step_coef = v;
    }
    s_env = getenv("CCUBES_LAGR_STEP_MIN");
    if (s_env) {
        double v = strtod(s_env, NULL);
        if (v > 0.0 && v < step_coef) step_min = v;
    }
    s_env = getenv("CCUBES_LAGR_HALVE_PERIOD");
    if (s_env) {
        long v = strtol(s_env, NULL, 10);
        if (v >= 1 && v < 10000) halve_period = (int)v;
    }
    for (int i = 0; i < rows; ++i) {
        t[i] = 1.0;
    }
    compute_reduced_and_lb(rows, cols, rowsCovered, rowsCoveredCount, t, ls);
    int sol_size = -1, sol_size2 = -1;
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
    if (local_passes > 0) {
        if (sol_size > 0) {
            local_search_drop1_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp, &sol_size, local_passes
            );
            local_search_drop2_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp, &sol_size, local_passes
            );
        }
        if (sol_size2 > 0) {
            local_search_drop1_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp2, &sol_size2, local_passes
            );
            local_search_drop2_repair(
                rows, cols, rowsCovered, rowsCoveredCount,
                colsCovering, colsCoveringCount,
                ls, weights, sol_tmp2, &sol_size2, local_passes
            );
        }
    }
    if (sol_size == -1 && sol_size2 == -1) {
        R_Free(t);
        R_Free(ls);
        R_Free(sol_tmp);
        R_Free(sol_tmp2);
        free_adjacency(rowsCovered, rowsCoveredCount, cols, colsCovering, colsCoveringCount, rows);
        *solmin = -1;
        return;
    }
    double w1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
    double w2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
    double bestUB;
    int bestTerms;
    double bestWeight;
    if (
        sol_size != -1 &&
        (sol_size2 == -1 || sol_size < sol_size2 || (sol_size == sol_size2 && w1 >= w2))
    ) {
        bestUB = (double)sol_size;
        bestTerms = sol_size;
        bestWeight = w1;
        best_sol_size = sol_size;
        Memcpy(solution, sol_tmp, sol_size);
    }
    else {
        bestUB = (double)sol_size2;
        bestTerms = sol_size2;
        bestWeight = w2;
        best_sol_size = sol_size2;
        Memcpy(solution, sol_tmp2, sol_size2);
    }
    double bestLB = -DBL_MAX;
    for (int it = 0; it < max_iter; ++it) {
        maybe_check_user_interrupt(it);
        double ZLB = compute_reduced_and_lb(rows, cols, rowsCovered, rowsCoveredCount, t, ls);
        double LBint = ceil(ZLB - 1e-12);
        if (LBint > bestLB + EPS) {
            bestLB = LBint;
            stuck_iter = 0;
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
            if (local_passes > 0) {
                if (sol_size > 0) {
                    local_search_drop1_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp, &sol_size, local_passes
                    );
                    local_search_drop2_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp, &sol_size, local_passes
                    );
                }
                if (sol_size2 > 0) {
                    local_search_drop1_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp2, &sol_size2, local_passes
                    );
                    local_search_drop2_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp2, &sol_size2, local_passes
                    );
                }
            }
            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
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
        }
        if (bestUB <= LBint + EPS) break;
        int updated = subgradient_update(
            rows, cols, colsCovering, colsCoveringCount,
            ls, bestUB, ZLB, t, &step_coef, step_min, &stuck_iter, halve_period
        );
        if (!updated || step_coef <= step_min + EPS) {
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
            if (local_passes > 0) {
                if (sol_size > 0) {
                    local_search_drop1_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp, &sol_size, local_passes
                    );
                    local_search_drop2_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp, &sol_size, local_passes
                    );
                }
                if (sol_size2 > 0) {
                    local_search_drop1_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp2, &sol_size2, local_passes
                    );
                    local_search_drop2_repair(
                        rows, cols, rowsCovered, rowsCoveredCount,
                        colsCovering, colsCoveringCount,
                        ls, weights, sol_tmp2, &sol_size2, local_passes
                    );
                }
            }
            double candW1 = (sol_size != -1) ? solution_total_weight(sol_tmp, sol_size, weights) : -DBL_MAX;
            double candW2 = (sol_size2 != -1) ? solution_total_weight(sol_tmp2, sol_size2, weights) : -DBL_MAX;
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
            break;
        }
    }
    *solmin = best_sol_size;
    R_Free(t);
    R_Free(ls);
    R_Free(sol_tmp);
    R_Free(sol_tmp2);
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
    solvePIchart_lagrangian(pichart, cols, rows, NULL, indices, &solmin);
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
