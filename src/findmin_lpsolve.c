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

#include <R.h>
#include <R_ext/Boolean.h>
#include <Rinternals.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "findmin_lpsolve.h"
#include "lp_solve/lp_lib.h"
#include "solvePIchart_lagrangian.h"
static int solve_lpsolve_from_int_matrix(
    const int *p_chart,
    const int nr,
    const int nc,
    double *p_out
) {
    lprec *lp = NULL;
    REAL *objective = NULL;
    REAL *row = NULL;
    lp = make_lp(0, nc);
    if (lp == NULL) {
        return 0;
    }
    objective = (REAL *) calloc((size_t) nc + 1, sizeof(REAL));
    row = (REAL *) calloc((size_t) nc + 1, sizeof(REAL));
    if (objective == NULL || row == NULL) {
        delete_lp(lp);
        free(objective);
        free(row);
        return 0;
    }
    set_verbose(lp, CRITICAL);
    set_minim(lp);
    set_add_rowmode(lp, TRUE);
    set_scaling(lp, 196);
    for (int c = 1; c <= nc; c++) {
        objective[c] = 1.0;
    }
    if (!set_obj_fn(lp, objective)) {
        delete_lp(lp);
        free(objective);
        free(row);
        return 0;
    }
    for (int r = 0; r < nr; r++) {
        memset(row, 0, ((size_t) nc + 1) * sizeof(REAL));
        for (int c = 0; c < nc; c++) {
            row[c + 1] = p_chart[c * nr + r] ? 1.0 : 0.0;
        }
        if (!add_constraint(lp, row, GE, 1.0)) {
            delete_lp(lp);
            free(objective);
            free(row);
            return 0;
        }
    }
    set_add_rowmode(lp, FALSE);
    for (int c = 1; c <= nc; c++) {
        if (!set_binary(lp, c, TRUE)) {
            delete_lp(lp);
            free(objective);
            free(row);
            return 0;
        }
    }
    if (solve(lp) != OPTIMAL) {
        delete_lp(lp);
        free(objective);
        free(row);
        return 0;
    }
    if (!get_variables(lp, p_out)) {
        delete_lp(lp);
        free(objective);
        free(row);
        return 0;
    }
    delete_lp(lp);
    free(objective);
    free(row);
    return 1;
}
Rboolean solvePIchart_lpsolve(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
) {
    int lagr_solmin = 0;
    double lagr_lb = -1e308;
    double *solution = (double *) calloc((size_t) ncols, sizeof(double));
    int *lagr_indices = (int *) calloc((size_t) ncols, sizeof(int));
    if (solution == NULL || lagr_indices == NULL) {
        free(solution);
        free(lagr_indices);
        return FALSE;
    }
    solvePIchart_lagrangian(
        (int *) chart,
        ncols,
        nrows,
        NULL,
        lagr_indices,
        &lagr_solmin,
        &lagr_lb,
        NULL
    );
    if (lagr_solmin > 0 && lagr_lb > -1e307) {
        double lagr_lb_int = ceil(lagr_lb - 1e-12);
        if ((double) lagr_solmin <= lagr_lb_int + 1e-12) {
            for (int i = 0; i < lagr_solmin; i++) {
                indices[i] = lagr_indices[i];
            }
            *solmin = lagr_solmin;
            free(solution);
            free(lagr_indices);
            return TRUE;
        }
    }
    if (!solve_lpsolve_from_int_matrix(chart, nrows, ncols, solution)) {
        free(solution);
        free(lagr_indices);
        return FALSE;
    }
    *solmin = 0;
    for (int c = 0; c < ncols; c++) {
        if (solution[c] > 0.5) {
            indices[*solmin] = c;
            (*solmin)++;
        }
    }
    free(solution);
    free(lagr_indices);
    return TRUE;
}
SEXP C_findminLpSolveInternal(SEXP chart) {
    if (!isMatrix(chart) || TYPEOF(chart) != LGLSXP) {
        error("C_findminLpSolveInternal expects a logical matrix.");
    }
    const int nr = nrows(chart);
    const int nc = ncols(chart);
    const int *p_chart = LOGICAL(chart);
    SEXP out = PROTECT(allocVector(REALSXP, nc));
    double *tmp_out = (double *) calloc((size_t) nc, sizeof(double));
    int lagr_solmin = 0;
    double lagr_lb = -1e308;
    int *lagr_indices = (int *) calloc((size_t) nc, sizeof(int));
    if (tmp_out == NULL || lagr_indices == NULL) {
        free(tmp_out);
        free(lagr_indices);
        UNPROTECT(1);
        error("Failed to allocate Lagrangian workspace for internal lp_solve.");
    }
    solvePIchart_lagrangian(
        (int *) p_chart,
        nc,
        nr,
        NULL,
        lagr_indices,
        &lagr_solmin,
        &lagr_lb,
        NULL
    );
    if (lagr_solmin > 0 && lagr_lb > -1e307) {
        double lagr_lb_int = ceil(lagr_lb - 1e-12);
        if ((double) lagr_solmin <= lagr_lb_int + 1e-12) {
            for (int c = 0; c < nc; c++) {
                tmp_out[c] = 0.0;
            }
            for (int i = 0; i < lagr_solmin; i++) {
                int col = lagr_indices[i];
                if (col >= 0 && col < nc) {
                    tmp_out[col] = 1.0;
                }
            }
            for (int c = 0; c < nc; c++) {
                SET_REAL_ELT(out, c, tmp_out[c]);
            }
            free(tmp_out);
            free(lagr_indices);
            UNPROTECT(1);
            return out;
        }
    }
    free(lagr_indices);
    if (!solve_lpsolve_from_int_matrix(p_chart, nr, nc, tmp_out)) {
        free(tmp_out);
        UNPROTECT(1);
        error("Internal lp_solve failed to find an optimal solution.");
    }
    for (int c = 0; c < nc; c++) {
        SET_REAL_ELT(out, c, tmp_out[c]);
    }
    free(tmp_out);
    UNPROTECT(1);
    return out;
}
