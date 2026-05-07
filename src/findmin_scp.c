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
#include <R_ext/RS.h>
#include <Rinternals.h>
#include <string.h>
#include "findmin_scp.h"
#include "scp_solver/scp_solver.h"
#include "solvePIchart_lagrangian.h"
static int solve_scp_from_int_matrix(
    const int *p_chart,
    const int nr,
    const int nc,
    int *solution
) {
    const int nwords_rows = (nr + 63) / 64;
    int nnz = 0;
    for (int c = 0; c < nc; c++) {
        for (int r = 0; r < nr; r++) {
            if (p_chart[c * nr + r]) {
                nnz++;
            }
        }
    }
    int *row_counts = (int *) R_Calloc((size_t) nr, int);
    int *row_starts = (int *) R_Calloc((size_t) nr + 1, int);
    int *row_cols = (int *) R_Calloc((size_t) nnz, int);
    unsigned long long *col_masks = (unsigned long long *) R_Calloc((size_t) nc * nwords_rows, unsigned long long);
    if (row_counts == NULL || row_starts == NULL || row_cols == NULL || col_masks == NULL || solution == NULL) {
        if (row_counts) R_Free(row_counts);
        if (row_starts) R_Free(row_starts);
        if (row_cols) R_Free(row_cols);
        if (col_masks) R_Free(col_masks);
        if (solution) R_Free(solution);
        error("Failed to allocate SCP solver workspace.");
    }
    for (int r = 0; r < nr; r++) {
        for (int c = 0; c < nc; c++) {
            if (p_chart[c * nr + r]) {
                row_counts[r]++;
                col_masks[c * nwords_rows + (r >> 6)] |= (1ULL << (r & 63));
            }
        }
    }
    row_starts[0] = 0;
    for (int r = 0; r < nr; r++) {
        row_starts[r + 1] = row_starts[r] + row_counts[r];
    }
    memset(row_counts, 0, (size_t) nr * sizeof(int));
    for (int r = 0; r < nr; r++) {
        for (int c = 0; c < nc; c++) {
            if (p_chart[c * nr + r]) {
                int pos = row_starts[r] + row_counts[r];
                row_cols[pos] = c;
                row_counts[r]++;
            }
        }
    }
    qca_scp_problem problem;
    problem.nr = nr;
    problem.nc = nc;
    problem.nwords_rows = nwords_rows;
    problem.row_starts = row_starts;
    problem.row_cols = row_cols;
    problem.col_masks = col_masks;
    problem.branch_priority = NULL;
    int solution_size = 0;
    int incumbent_size = -1;
    int *incumbent_indices = (int *) R_Calloc((size_t) nc, int);
    int *incumbent_solution = (int *) R_Calloc((size_t) nc, int);
    double *lagr_scores = (double *) R_Calloc((size_t) nc, double);
    double lagr_lb = -1e308;
    int proof_target_size = -1;
    int ok = 0;
    if (incumbent_indices == NULL || incumbent_solution == NULL || lagr_scores == NULL) {
        R_Free(row_counts);
        R_Free(row_starts);
        R_Free(row_cols);
        R_Free(col_masks);
        if (incumbent_indices) R_Free(incumbent_indices);
        if (incumbent_solution) R_Free(incumbent_solution);
        if (lagr_scores) R_Free(lagr_scores);
        if (solution) R_Free(solution);
        error("Failed to allocate SCP incumbent workspace.");
    }
    solvePIchart_lagrangian((int *) p_chart, nc, nr, NULL, incumbent_indices, &incumbent_size, &lagr_lb, lagr_scores);
    if (incumbent_size > 0 && incumbent_size <= nc) {
        for (int i = 0; i < incumbent_size; i++) {
            int col = incumbent_indices[i];
            if (col >= 0 && col < nc) {
                incumbent_solution[col] = 1;
            }
        }
    }
    if (incumbent_size > 0 && lagr_lb > -1e307) {
        double lagr_lb_int = ceil(lagr_lb - 1e-12);
        if ((double) incumbent_size <= lagr_lb_int + 1e-12) {
            memcpy(solution, incumbent_solution, (size_t) nc * sizeof(int));
            R_Free(row_counts);
            R_Free(row_starts);
            R_Free(row_cols);
            R_Free(col_masks);
            R_Free(incumbent_indices);
            R_Free(incumbent_solution);
            R_Free(lagr_scores);
            return 1;
        }
        if ((double) incumbent_size <= lagr_lb_int + 1.0 + 1e-12) {
            proof_target_size = incumbent_size - 1;
        }
    }
    problem.branch_priority = lagr_scores;
    ok = qca_scp_solve_exact_with_incumbent(
        &problem,
        solution,
        &solution_size,
        (incumbent_size > 0 && incumbent_size <= nc) ? incumbent_solution : NULL,
        (incumbent_size > 0 && incumbent_size <= nc) ? incumbent_size : 0,
        proof_target_size
    );
    R_Free(row_counts);
    R_Free(row_starts);
    R_Free(row_cols);
    R_Free(col_masks);
    R_Free(incumbent_indices);
    R_Free(incumbent_solution);
    R_Free(lagr_scores);
    return ok;
}
Rboolean solvePIchart_scp(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
) {
    int *solution = (int *) R_Calloc((size_t) ncols, int);
    if (solution == NULL) {
        return FALSE;
    }
    if (!solve_scp_from_int_matrix(chart, nrows, ncols, solution)) {
        R_Free(solution);
        return FALSE;
    }
    *solmin = 0;
    for (int c = 0; c < ncols; c++) {
        if (solution[c] > 0) {
            indices[*solmin] = c;
            (*solmin)++;
        }
    }
    R_Free(solution);
    return TRUE;
}
SEXP C_getScpProfile(void) {
    qca_scp_profile profile = qca_scp_profile_get();
    SEXP out = PROTECT(allocVector(VECSXP, 10));
    SEXP names = PROTECT(allocVector(STRSXP, 10));
    SET_STRING_ELT(names, 0, mkChar("total_seconds"));
    SET_STRING_ELT(names, 1, mkChar("reductions_seconds"));
    SET_STRING_ELT(names, 2, mkChar("lower_bound_seconds"));
    SET_STRING_ELT(names, 3, mkChar("greedy_seconds"));
    SET_STRING_ELT(names, 4, mkChar("branching_seconds"));
    SET_STRING_ELT(names, 5, mkChar("nodes"));
    SET_STRING_ELT(names, 6, mkChar("leaves"));
    SET_STRING_ELT(names, 7, mkChar("reduction_calls"));
    SET_STRING_ELT(names, 8, mkChar("lower_bound_calls"));
    SET_STRING_ELT(names, 9, mkChar("greedy_calls"));
    SET_VECTOR_ELT(out, 0, ScalarReal(profile.total_seconds));
    SET_VECTOR_ELT(out, 1, ScalarReal(profile.reductions_seconds));
    SET_VECTOR_ELT(out, 2, ScalarReal(profile.lower_bound_seconds));
    SET_VECTOR_ELT(out, 3, ScalarReal(profile.greedy_seconds));
    SET_VECTOR_ELT(out, 4, ScalarReal(profile.branching_seconds));
    SET_VECTOR_ELT(out, 5, ScalarReal((double) profile.nodes));
    SET_VECTOR_ELT(out, 6, ScalarReal((double) profile.leaves));
    SET_VECTOR_ELT(out, 7, ScalarReal((double) profile.reduction_calls));
    SET_VECTOR_ELT(out, 8, ScalarReal((double) profile.lower_bound_calls));
    SET_VECTOR_ELT(out, 9, ScalarReal((double) profile.greedy_calls));
    setAttrib(out, R_NamesSymbol, names);
    UNPROTECT(2);
    return out;
}
SEXP C_resetScpProfile(void) {
    qca_scp_profile_reset();
    return R_NilValue;
}
SEXP C_findminScpInternal(SEXP chart) {
    if (!isMatrix(chart) || TYPEOF(chart) != LGLSXP) {
        error("C_findminScpInternal expects a logical matrix.");
    }
    const int nr = nrows(chart);
    const int nc = ncols(chart);
    const int *p_chart = LOGICAL(chart);
    int *solution = (int *) R_Calloc((size_t) nc, int);
    if (solution == NULL) {
        error("Failed to allocate SCP solver solution vector.");
    }
    if (!solve_scp_from_int_matrix(p_chart, nr, nc, solution)) {
        R_Free(solution);
        error("Internal SCP solver failed to find a feasible exact cover.");
    }
    SEXP out = PROTECT(allocVector(REALSXP, nc));
    double *p_out = REAL(out);
    for (int c = 0; c < nc; c++) {
        p_out[c] = (double) solution[c];
    }
    R_Free(solution);
    UNPROTECT(1);
    return out;
}
