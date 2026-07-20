#include "qca_rinternals.h"
#include <stdlib.h>
#include <string.h>

#include "findmin_lpsolve.h"
#include "pichart_presolve.h"
#include "lp_solve/lp_lib.h"

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
    double *solution = (double *) calloc((size_t) ncols, sizeof(double));

    if (solution == NULL) {
        return FALSE;
    }

    if (!solve_lpsolve_from_int_matrix(chart, nrows, ncols, solution)) {
        free(solution);
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
    return TRUE;
}

Rboolean solvePIchart_lpsolve_active(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    int *indices,
    int *solmin
) {
    int compact_ncols = 0;
    for (int c = 0; c < ncols; ++c) {
        compact_ncols += active == NULL || active[c] != 0;
    }
    if (compact_ncols <= 0) return FALSE;
    if (compact_ncols == ncols) {
        return solvePIchart_lpsolve(chart, nrows, ncols, indices, solmin);
    }

    int *compact_map = (int *)malloc((size_t)compact_ncols * sizeof(int));
    int *compact_chart = (int *)malloc(
        (size_t)nrows * (size_t)compact_ncols * sizeof(int)
    );
    int *compact_indices = (int *)malloc((size_t)compact_ncols * sizeof(int));
    if (!compact_map || !compact_chart || !compact_indices) {
        free(compact_map);
        free(compact_chart);
        free(compact_indices);
        return FALSE;
    }

    for (int c = 0, cc = 0; c < ncols; ++c) {
        if (active != NULL && !active[c]) continue;
        compact_map[cc] = c;
        memcpy(
            &compact_chart[(size_t)cc * nrows],
            &chart[(size_t)c * nrows],
            (size_t)nrows * sizeof(int)
        );
        ++cc;
    }

    int compact_solmin = 0;
    Rboolean ok = solvePIchart_lpsolve(
        compact_chart,
        nrows,
        compact_ncols,
        compact_indices,
        &compact_solmin
    );
    if (ok) {
        *solmin = compact_solmin;
        for (int i = 0; i < compact_solmin; ++i) {
            int cc = compact_indices[i];
            if (cc < 0 || cc >= compact_ncols) {
                ok = FALSE;
                *solmin = 0;
                break;
            }
            indices[i] = compact_map[cc];
        }
    }

    free(compact_map);
    free(compact_chart);
    free(compact_indices);
    return ok;
}

SEXP C_findminLpSolveInternal(SEXP chart) {
    if (!isMatrix(chart) || TYPEOF(chart) != LGLSXP) {
        error("C_findminLpSolveInternal expects a logical matrix.");
    }

    const int nr = nrows(chart);
    const int nc = ncols(chart);
    const int *p_chart = LOGICAL(chart);
    SEXP out = PROTECT(allocVector(REALSXP, nc));
    unsigned char *active = (unsigned char *)malloc((size_t)nc);
    int *indices = (int *)calloc((size_t)nc, sizeof(int));

    if (active == NULL || indices == NULL) {
        free(active);
        free(indices);
        UNPROTECT(1);
        error("Failed to allocate internal lp_solve workspace.");
    }
    memset(active, 1, (size_t)nc);
    qca_reduce_active_columns(p_chart, nr, nc, active);

    int solmin = 0;
    if (!solvePIchart_lpsolve_active(
        p_chart, nr, nc, active, indices, &solmin
    )) {
        free(active);
        free(indices);
        UNPROTECT(1);
        error("Internal lp_solve failed to find an optimal solution.");
    }

    for (int c = 0; c < nc; c++) {
        SET_REAL_ELT(out, c, 0.0);
    }
    for (int i = 0; i < solmin; ++i) {
        SET_REAL_ELT(out, indices[i], 1.0);
    }
    free(active);
    free(indices);

    UNPROTECT(1);
    return out;
}
