#include "findmin_gurobi.h"
#include "pichart_presolve.h"
#include <string.h>

#ifdef HAVE_GUROBI

#include <R_ext/RS.h>
#include <fcntl.h>
#include <stdio.h>
#include <unistd.h>
#include "gurobi_c.h"

static GRBenv *qca_gurobi_env = NULL;

static bool qca_gurobi_env_ready(void) {
    if (qca_gurobi_env != NULL) {
        return true;
    }

    int error = 0;
    FILE *nullout = fopen("/dev/null", "w");
    int saved_stdout = -1;
    int saved_stderr = -1;
    int null_fd = -1;

    if (nullout) {
        null_fd = fileno(nullout);
        saved_stdout = dup(STDOUT_FILENO);
        saved_stderr = dup(STDERR_FILENO);
        dup2(null_fd, STDOUT_FILENO);
        dup2(null_fd, STDERR_FILENO);
    }

    error = GRBloadenv(&qca_gurobi_env, "/dev/null");

    if (nullout) {
        fflush(NULL);

        if (saved_stdout != -1) {
            dup2(saved_stdout, STDOUT_FILENO);
            close(saved_stdout);
        }

        if (saved_stderr != -1) {
            dup2(saved_stderr, STDERR_FILENO);
            close(saved_stderr);
        }

        fclose(nullout);
    }

    if (error) {
        if (qca_gurobi_env) {
            GRBfreeenv(qca_gurobi_env);
            qca_gurobi_env = NULL;
        }
        return false;
    }

    error = GRBsetintparam(qca_gurobi_env, "OutputFlag", 0);
    if (error) {
        GRBfreeenv(qca_gurobi_env);
        qca_gurobi_env = NULL;
        return false;
    }

    return true;
}

bool gurobi_runtime_available(void) {
    return qca_gurobi_env_ready();
}

void gurobi_release_env(void) {
    if (qca_gurobi_env != NULL) {
        GRBfreeenv(qca_gurobi_env);
        qca_gurobi_env = NULL;
    }
}

SEXP C_gurobiRuntimeAvailable(void) {
    SEXP out = PROTECT(allocVector(LGLSXP, 1));
    LOGICAL(out)[0] = gurobi_runtime_available();
    UNPROTECT(1);
    return out;
}

SEXP C_findminExact(SEXP chart) {
    if (!isMatrix(chart)) {
        return R_NilValue;
    }

    SEXP dims = getAttrib(chart, R_DimSymbol);
    int on_minterms = INTEGER(dims)[0];
    int foundPI = INTEGER(dims)[1];

    if (on_minterms < 1 || foundPI < 1) {
        return R_NilValue;
    }

    int *pichart = (int *) R_Calloc((size_t) on_minterms * (size_t) foundPI, int);
    int *indices = (int *) R_Calloc((size_t) foundPI, int);
    unsigned char *active = (unsigned char *)R_Calloc((size_t)foundPI, unsigned char);
    if (!pichart || !indices || !active) {
        if (pichart) R_Free(pichart);
        if (indices) R_Free(indices);
        if (active) R_Free(active);
        return R_NilValue;
    }

    if (TYPEOF(chart) == LGLSXP) {
        int *src = LOGICAL(chart);
        for (int i = 0; i < on_minterms * foundPI; i++) {
            pichart[i] = src[i] != 0;
        }
    }
    else if (TYPEOF(chart) == INTSXP) {
        int *src = INTEGER(chart);
        for (int i = 0; i < on_minterms * foundPI; i++) {
            pichart[i] = src[i] != 0;
        }
    }
    else {
        R_Free(pichart);
        R_Free(indices);
        R_Free(active);
        return R_NilValue;
    }

    memset(active, 1, (size_t)foundPI);
    qca_reduce_active_columns(pichart, on_minterms, foundPI, active);

    int solmin = 0;
    SEXP out = R_NilValue;

    if (solvePIchart_gurobi_active(
        pichart, foundPI, on_minterms, active, indices, &solmin
    )) {
        out = PROTECT(allocVector(REALSXP, foundPI));
        for (int j = 0; j < foundPI; j++) {
            REAL(out)[j] = 0.0;
        }
        for (int j = 0; j < solmin; j++) {
            REAL(out)[indices[j]] = 1.0;
        }
        UNPROTECT(1);
    }

    R_Free(pichart);
    R_Free(indices);
    R_Free(active);
    return out;
}

bool solvePIchart_gurobi(
    const int pichart[],
    int foundPI,
    int on_minterms,
    int indices[],
    int *solmin
) {
    int error = 0;
    GRBmodel *model = NULL;
    int *ind = NULL;
    double *coeffs = NULL;
    double *solution = NULL;

    if (solmin) {
        *solmin = 0;
    }

    ind = (int *) R_Calloc((size_t) foundPI, int);
    coeffs = (double *) R_Calloc((size_t) foundPI, double);
    solution = (double *) R_Calloc((size_t) foundPI, double);

    if (!ind || !coeffs || !solution) {
        error = 1;
        goto QUIT;
    }

    if (!qca_gurobi_env_ready()) goto QUIT;

    error = GRBnewmodel(
        qca_gurobi_env,
        &model,
        "QCASetCover",
        foundPI,
        NULL, NULL, NULL, NULL, NULL
    );
    if (error) goto QUIT;

    for (int j = 0; j < foundPI; j++) {
        error = GRBsetcharattrelement(model, GRB_CHAR_ATTR_VTYPE, j, GRB_BINARY);
        if (error) goto QUIT;
    }

    for (int i = 0; i < on_minterms; i++) {
        int nz = 0;
        for (int j = 0; j < foundPI; j++) {
            if (pichart[i + on_minterms * j] == 1) {
                ind[nz] = j;
                coeffs[nz] = 1.0;
                nz++;
            }
        }

        error = GRBaddconstr(model, nz, ind, coeffs, GRB_GREATER_EQUAL, 1.0, NULL);
        if (error) goto QUIT;
    }

    for (int j = 0; j < foundPI; j++) {
        ind[j] = j;
        coeffs[j] = 1.0;
    }

    error = GRBsetobjectiven(
        model,
        0,
        1,
        1.0,
        0.0,
        0.0,
        "mincols",
        0.0,
        foundPI,
        ind,
        coeffs
    );
    if (error) goto QUIT;

    error = GRBoptimize(model);
    if (error) goto QUIT;

    double objval = 0.0;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;

    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, foundPI, solution);
    if (error) goto QUIT;

    if (solmin) {
        *solmin = (int) objval;
    }

    int pos = 0;
    for (int j = 0; j < foundPI; j++) {
        if (solution[j] > 0.9) {
            indices[pos] = j;
            pos++;
        }
    }

QUIT:
    if (ind) R_Free(ind);
    if (coeffs) R_Free(coeffs);
    if (solution) R_Free(solution);
    if (model) GRBfreemodel(model);

    return error == 0 && solmin != NULL && *solmin > 0;
}

static bool solvePIchart_gurobi_with_incumbent(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const int initial_indices[],
    int initial_solmin,
    int indices[],
    int *solmin
) {
    int error = 0;
    GRBmodel *model = NULL;
    int *ind = NULL;
    double *coeffs = NULL;
    double *solution = NULL;
    double *start = NULL;
    if (solmin) *solmin = 0;

    ind = (int *)R_Calloc((size_t)foundPI, int);
    coeffs = (double *)R_Calloc((size_t)foundPI, double);
    solution = (double *)R_Calloc((size_t)foundPI, double);
    start = (double *)R_Calloc((size_t)foundPI, double);
    if (!ind || !coeffs || !solution || !start || !qca_gurobi_env_ready()) {
        error = 1;
        goto QUIT;
    }

    error = GRBnewmodel(qca_gurobi_env, &model, "QCASetCover",
                        foundPI, NULL, NULL, NULL, NULL, NULL);
    if (error) goto QUIT;
    for (int j = 0; j < foundPI; ++j) {
        error = GRBsetcharattrelement(model, GRB_CHAR_ATTR_VTYPE, j, GRB_BINARY);
        if (error) goto QUIT;
    }
    for (int i = 0; i < on_minterms; ++i) {
        int nz = 0;
        for (int j = 0; j < foundPI; ++j) {
            if (pichart[i + on_minterms * j]) {
                ind[nz] = j;
                coeffs[nz++] = 1.0;
            }
        }
        error = GRBaddconstr(model, nz, ind, coeffs, GRB_GREATER_EQUAL, 1.0, NULL);
        if (error) goto QUIT;
    }
    for (int j = 0; j < foundPI; ++j) {
        ind[j] = j;
        coeffs[j] = 1.0;
        start[j] = 0.0;
    }
    error = GRBsetobjectiven(model, 0, 1, 1.0, 0.0, 0.0,
                             "mincols", 0.0, foundPI, ind, coeffs);
    if (error) goto QUIT;
    for (int i = 0; i < initial_solmin; ++i) start[initial_indices[i]] = 1.0;
    error = GRBsetdblattrarray(model, GRB_DBL_ATTR_START, 0, foundPI, start);
    if (error) goto QUIT;
    error = GRBoptimize(model);
    if (error) goto QUIT;
    double objval = 0.0;
    error = GRBgetdblattr(model, GRB_DBL_ATTR_OBJVAL, &objval);
    if (error) goto QUIT;
    error = GRBgetdblattrarray(model, GRB_DBL_ATTR_X, 0, foundPI, solution);
    if (error) goto QUIT;
    *solmin = (int)objval;
    for (int j = 0, pos = 0; j < foundPI; ++j) {
        if (solution[j] > 0.9) indices[pos++] = j;
    }

QUIT:
    if (ind) R_Free(ind);
    if (coeffs) R_Free(coeffs);
    if (solution) R_Free(solution);
    if (start) R_Free(start);
    if (model) GRBfreemodel(model);
    return error == 0 && solmin && *solmin > 0;
}

bool solvePIchart_gurobi_active(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const unsigned char active[],
    int indices[],
    int *solmin
) {
    int compact_foundPI = 0;
    for (int c = 0; c < foundPI; ++c) {
        compact_foundPI += active == NULL || active[c] != 0;
    }
    if (compact_foundPI <= 0) return false;
    if (compact_foundPI == foundPI) {
        return solvePIchart_gurobi(
            pichart, foundPI, on_minterms, indices, solmin
        );
    }

    int *compact_map = (int *)R_Calloc((size_t)compact_foundPI, int);
    int *compact_chart = (int *)R_Calloc(
        (size_t)on_minterms * (size_t)compact_foundPI, int
    );
    int *compact_indices = (int *)R_Calloc((size_t)compact_foundPI, int);
    if (!compact_map || !compact_chart || !compact_indices) {
        if (compact_map) R_Free(compact_map);
        if (compact_chart) R_Free(compact_chart);
        if (compact_indices) R_Free(compact_indices);
        return false;
    }

    for (int c = 0, cc = 0; c < foundPI; ++c) {
        if (active != NULL && !active[c]) continue;
        compact_map[cc] = c;
        Memcpy(
            &compact_chart[(size_t)cc * on_minterms],
            &pichart[(size_t)c * on_minterms],
            on_minterms
        );
        ++cc;
    }

    int compact_solmin = 0;
    bool ok = solvePIchart_gurobi(
        compact_chart,
        compact_foundPI,
        on_minterms,
        compact_indices,
        &compact_solmin
    );
    if (ok) {
        *solmin = compact_solmin;
        for (int i = 0; i < compact_solmin; ++i) {
            int cc = compact_indices[i];
            if (cc < 0 || cc >= compact_foundPI) {
                ok = false;
                *solmin = 0;
                break;
            }
            indices[i] = compact_map[cc];
        }
    }

    R_Free(compact_map);
    R_Free(compact_chart);
    R_Free(compact_indices);
    return ok;
}

bool solvePIchart_gurobi_active_with_incumbent(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const unsigned char active[],
    const int initial_indices[],
    int initial_solmin,
    int indices[],
    int *solmin
) {
    if (!initial_indices || initial_solmin <= 0) {
        return solvePIchart_gurobi_active(
            pichart, foundPI, on_minterms, active, indices, solmin
        );
    }
    int compact_foundPI = 0;
    for (int c = 0; c < foundPI; ++c) {
        compact_foundPI += active == NULL || active[c] != 0;
    }
    int *map = (int *)R_Calloc((size_t)compact_foundPI, int);
    int *reverse = (int *)R_Calloc((size_t)foundPI, int);
    int *compact_chart = (int *)R_Calloc(
        (size_t)on_minterms * (size_t)compact_foundPI, int
    );
    int *compact_initial = (int *)R_Calloc((size_t)initial_solmin, int);
    int *compact_indices = (int *)R_Calloc((size_t)compact_foundPI, int);
    bool ok = false;
    if (!map || !reverse || !compact_chart || !compact_initial || !compact_indices) goto DONE;
    for (int c = 0; c < foundPI; ++c) reverse[c] = -1;
    for (int c = 0, cc = 0; c < foundPI; ++c) {
        if (active && !active[c]) continue;
        map[cc] = c;
        reverse[c] = cc;
        Memcpy(&compact_chart[(size_t)cc * on_minterms],
               &pichart[(size_t)c * on_minterms], on_minterms);
        ++cc;
    }
    for (int i = 0; i < initial_solmin; ++i) {
        int col = initial_indices[i];
        if (col < 0 || col >= foundPI || reverse[col] < 0) goto DONE;
        compact_initial[i] = reverse[col];
    }
    int compact_solmin = 0;
    ok = solvePIchart_gurobi_with_incumbent(
        compact_chart, compact_foundPI, on_minterms,
        compact_initial, initial_solmin, compact_indices, &compact_solmin
    );
    if (ok) {
        *solmin = compact_solmin;
        for (int i = 0; i < compact_solmin; ++i) indices[i] = map[compact_indices[i]];
    }
DONE:
    if (map) R_Free(map);
    if (reverse) R_Free(reverse);
    if (compact_chart) R_Free(compact_chart);
    if (compact_initial) R_Free(compact_initial);
    if (compact_indices) R_Free(compact_indices);
    return ok;
}

#else

bool gurobi_runtime_available(void) {
    return false;
}

void gurobi_release_env(void) {
}

SEXP C_gurobiRuntimeAvailable(void) {
    SEXP out = PROTECT(allocVector(LGLSXP, 1));
    LOGICAL(out)[0] = 0;
    UNPROTECT(1);
    return out;
}

SEXP C_findminExact(SEXP chart) {
    (void) chart;
    return R_NilValue;
}

bool solvePIchart_gurobi(
    const int pichart[],
    int foundPI,
    int on_minterms,
    int indices[],
    int *solmin
) {
    (void) pichart;
    (void) foundPI;
    (void) on_minterms;
    (void) indices;
    if (solmin) {
        *solmin = 0;
    }
    return false;
}

bool solvePIchart_gurobi_active(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const unsigned char active[],
    int indices[],
    int *solmin
) {
    (void) pichart;
    (void) foundPI;
    (void) on_minterms;
    (void) active;
    (void) indices;
    if (solmin) *solmin = 0;
    return false;
}

bool solvePIchart_gurobi_active_with_incumbent(
    const int pichart[], int foundPI, int on_minterms,
    const unsigned char active[], const int initial_indices[],
    int initial_solmin, int indices[], int *solmin
) {
    (void)initial_indices;
    (void)initial_solmin;
    return solvePIchart_gurobi_active(
        pichart, foundPI, on_minterms, active, indices, solmin
    );
}

#endif
