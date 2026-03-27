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

#include "solvePIchart_gurobi.h"
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
    if (!pichart || !indices) {
        if (pichart) R_Free(pichart);
        if (indices) R_Free(indices);
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
        return R_NilValue;
    }
    int solmin = 0;
    SEXP out = R_NilValue;
    if (solvePIchart_gurobi(pichart, foundPI, on_minterms, indices, &solmin)) {
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
#endif
