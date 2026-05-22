#include <R_ext/RS.h>
#include "qca_rinternals.h"
#include <float.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "qca_rinternals.h"
#include "qca_rinternals.h"
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include "utils.h"
#include "find_models.h"
#include "generate_matrix.h"
#include "qca_threads.h"
#include "sort_matrix.h"
#include "CCubes.h"

// IMPORTANT!
// memset() works <ONLY> with FALSE logical vectors...!!


static R_INLINE SEXP Rtranspose(SEXP matrix) {

    SEXPTYPE type = TYPEOF(matrix);
    int nr = nrows(matrix);
    int nc = ncols(matrix);

    // some matrices are preallocated and only the first part is filled
    // only those columns which are filled should be transposed
    SEXP lastcol = PROTECT(mkString("last_column"));
    SEXP ncm = PROTECT(getAttrib(matrix, lastcol));
    if (!Rf_isNull(ncm)) {
        nc = INTEGER(ncm)[0];
    }

    SEXP out = PROTECT(allocMatrix(type, nc, nr));
    R_xlen_t len = nr * nc;
    R_xlen_t i, j, l_1 = len - 1;

    if (type == INTSXP) {

        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            INTEGER(out)[i] = INTEGER(matrix)[j];
        }

        // int *p_out = INTEGER(out);
        // for (int r = 0; r < nr; r++) {
        //     for (int c = 0; c < nc; c++) {
        //         p_out[r * nc + c] = INTEGER(matrix)[c * nr + r];
        //     }
        // }
    }
    else if (type == LGLSXP) {

        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            LOGICAL(out)[i] = LOGICAL(matrix)[j];
        }

        // int *p_out = LOGICAL(out);
        // for (int r = 0; r < nr; r++) {
        //     for (int c = 0; c < nc; c++) {
        //         p_out[r * nc + c] = INTEGER(matrix)[c * nr + r];
        //     }
        // }
    }
    else if (type == REALSXP) {

        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            REAL(out)[i] = REAL(matrix)[j];
        }

        // double *p_out = REAL(out);
        // for (int r = 0; r < nr; r++) {
        //     for (int c = 0; c < nc; c++) {
        //         p_out[r * nc + c] = REAL(matrix)[c * nr + r];
        //     }
        // }
    }

    UNPROTECT(3);
    return(out);
}


static R_INLINE Rboolean hasColnames(SEXP matrix) {

    if (Rf_isNull(getAttrib(matrix, R_DimNamesSymbol))) {
        // does not have dimnames
        return(FALSE);
    }

    return !Rf_isNull(VECTOR_ELT(getAttrib(matrix, R_DimNamesSymbol), 1));

}

static R_INLINE Rboolean getpos(SEXP list, const char *str) {
    SEXP names = getAttrib(list, R_NamesSymbol);
    int pos = -1;

    if (!Rf_isNull(names)) {
        for (int i = 0; i < length(list); i++) {
            if (strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
                pos = i;
                break;
            }
        }
    }

    return(pos);
}

static R_INLINE SEXP Rresize(SEXP obj, int len) {

    SEXP usage = PROTECT(allocVector(VECSXP, 2));
    SEXP copy;

    int oldlen = length(obj);
    int copylen = (oldlen < len) ? oldlen : len;

    Rboolean objlogical = isLogical(obj); // TYPEOF(obj) == LGLSXP

    SET_VECTOR_ELT(usage, 0, copy = duplicate(obj));
    int *p_copy = INTEGER(copy);

    if (isMatrix(obj)) {
        int rows = nrows(obj);
        int cols = len / rows;
        SET_VECTOR_ELT(usage, 1, obj = allocMatrix(objlogical ? LGLSXP : INTSXP, rows, cols));
    }
    else {
        SET_VECTOR_ELT(usage, 1, obj = allocVector(objlogical ? LGLSXP : INTSXP, len));
    }

    int *p_obj = objlogical ? LOGICAL(obj) : INTEGER(obj);

    if (len > oldlen) {
        memset(p_obj, objlogical ? FALSE : 0, len * sizeof(int));
    }

    Memcpy(p_obj, p_copy, copylen);

    UNPROTECT(1);
    return(obj);
}

static int qca_compare_int(const void *a, const void *b) {
    int ia = *(const int *) a;
    int ib = *(const int *) b;

    return (ia > ib) - (ia < ib);
}

static Rboolean qca_int_in_sorted(const int value, const int values[], const int n) {
    int lo = 0;
    int hi = n - 1;

    while (lo <= hi) {
        int mid = lo + (hi - lo) / 2;

        if (values[mid] == value) {
            return TRUE;
        }
        if (values[mid] < value) {
            lo = mid + 1;
        }
        else {
            hi = mid - 1;
        }
    }

    return FALSE;
}

static R_INLINE SEXP Runique(SEXP mat) {
    // This function implicitly assumes that mat has at least 2 columns
    // and it is either integer or logical
    int nc = ncols(mat);
    int nr = nrows(mat);
    Rboolean logmat = isLogical(mat);
    int *p_mat = logmat ? LOGICAL(mat) : INTEGER(mat);
    SEXP umat;
    SEXP usage = PROTECT(allocVector(VECSXP, 1));

    Rboolean survived[nc];
    // memset() does NOT work with TRUE at Rboolean;
    // memset(survived, TRUE, nc * sizeof(int));

    int all_cols = nc;

    for (int i = 0; i < all_cols; i++) {
        survived[i] = TRUE;
    }

    for (int i = 0; i < all_cols; i++) {
        if (survived[i]) {
            for (int j = i + 1; j < all_cols; j++) {
                if (survived[j]) {
                    int r = 0;
                    Rboolean same = TRUE;
                    while (r < nr && same) {
                        same = p_mat[i * nr + r] == p_mat[j * nr + r];
                        r++;
                    }

                    if (same) {
                        survived[j] = FALSE;
                        --(nc);
                    }
                }
            }
        }
    }

    SET_VECTOR_ELT(usage, 0, umat = allocMatrix(logmat ? LGLSXP : INTSXP, nr, nc));
    int *p_umat = logmat ? LOGICAL(umat) : INTEGER(umat);

    int ci = 0;
    for (int c = 0; c < all_cols; c++) {
        if (survived[c]) {
            for (int r = 0; r < nr; r++) {
                p_umat[ci * nr + r] = p_mat[c * nr + r];
            }
            ci++;
        }
    }

    UNPROTECT(1);
    return(umat);
}

SEXP C_solveChart(SEXP pichart, SEXP allsol, SEXP vdepth, SEXP k, SEXP maxcomb, SEXP firstmin) {

    SEXP models = R_NilValue;
    SEXP usage = PROTECT(allocVector(VECSXP, 1));
    SEXP out = PROTECT(allocVector(VECSXP, 2));

    SET_VECTOR_ELT(usage, 0, pichart = Rtranspose(pichart));
    int *p_pichart = LOGICAL(pichart);

    int posrows = nrows(pichart);
    int foundPI = ncols(pichart);

    int *p_solutions = R_Calloc(1, int);
    int nr = 0;
    int nc = 0;

    find_models(
        p_pichart,
        posrows,
        foundPI,
        LOGICAL(allsol)[0],
        INTEGER(k)[0],
        REAL(maxcomb)[0],
        LOGICAL(firstmin)[0],
        &p_solutions,
        &nr,
        &nc
    );

    if (nr > 0 && nc > 0) {
        SET_VECTOR_ELT(out, 0, models = allocMatrix(INTSXP, nr, nc));
        // int *p_models = INTEGER(models);
        // for (int i = 0; i < nr * nc; i++) {
        //     p_models[i] = p_solutions[i];
        // }
        Memcpy(INTEGER(models), p_solutions, nr * nc);

        SEXP toocomplex;
        SET_VECTOR_ELT(out, 1, toocomplex = allocVector(LGLSXP, 1));
        LOGICAL(toocomplex)[0] = too_complex(foundPI, INTEGER(k)[0], REAL(maxcomb)[0]);
    }

    R_Free(p_solutions);
    UNPROTECT(2);
    return(out);

}

// void printfarray(int* arr, int size)
// {
//     for (int i = 0; i < size; i++)
//     {
//         Rprintf("%d ", arr[i]);
//     }
//     Rprintf("\n");
// }

//#define SHOW_DEBUG_OUTPUT        // This enable some textual output
//#define SHOW_DEBUG_PROFILE         // This enable time debugging

SEXP C_getRow(SEXP input) {
    PROTECT(input);
    SEXP rowno, noflevels, mbase, matrix;

    SEXP usage = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(usage, 0, rowno = coerceVector(VECTOR_ELT(input, 0), INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(VECTOR_ELT(input, 1), INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase   = coerceVector(VECTOR_ELT(input, 2), INTSXP));
    int *p_rowno = INTEGER(rowno);
    int *p_noflevels = INTEGER(noflevels);
    int *p_mbase = INTEGER(mbase);

    int nrows = length(rowno);
    int ncols = length(noflevels);

    SET_VECTOR_ELT(usage, 3, matrix = allocMatrix(INTSXP, nrows, ncols));
    int *p_matrix = INTEGER(matrix);

    // Rf_PrintValue(usage);

    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            p_matrix[c * nrows + r] = (p_rowno[r] / p_mbase[c]) % p_noflevels[c];
        }
    }

    UNPROTECT(2);
    return(matrix);
}

SEXP C_getSA(
    SEXP solution_list,
    SEXP expressions,
    SEXP noflevels,
    SEXP mbaseexpr,
    SEXP inputt,
    SEXP mbaseplus,
    SEXP mbase
) {

    SEXP usage = PROTECT(allocVector(VECSXP, 7));

    SET_VECTOR_ELT(usage, 0, expressions = coerceVector(expressions, INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbaseexpr = coerceVector(mbaseexpr, INTSXP));
    SET_VECTOR_ELT(usage, 3, inputt = coerceVector(inputt, INTSXP));
    SET_VECTOR_ELT(usage, 4, mbaseplus = coerceVector(mbaseplus, INTSXP));
    SET_VECTOR_ELT(usage, 5, mbase = coerceVector(mbase, INTSXP));

    if (!isNewList(solution_list)) {
        error("C_getSA expects a list of solutions.");
    }

    int *p_expr = INTEGER(expressions);
    int *p_noflevels = INTEGER(noflevels);
    int *p_mbaseexpr = INTEGER(mbaseexpr);
    int *p_inputt = INTEGER(inputt);
    int *p_mbaseplus = INTEGER(mbaseplus);
    int *p_mbase = INTEGER(mbase);

    int nsol = length(solution_list);
    int expr_rows = nrows(expressions);
    int expr_cols = ncols(expressions);
    int input_rows = nrows(inputt);
    int input_cols = ncols(inputt);

    SEXP out = PROTECT(allocVector(VECSXP, nsol));

    int *observed = (int *) R_Calloc((size_t) input_rows, int);
    if (observed == NULL) {
        error("Memory allocation failed in C_getSA().");
    }

    for (int r = 0; r < input_rows; r++) {
        int key = 0;
        for (int c = 0; c < input_cols; c++) {
            key += p_inputt[c * input_rows + r] * p_mbaseplus[c];
        }
        observed[r] = key;
    }

    qsort(observed, (size_t) input_rows, sizeof(int), qca_compare_int);

    SEXP input_dimnames = getAttrib(inputt, R_DimNamesSymbol);
    SEXP input_colnames = R_NilValue;
    if (!Rf_isNull(input_dimnames)) {
        input_colnames = VECTOR_ELT(input_dimnames, 1);
    }

    SEXP expr_dimnames = getAttrib(expressions, R_DimNamesSymbol);
    SEXP expr_rownames = R_NilValue;
    if (!Rf_isNull(expr_dimnames)) {
        expr_rownames = VECTOR_ELT(expr_dimnames, 0);
    }

    for (int s = 0; s < nsol; s++) {
        SEXP sol = VECTOR_ELT(solution_list, s);
        int sol_len = length(sol);

        size_t keys_capacity = 1024;
        size_t keys_count = 0;
        int *keys = (int *) R_Calloc(keys_capacity, int);
        int *zero_pos = (int *) R_Calloc((size_t) expr_cols, int);

        if (keys == NULL || zero_pos == NULL) {
            if (keys != NULL) R_Free(keys);
            if (zero_pos != NULL) R_Free(zero_pos);
            R_Free(observed);
            error("Memory allocation failed in C_getSA().");
        }

        for (int i = 0; i < sol_len; i++) {
            int expr_row = -1;

            if (TYPEOF(sol) == INTSXP) {
                expr_row = INTEGER(sol)[i] - 1;
            }
            else if (TYPEOF(sol) == REALSXP) {
                expr_row = (int) REAL(sol)[i] - 1;
            }
            else if (TYPEOF(sol) == STRSXP && !Rf_isNull(expr_rownames)) {
                const char *target = CHAR(STRING_ELT(sol, i));

                for (int r = 0; r < expr_rows && expr_row < 0; r++) {
                    if (strcmp(target, CHAR(STRING_ELT(expr_rownames, r))) == 0) {
                        expr_row = r;
                    }
                }
            }
            else {
                R_Free(keys);
                R_Free(zero_pos);
                R_Free(observed);
                error("Unsupported solution index type in C_getSA().");
            }

            if (expr_row < 0 || expr_row >= expr_rows) {
                R_Free(keys);
                R_Free(zero_pos);
                R_Free(observed);
                error("Solution index outside the expression matrix.");
            }

            int base = 0;
            int zero_count = 0;

            for (int c = 0; c < expr_cols; c++) {
                int value = p_expr[c * expr_rows + expr_row];
                base += p_mbaseexpr[c] * value;

                if (value == 0) {
                    base += p_mbaseexpr[c];
                    zero_pos[zero_count] = c;
                    zero_count++;
                }
            }

            if (zero_count == 0) {
                continue;
            }

            size_t combinations = 1;
            for (int z = 0; z < zero_count; z++) {
                int c = zero_pos[z];
                if (p_noflevels[c] <= 0 || combinations > SIZE_MAX / (size_t) p_noflevels[c]) {
                    R_Free(keys);
                    R_Free(zero_pos);
                    R_Free(observed);
                    error("Too many simplifying assumptions to enumerate.");
                }
                combinations *= (size_t) p_noflevels[c];
            }

            if (keys_count > SIZE_MAX - combinations) {
                R_Free(keys);
                R_Free(zero_pos);
                R_Free(observed);
                error("Too many simplifying assumptions to enumerate.");
            }

            if (keys_count + combinations > keys_capacity) {
                size_t new_capacity = keys_capacity;
                while (keys_count + combinations > new_capacity) {
                    if (new_capacity > SIZE_MAX / 2) {
                        R_Free(keys);
                        R_Free(zero_pos);
                        R_Free(observed);
                        error("Too many simplifying assumptions to enumerate.");
                    }
                    new_capacity *= 2;
                }

                keys = (int *) R_Realloc(keys, new_capacity, int);
                keys_capacity = new_capacity;
            }

            for (size_t comb = 0; comb < combinations; comb++) {
                size_t rest = comb;
                int key = base;

                for (int z = 0; z < zero_count; z++) {
                    int c = zero_pos[z];
                    int level = (int) (rest % (size_t) p_noflevels[c]);
                    rest /= (size_t) p_noflevels[c];
                    key += p_mbaseexpr[c] * level;
                }

                keys[keys_count] = key;
                keys_count++;
            }
        }

        qsort(keys, keys_count, sizeof(int), qca_compare_int);

        size_t kept_count = 0;
        int previous = 0;
        Rboolean have_previous = FALSE;

        for (size_t i = 0; i < keys_count; i++) {
            int key = keys[i];

            if (have_previous && key == previous) {
                continue;
            }
            have_previous = TRUE;
            previous = key;

            if (!qca_int_in_sorted(key, observed, input_rows)) {
                keys[kept_count] = key;
                kept_count++;
            }
        }

        if (kept_count > INT_MAX) {
            R_Free(keys);
            R_Free(zero_pos);
            R_Free(observed);
            error("Too many simplifying assumptions to return.");
        }

        int kept_rows = (int) kept_count;
        SEXP matrix;
        SET_VECTOR_ELT(usage, 6, matrix = allocMatrix(INTSXP, kept_rows, input_cols));
        int *p_matrix = INTEGER(matrix);

        for (int r = 0; r < kept_rows; r++) {
            int row_name = 0;

            for (int c = 0; c < input_cols; c++) {
                int value = ((keys[r] / p_mbaseplus[c]) % (p_noflevels[c] + 1)) - 1;
                p_matrix[c * kept_rows + r] = value;
                row_name += value * p_mbase[c];
            }

            keys[r] = row_name + 1;
        }

        SEXP dimnames = PROTECT(allocVector(VECSXP, 2));
        if (kept_rows > 0) {
            SEXP rownames = PROTECT(allocVector(STRSXP, kept_rows));
            for (int r = 0; r < kept_rows; r++) {
                char buffer[32];
                snprintf(buffer, sizeof(buffer), "%d", keys[r]);
                SET_STRING_ELT(rownames, r, mkChar(buffer));
            }
            SET_VECTOR_ELT(dimnames, 0, rownames);
            UNPROTECT(1);
        }
        if (!Rf_isNull(input_colnames)) {
            SET_VECTOR_ELT(dimnames, 1, input_colnames);
        }
        setAttrib(matrix, R_DimNamesSymbol, dimnames);
        SET_VECTOR_ELT(out, s, matrix);

        UNPROTECT(1);
        R_Free(keys);
        R_Free(zero_pos);
    }

    R_Free(observed);
    UNPROTECT(2);
    return(out);
}

SEXP C_createMatrix(SEXP input) {
    PROTECT(input);
    SEXP matrix, noflevels, arrange, depth;


    SEXP usage = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(usage, 0, noflevels = coerceVector(VECTOR_ELT(input, 0), INTSXP));
    SET_VECTOR_ELT(usage, 1, arrange   = coerceVector(VECTOR_ELT(input, 1), INTSXP));
    SET_VECTOR_ELT(usage, 2, depth   = coerceVector(VECTOR_ELT(input, 2), INTSXP));
    int *p_noflevels = INTEGER(noflevels);
    int *p_arrange = INTEGER(arrange);
    int *p_depth = INTEGER(depth);


    int ncols = length(noflevels);

    int nofl[ncols];
    for (int c = 0; c < ncols; c++) {
        nofl[c] = p_noflevels[c];
    }

    if (p_depth[0] > ncols) {
        p_depth[0] = ncols;
    }

    int intarrange = p_arrange[0];
    int intdepth = p_depth[0];

    int nrows;
    calculate_rows(&nrows, ncols, nofl, intarrange, intdepth);

    SET_VECTOR_ELT(usage, 3, matrix = allocMatrix(INTSXP, nrows, ncols));


    generate_matrix(nrows, ncols, nofl, intarrange, intdepth, INTEGER(matrix));

    if (length(input) > 3) { // the colnms is supplied
        SEXP dimnames;
        SET_VECTOR_ELT(usage, 3, dimnames = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(input, 3)); // transposed!
        setAttrib(matrix, R_DimNamesSymbol, dimnames);
    }

    UNPROTECT(2);
    return(matrix);
}

SEXP C_superSubset(SEXP x, SEXP noflevels, SEXP fuz, SEXP vo,
                    SEXP nec, SEXP inclcut, SEXP covcut, SEXP depth) {

    SEXP usage = PROTECT(allocVector(VECSXP, 19));
    SET_VECTOR_ELT(usage,  0, x         = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage,  1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage,  2, fuz       = coerceVector(fuz, INTSXP));
    SET_VECTOR_ELT(usage,  3, vo        = coerceVector(vo, REALSXP));
    SET_VECTOR_ELT(usage,  4, nec       = coerceVector(nec, INTSXP));
    SET_VECTOR_ELT(usage,  5, inclcut   = coerceVector(inclcut, REALSXP));
    SET_VECTOR_ELT(usage,  6, covcut    = coerceVector(covcut, REALSXP));
    SET_VECTOR_ELT(usage,  7, depth     = coerceVector(depth, INTSXP));

    double *p_x = REAL(x);
    int *p_noflevels = INTEGER(noflevels);
    int *p_fuz = INTEGER(fuz);
    double *p_vo = REAL(vo);
    int *p_nec = INTEGER(nec);
    double *p_inclcut = REAL(inclcut);
    double *p_covcut = REAL(covcut);
    int *p_depth = INTEGER(depth);

    int estim1 = 1000;
    int estim2 = 1000; // rows for each of the solution types (conjunctions and disjunctions)
    int xrows = nrows(x);
    int xcols = ncols(x);
    int nconds = ncols(x); // redundant, to preserve some code similarity

    if (p_depth[0] == 0) {
        p_depth[0] = nconds;
    }

    SEXP tmconj, tmdisj, ticpr_conj, ticpr_disj, combkl, tcoms_conj, tcoms_disj,
         indx_conj, indx_disj, ck_conj, ck_disj;
    // Another idea is _not_ to store the combinations, but instead produce the rownames
    // for the icp matrices 1 and 2, using the dimnames attribute
    SET_VECTOR_ELT(usage,  8, tmconj     = allocVector(INTSXP, nconds * estim1));
    SET_VECTOR_ELT(usage,  9, tmdisj     = allocVector(INTSXP, nconds * estim2));
    SET_VECTOR_ELT(usage, 10, ticpr_conj = allocVector(REALSXP, 3 * estim1));    // temporary incl cov pri ron for conjunctions
    SET_VECTOR_ELT(usage, 11, ticpr_disj = allocVector(REALSXP, 3 * estim2));    // temporary incl cov pri ron for disjunctions
    SET_VECTOR_ELT(usage, 12, tcoms_conj = allocMatrix(REALSXP, xrows, estim1)); // temporary component membership scores for conjunctions
    SET_VECTOR_ELT(usage, 13, tcoms_disj = allocMatrix(REALSXP, xrows, estim1)); // temporary component membership scores for disjunctions
    SET_VECTOR_ELT(usage, 14, indx_conj  = allocVector(INTSXP, p_depth[0] * estim1)); // combinations of c out of k for conjunctions
    SET_VECTOR_ELT(usage, 15, indx_disj  = allocVector(INTSXP, p_depth[0] * estim2)); // combinations of c out of k for disjunctions
    SET_VECTOR_ELT(usage, 16, ck_conj    = allocVector(INTSXP, estim1)); // how many conditions in a conjunction
    SET_VECTOR_ELT(usage, 17, ck_disj    = allocVector(INTSXP, estim2)); // how many conditions in a disjunction
    int    *p_tmconj     = INTEGER(tmconj);
    int    *p_tmdisj     = INTEGER(tmdisj);
    double *p_ticpr_conj = REAL(ticpr_conj);
    double *p_ticpr_disj = REAL(ticpr_disj);
    double *p_tcoms_conj = REAL(tcoms_conj);
    double *p_tcoms_disj = REAL(tcoms_disj);
    int    *p_indx_conj  = INTEGER(indx_conj);
    int    *p_indx_disj  = INTEGER(indx_disj);
    int    *p_ck_conj    = INTEGER(ck_conj);
    int    *p_ck_disj    = INTEGER(ck_disj);
    double minx[xrows], maxx[xrows];
    double incovpron[6];
    double so = 0.0,
           sum_minx,
           sum_maxx,
           sum_1_minx,
           sum_1_maxx,
           sum_1_min_y_minx,
           sum_1_min_y_maxx,
           sum_min_y_minx,
           sum_min_y_maxx,
           prisum_minx,
           prisum_maxx,
           tmpv11, tmpv12, tmpv21, tmpv22;

    int found1 = 0;
    int found = 0;
    int foundk1 = 0;
    int foundk2 = 0;

    if (nconds < p_depth[0]) {
        p_depth[0] = nconds;
    }


    // sum of the outcome variable
    for (int i = 0; i < length(vo); i++) {
        so += p_vo[i];
    }

    int chkred[nconds], inclcov;

    int k = 1;
    int foundk = 1;
    while (k <= p_depth[0] && foundk) {

        if (found1 + found > 0 && k > 3) {
            /*
            k > 3 is arbitrary, there are situations when a
            necessary conjunction with a single condition doesn't exclude
            disjunctions with more than 2 conditions, and normally after
            k = 1 it would stop
            */
            foundk = 0;
        }

        int klnofl[k];

        // calculate all combinations of k conditions out of nconds
        unsigned long long int maxtasks = nchoosek(nconds, k);

        for (unsigned long long int task = 0; task < maxtasks; task++) {
            if (task > 0 && task % 1024 == 0) {
                R_CheckUserInterrupt();
            }

            int tempk[k];
            unsigned long long int combination = task;
            int x = 0;

            for (int i = 0; i < k; i++) {
                while (1) {
                    unsigned long long int cval = nchoosek(nconds - (x + 1), k - (i + 1));
                    if (cval == 0 || cval > combination) {
                        break;
                    }
                    combination -= cval;
                    x++;
                }

                if (x < 0) {
                    x = 0;
                }
                if (x >= nconds) {
                    x = nconds - 1;
                }

                tempk[i] = x;
                x++;
            }

            int klcols[k];
            int klrows = 1;
            for (int j = 0; j < k; j++) {
                klnofl[j] = p_noflevels[tempk[j]];
                klrows *= klnofl[j];
                klcols[j] = j;
            }

            SET_VECTOR_ELT(usage, 18, combkl = allocVector(INTSXP, klrows * k));
            int *p_combkl = INTEGER(combkl);
            fill_matrix(klrows, k, klnofl, p_combkl, 0, klcols, 0);

            for (int kli = 0; kli < klrows; kli++) {

                for (int c = 0; c < nconds; c++) {
                    chkred[c] = 0;
                }

                for (int j = 0; j < k; j++) {
                    chkred[tempk[j]] = p_combkl[j * klrows + kli] + 1; // + 1 to get R 1-based numbers
                }

                // chkred is now an implicant matrix line

                sum_minx = 0;         // sum(min(x[i,]))
                sum_maxx = 0;         // sum(max(x[i,]))
                sum_min_y_minx = 0;   // sum(min(min(x[i,]), y))
                sum_min_y_maxx = 0;   // sum(min(max(x[i,]), y))
                prisum_minx = 0;      // PRI for conjunctions
                prisum_maxx = 0;      // PRI for disjunctions
                sum_1_minx = 0;       // sum(1 - min(x[i,]))
                sum_1_min_y_minx = 0; // sum(1 - min(y, min(x[i,])))
                sum_1_maxx = 0;       // sum(1 - max(x[i,]))
                sum_1_min_y_maxx = 0; // sum(1 - min(y, max(x[i,])))

                // =======================
                // this part is similar to the truth table calculations, maybe a separate function
                //
                for (int r = 0; r < xrows; r++) { // loop over every line of the data matrix

                    double min_local = 1000000;
                    double max_local = 0;

                    for (int c = 0; c < xcols; c++) { // loop over each column of the data matrix
                        double value = p_x[c * xrows + r];

                        if (p_fuz[c]) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
                            if (chkred[c] == 1) {
                                value = 1 - value;
                            }
                        }
                        else {
                            if (chkred[c] == (value + 1)) {
                                value = 1; // consistent
                            }
                            else {
                                value = 0; // inconsistent
                            }
                        }

                        if (chkred[c] != 0) {
                            if (value < min_local) {
                                min_local = value;
                            }

                            if (value > max_local) {
                                max_local = value;
                            }
                        }

                    } // end of loop over data columns

                    minx[r] = min_local;   // min(x[i, ])
                    maxx[r] = max_local;   // max(x[i, ])

                    sum_minx += min_local; // sum(min(x[i, ]))
                    sum_maxx += max_local; // sum(max(x[i, ]))
                    sum_min_y_minx += (min_local < p_vo[r]) ? min_local : p_vo[r];
                    sum_min_y_maxx += (max_local < p_vo[r]) ? max_local : p_vo[r];

                    if (p_nec[0]) {  // necessity, calculate RoN
                        sum_1_minx += 1 - min_local;
                        sum_1_maxx += 1 - max_local;
                        sum_1_min_y_minx += 1 - ((min_local < p_vo[r]) ? min_local : p_vo[r]);
                        sum_1_min_y_maxx += 1 - ((max_local < p_vo[r]) ? max_local : p_vo[r]);

                    }
                    else {           // sufficiency, calculate PRI
                        double tmpv11_local = (min_local < p_vo[r]) ? min_local : p_vo[r];
                        double tmpv12_local = 1 - p_vo[r];
                        prisum_minx += (tmpv11_local < tmpv12_local) ? tmpv11_local : tmpv12_local;
                        double tmpv21_local = (max_local < p_vo[r]) ? max_local : p_vo[r];
                        double tmpv22_local = 1 - max_local;
                        prisum_maxx += (tmpv21_local < tmpv22_local) ? tmpv21_local : tmpv22_local;
                    }


    /*
    if (k == 1) {
        Rprintf("       l: %d; min: %f; max: %f; t11: %f; t12: %f; t21: %f; t22: %f; vo: %f\n", i, min, max, tmpv11, tmpv12, tmpv21, tmpv22, p_vo[i]);
    }
    */

                } // end of loop over data rows

                //
                // =======================


                // Calculate inclusion, coverage, PRI (or RoN)
                /*
                0. incl suf / cov nec \
                1. cov suf / incl nec / for conjunctions
                2. incl nec \
                3. cov nec  / for disjunctions
                4. PRI / RoN  for conjunctions
                5. PRI / RoN  for disjunctions
                */

                incovpron[0] = (sum_min_y_minx == 0 && sum_minx == 0)?0:(sum_min_y_minx/sum_minx);
                incovpron[1] = (sum_min_y_minx == 0 && so == 0)?0:(sum_min_y_minx/so);
                incovpron[2] = (sum_min_y_maxx == 0 && so == 0)?0:(sum_min_y_maxx/so);
                incovpron[3] = (sum_min_y_maxx == 0 && sum_maxx == 0)?0:(sum_min_y_maxx/sum_maxx);
                if (p_nec[0]) {
                    incovpron[4] = (sum_1_minx == 0 && sum_1_min_y_minx == 0)?0:(sum_1_minx/sum_1_min_y_minx);
                    incovpron[5] = (sum_1_maxx == 0 && sum_1_min_y_maxx == 0)?0:(sum_1_maxx/sum_1_min_y_maxx);
                }
                else {
                    tmpv11 = sum_min_y_minx - prisum_minx;
                    tmpv12 = (p_nec[0]?so:sum_minx) - prisum_minx;
                    incovpron[4] = (tmpv11 == 0 && tmpv12 == 0)?0:(tmpv11/tmpv12);

                    tmpv21 = sum_min_y_maxx - prisum_maxx;
                    tmpv22 = so - prisum_maxx;
                    incovpron[5] = (tmpv21 == 0 && tmpv22 == 0)?0:(tmpv21/tmpv22);
                }

                // inclcov for conjunctions
                inclcov = incovpron[p_nec[0]] >= p_inclcut[0] && incovpron[1 - p_nec[0]] >= p_covcut[0];

    /*
    if (k == 1) {
        for (int c = 0; c < nconds; c++) {
            Rprintf("%d ", chkred[c]);
        }
    }


    if (k == 1) {
        Rprintf("         ");
        for (int i = 0; i < 6; i++) {
            Rprintf("%f ", incovpron[i]);
        }
        Rprintf("\n");
    }


    if (k == 1) {
        Rprintf("       spmin: %f; sminx: %f; spmax: %f; smaxx: %f; so: %f; t11: %f; t12: %f; t21: %f; t22: %f; psmin: %f; psmax: %f\n",
                        sum_min_y_minx, sum_minx, sum_min_y_maxx, sum_maxx, so, tmpv11, tmpv12, tmpv21, tmpv22, prisum_minx, prisum_maxx);
    }

    for (int c = 0; c < nconds; c++) {
        Rprintf("%d ", chkred[c]);
    }
    Rprintf(";  ");
    for (int c = 0; c < nconds; c++) {
        Rprintf("%f ", copyline[c]);
    }
    Rprintf("  ;   %d   ;   ", inclcov);
    for (int c = 0; c < 6; c++) {
        Rprintf("%f ", incovpron[c]);
    }
    Rprintf("\n");
    */

                int redundant = 0;

                if (inclcov) {
                    if (foundk1 > 0 && !p_nec[0]) { // check if redundant, only for sufficiency
                        int i = 0;
                        while (i < foundk1 && !redundant) {
                            // int sumnonz = 0; // nu mai este necesar daca merge
                            int sumeq = 0;
                            int v = 0;

                            while (sumeq == v && v < p_ck_conj[i]) {
                                for (int c = 0; c < k; c++) {
                                    if (p_indx_conj[i * p_depth[0] + v] == tempk[c] + 1) {
                                        sumeq += (p_tmconj[i * nconds + p_indx_conj[i * p_depth[0] + v] - 1] == chkred[tempk[c]]);
                                    }
                                }
                                v += 1;
                            }

                            if (sumeq == v) {
                                redundant = 1;
                            }

                            /*
                            for (int c = 0; c < nconds; c++) {

                                // tmconj is a horizontal matrix: i * conds + c
                                if (p_tmconj[i * nconds + c] != 0) {
                                    sumnonz += 1;
                                    if (chkred[c] == p_tmconj[i * nconds + c]) {
                                        sumeq += 1;
                                    }
                                }
                            }

                            if (sumnonz == sumeq) {
                                redundant = 1;
                            }
                            */
                            i += 1;
                        }
                    }

                    if (!redundant) { // we already know inclcov = 1
                        for (int c = 0; c < nconds; c++) {
                            p_tmconj[found1 * nconds + c] = chkred[c];
                        }

                        // ADD inclusion, coverage and PRI into ticpr_conj
                        p_ticpr_conj[found1 * 3 + 0] = incovpron[p_nec[0]];
                        p_ticpr_conj[found1 * 3 + 1] = incovpron[4];
                        p_ticpr_conj[found1 * 3 + 2] = incovpron[1 - p_nec[0]];

                        for (int r = 0; r < xrows; r++) {
                            p_tcoms_conj[found1 * xrows + r] = minx[r];
                        }

                        for (int c = 0; c < k; c++) {
                            p_indx_conj[p_depth[0] * found1 + c] = tempk[c] + 1;
                        }

                        p_ck_conj[found1] = k;

                        foundk += 1;
                        found1 += 1;

                        // double the matrices if needed
                        if (found1 == estim1) {
                            int copytm_conj[nconds * found1];
                            int tindx_conj[p_depth[0] * found1];
                            int tck_conj[found1];
                            double copyticpr_conj[3 * found1];
                            double copytcoms_conj[xrows * found1];
                            for (int i = 0; i < nconds * found1; i++) {
                                copytm_conj[i] = p_tmconj[i];
                            }

                            for (int i = 0; i < 3 * found1; i++) {
                                copyticpr_conj[i] = p_ticpr_conj[i];
                            }

                            for (int i = 0; i < xrows * found1; i++) {
                                copytcoms_conj[i] = p_tcoms_conj[i];
                            }

                            for (int i = 0; i < p_depth[0] * found1; i++) {
                                tindx_conj[i] = p_indx_conj[i];
                            }

                            for (int i = 0; i < found1; i++) {
                                tck_conj[i] = p_ck_conj[i];
                            }

                            estim1 *= 2;
                            SET_VECTOR_ELT(usage, 8,  tmconj     = allocVector(INTSXP, nconds * estim1));
                            p_tmconj = INTEGER(tmconj);
                            SET_VECTOR_ELT(usage, 10, ticpr_conj = allocVector(REALSXP, 3 * estim1));
                            p_ticpr_conj = REAL(ticpr_conj);
                            SET_VECTOR_ELT(usage, 12, tcoms_conj = allocMatrix(REALSXP, xrows, estim1));
                            p_tcoms_conj = REAL(tcoms_conj);
                            SET_VECTOR_ELT(usage, 14, indx_conj  = allocVector(INTSXP, p_depth[0] * estim1));
                            p_indx_conj = INTEGER(indx_conj);
                            SET_VECTOR_ELT(usage, 16, ck_conj    = allocVector(INTSXP, estim1));
                            p_ck_conj = INTEGER(ck_conj);

                            for (int i = 0; i < nconds * found1; i++) {
                                p_tmconj[i] = copytm_conj[i];
                            }

                            for (int i = 0; i < 3 * found1; i++) {
                                p_ticpr_conj[i] = copyticpr_conj[i];
                            }

                            for (int i = 0; i < xrows * found1; i++) {
                                p_tcoms_conj[i] = copytcoms_conj[i];
                            }

                            for (int i = 0; i < p_depth[0] * found1; i++) {
                                p_indx_conj[i] = tindx_conj[i];
                            }

                            for (int i = 0; i < found1; i++) {
                                p_ck_conj[i] = tck_conj[i];
                            }

                        }
                    }
                }
                else {
                    if (p_nec[0]) {
                        // inclcov for necessary disjunctions
                        inclcov = incovpron[2] >= p_inclcut[0] && incovpron[3] >= p_covcut[0];

                        redundant = 0;

                        if (inclcov && foundk1 > 0) {
                            int i = 0;
                            while (i < foundk1 && !redundant) {
                                // int sumnonz = 0;
                                int sumeq = 0;
                                int v = 0;

                                while (sumeq == v && v < p_ck_conj[i]) {
                                    for (int c = 0; c < k; c++) {
                                        if (p_indx_conj[i * p_depth[0] + v] == tempk[c] + 1) {
                                            sumeq += (p_tmconj[i * nconds + p_indx_conj[i * p_depth[0] + v] - 1] == chkred[tempk[c]]);
                                        }
                                    }
                                    v += 1;
                                }

                                if (sumeq == v) {
                                    redundant = 1;
                                }

                                /*
                                for (int c = 0; c < nconds; c++) {

                                    // tmconj is also a horizontal matrix: i * conds + c
                                    if (p_tmconj[i * nconds + c] != 0) {
                                        sumnonz += 1;
                                        if (chkred[c] == p_tmconj[i * nconds + c]) {
                                            sumeq += 1;
                                        }
                                    }
                                }

                                if (sumnonz == sumeq) {
                                    redundant = 1;
                                }
                                */

                                i += 1;
                            }
                        }

                        if (inclcov && foundk2 > 0 && !redundant) {
                            int i = 0;
                            while (i < foundk2 && !redundant) {
                                // int sumnonz = 0;
                                int sumeq = 0;
                                int v = 0;

                                while (sumeq == v && v < p_ck_disj[i]) {
                                    for (int c = 0; c < k; c++) {
                                        if (p_indx_disj[i * p_depth[0] + v] == tempk[c] + 1) {
                                            sumeq += (p_tmdisj[i * nconds + p_indx_disj[i * p_depth[0] + v] - 1] == chkred[tempk[c]]);
                                        }
                                    }
                                    v += 1;
                                }

                                if (sumeq == v) {
                                    redundant = 1;
                                }

                                /*
                                for (int c = 0; c < nconds; c++) {

                                    // tmconj is also a horizontal matrix: i * conds + c
                                    if (p_tmdisj[i * nconds + c] != 0) {
                                        sumnonz += 1;
                                        if (chkred[c] == p_tmdisj[i * nconds + c]) {
                                            sumeq += 1;
                                        }
                                    }
                                }

                                if (sumnonz == sumeq) {
                                    redundant = 1;
                                }
                                */

                                i += 1;
                            }
                        }

                        if (inclcov && !redundant) {
                            for (int c = 0; c < nconds; c++) {
                                p_tmdisj[found * nconds + c] = chkred[c];
                            }

                            // ADD inclusion, coverage and PRI into ticpr_disj
                            p_ticpr_disj[found * 3 + 0] = incovpron[2];
                            p_ticpr_disj[found * 3 + 1] = incovpron[5];
                            p_ticpr_disj[found * 3 + 2] = incovpron[3];

                            for (int r = 0; r < xrows; r++) {
                                p_tcoms_disj[found * xrows + r] = maxx[r];
                            }

                            p_ck_disj[found] = k;

                            for (int c = 0; c < k; c++) {
                                p_indx_disj[p_depth[0] * found + c] = tempk[c] + 1;
                            }

                            foundk += 1;
                            found += 1;

                            // double the matrices if needed
                            if (found == estim2) {
                                int copytmdisj[nconds * found];
                                int tindx_disj[p_depth[0] * found];
                                int tck_disj[found];
                                double copyticpr_disj[3 * found];
                                double copytcoms_disj[xrows * found];

                                for (int i = 0; i < nconds * found; i++) {
                                    copytmdisj[i] = p_tmdisj[i];
                                }

                                for (int i = 0; i < 3 * found; i++) {
                                    copyticpr_disj[i] = p_ticpr_disj[i];
                                }

                                for (int i = 0; i < xrows * found; i++) {
                                    copytcoms_disj[i] = p_tcoms_disj[i];
                                }

                                for (int i = 0; i < p_depth[0] * found; i++) {
                                    tindx_disj[i] = p_indx_disj[i];
                                }

                                for (int i = 0; i < found; i++) {
                                    tck_disj[i] = p_ck_disj[i];
                                }

                                estim2 *= 2;
                                SET_VECTOR_ELT(usage,  9, tmdisj     = allocVector(INTSXP, nconds * estim2));
                                p_tmdisj = INTEGER(tmdisj);
                                SET_VECTOR_ELT(usage, 11, ticpr_disj = allocVector(REALSXP, 3 * estim2));
                                p_ticpr_disj = REAL(ticpr_disj);
                                SET_VECTOR_ELT(usage, 13, tcoms_disj = allocMatrix(REALSXP, xrows, estim2));
                                p_tcoms_disj = REAL(tcoms_disj);
                                SET_VECTOR_ELT(usage, 15, indx_disj  = allocVector(REALSXP, p_depth[0] * estim2));
                                p_indx_disj = INTEGER(indx_disj);
                                SET_VECTOR_ELT(usage, 17, ck_disj    = allocVector(INTSXP, estim2));
                                p_ck_disj = INTEGER(ck_disj);

                                for (int i = 0; i < nconds * found; i++) {
                                    p_tmdisj[i] = copytmdisj[i];
                                }

                                for (int i = 0; i < 3 * found; i++) {
                                    p_ticpr_disj[i] = copyticpr_disj[i];
                                }

                                for (int i = 0; i < xrows * found; i++) {
                                    p_tcoms_disj[i] = copytcoms_disj[i];
                                }

                                for (int i = 0; i < p_depth[0] * found; i++) {
                                    p_indx_disj[i] = tindx_disj[i];
                                }

                                for (int i = 0; i < found; i++) {
                                    p_ck_disj[i] = tck_disj[i];
                                }
                            }
                        }
                    }
                } // end verifying if the row passes the cutoffs, for both conjunctions and disjunctions
            } // next kli from klrows
        }

        foundk1 = found1;
        foundk2 = found;
        k++;
    }



    SEXP icpr_conj, icpr_disj, mconj, mdisj, coms_conj, coms_disj;
    SEXP result = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(result, 0, icpr_conj = allocMatrix(REALSXP, found1, 3)); // incl cov pri ron for conjunctions
    SET_VECTOR_ELT(result, 1, icpr_disj = allocMatrix(REALSXP, found, 3)); // incl cov pri ron for disjunctions
    SET_VECTOR_ELT(result, 2, mconj = allocMatrix(INTSXP, found1, nconds)); // base matrix for conjunctions
    SET_VECTOR_ELT(result, 3, mdisj = allocMatrix(INTSXP, found, nconds)); // base matrix for disjunctions
    SET_VECTOR_ELT(result, 4, coms_conj = allocMatrix(REALSXP, xrows, found1)); // component membership scores for conjunctions
    SET_VECTOR_ELT(result, 5, coms_disj = allocMatrix(REALSXP, xrows, found)); // component membership scores for disjunctions
    double *p_icpr_conj = REAL(icpr_conj);
    double *p_icpr_disj = REAL(icpr_disj);
    int *p_mconj = INTEGER(mconj);
    int *p_mdisj = INTEGER(mdisj);
    double *p_coms_conj = REAL(coms_conj);
    double *p_coms_disj = REAL(coms_disj);


    for (int r = 0; r < found1; r++) { // transpose from a horizontal to a vertical matrix
        for (int c = 0; c < 3; c++) {
            p_icpr_conj[c * found1 + r] = p_ticpr_conj[r * 3 + c];
        }
    }

    for (int r = 0; r < found1; r++) { // transpose from a horizontal to a vertical matrix
        for (int c = 0; c < nconds; c++) {
            p_mconj[c * found1 + r] = p_tmconj[r * nconds + c];
        }
    }

    for (int r = 0; r < found; r++) { // transpose from a horizontal to a vertical matrix
        for (int c = 0; c < 3; c++) {
            p_icpr_disj[c * found + r] = p_ticpr_disj[r * 3 + c];
        }
    }

    for (int r = 0; r < found; r++) { // transpose from a horizontal to a vertical matrix
        for (int c = 0; c < nconds; c++) {
            p_mdisj[c * found + r] = p_tmdisj[r * nconds + c];
        }
    }

    for (int i = 0; i < xrows * found1; i++) {
        p_coms_conj[i] = p_tcoms_conj[i];
    }

    for (int i = 0; i < xrows * found; i++) {
        p_coms_disj[i] = p_tcoms_disj[i];
    }

    UNPROTECT(2);

    return(result);
}

SEXP C_QMC(SEXP tt, SEXP noflevels) {
    // because this is the classical Quine-McCluskey, all configurations are positive
    // there are no negative configurations to worry about, therefore
    // the number of causal conditions is equal to the number of columns in the truth table

    SEXP tempmat, copymat, order, cl; // cl = complexity level
    int *p_tt, *p_noflevels, *p_pimat, *p_tempmat, *p_minimized, *p_copymat, *p_order,  *p_cl;
    SEXP usage = PROTECT(allocVector(VECSXP, 10));

    SEXP dimnames, colnms;
    SET_VECTOR_ELT(usage, 8, dimnames = getAttrib(tt, R_DimNamesSymbol));
    // fix CRAN warning colnms may be used uninitialised
    // AND protect it at the same time
    SET_VECTOR_ELT(usage, 9, colnms = getAttrib(tt, R_DimNamesSymbol));

    if (!Rf_isNull(dimnames)) {
        colnms = VECTOR_ELT(getAttrib(tt, R_DimNamesSymbol), 1);
    }

    SET_VECTOR_ELT(usage, 0, tt = coerceVector(tt, INTSXP));
    p_tt = INTEGER(tt);

    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    p_noflevels = INTEGER(noflevels);

    int nimplicants = nrows(tt);
    int nconds = ncols(tt);

    SET_VECTOR_ELT(usage, 2, tempmat = allocMatrix(INTSXP, nconds, nimplicants));
    p_tempmat = INTEGER(tempmat);

    for (int r = 0; r < nimplicants; r++) {
        for (int c = 0; c < nconds; c++) {
            p_tempmat[r * nconds + c] = p_tt[c * nimplicants + r];  // initiate with the transposed truth table
        }
    }


    int found = 1;

    // int iteration = 0;

    while (found > 0 && nimplicants > 1) {

        // while (blabla < 3) {

        // iteration++;
        // Rprintf("iteration: %d; found: %d\n", iteration, found);

        // reset the counter to zero for the next iteration
        found = 0;

        p_minimized = (int *) calloc((size_t) nimplicants, sizeof(int));
        int estimpi = 10000;
        p_pimat = (int *) calloc((size_t) nconds * (size_t) estimpi, sizeof(int));
        Rboolean qmc_alloc_failed = false;

        if (p_minimized == NULL || p_pimat == NULL) {
            free(p_minimized);
            free(p_pimat);
            error("Memory allocation failed during QMC initialization.");
        }

        unsigned long long int maxtasks = nchoosek(nimplicants, 2);

        for (unsigned long long int task = 0; task < maxtasks; task++) {
            if (task > 0 && task % 1024 == 0) {
                R_CheckUserInterrupt();
            }

            int combs_local[2];
            unsigned long long int combination = task;
            int x = 0;

            for (int i = 0; i < 2; i++) {
                while (1) {
                    unsigned long long int cval = nchoosek(nimplicants - (x + 1), 2 - (i + 1));
                    if (cval == 0 || cval > combination) {
                        break;
                    }
                    combination -= cval;
                    x++;
                }

                if (x < 0) {
                    x = 0;
                }
                if (x >= nimplicants) {
                    x = nimplicants - 1;
                }

                combs_local[i] = x;
                x++;
            }

            int temp_local[nconds];
            int r = 0;
            int diffs = 0;
            int which = 0;
            Rboolean comparable = TRUE;
            while (diffs < 2 && r < nconds && comparable) {
                temp_local[r] = p_tempmat[combs_local[0] * nconds + r];
                if (temp_local[r] != p_tempmat[combs_local[1] * nconds + r]) {
                    comparable = temp_local[r] > 0 && p_tempmat[combs_local[1] * nconds + r] > 0;
                    diffs++;
                    which = r;
                    temp_local[r] = 0;
                }
                r++;
            }

            if (diffs == 1 && comparable) {
                int minrows[p_noflevels[which]];
                minrows[0] = combs_local[0];
                minrows[1] = combs_local[1];
                int tominimize = 2;
                int c = 0;

                while (c < nimplicants && tominimize < p_noflevels[which]) {
                    if (c != combs_local[0] && c != combs_local[1]) {
                        Rboolean equal = TRUE;
                        int rr = 0;
                        while (rr < nconds && equal) {
                            if (rr != which) {
                                equal = temp_local[rr] == p_tempmat[c * nconds + rr];
                            }
                            rr++;
                        }

                        if (equal) {
                            minrows[tominimize] = c;
                            tominimize++;
                        }
                    }

                    c++;
                }

                if (tominimize == p_noflevels[which]) {
                    {
                        if (!qmc_alloc_failed) {
                            for (int i = 0; i < tominimize; i++) {
                                p_minimized[minrows[i]] = TRUE;
                            }

                            int f = 0;
                            Rboolean Runique = TRUE;
                            while (f < found && Runique) {
                                Rboolean equal = TRUE;
                                int rr = 0;
                                while (rr < nconds && equal) {
                                    equal = temp_local[rr] == p_pimat[f * nconds + rr];
                                    rr++;
                                }

                                Runique = !equal;
                                f++;
                            }

                            if (Runique) {
                                if (found == estimpi) {
                                    int new_estim = estimpi * 2;
                                    int *tmp = (int *) realloc(
                                        p_pimat,
                                        (size_t) nconds * (size_t) new_estim * sizeof(int)
                                    );

                                    if (tmp == NULL) {
                                        qmc_alloc_failed = true;
                                    } else {
                                        p_pimat = tmp;
                                        estimpi = new_estim;
                                    }
                                }

                                if (!qmc_alloc_failed) {
                                    for (int rr = 0; rr < nconds; rr++) {
                                        p_pimat[found * nconds + rr] = temp_local[rr];
                                    }
                                    found++;
                                }
                            }
                        }
                    }
                }
            }
        }

        R_CheckUserInterrupt();
        if (qmc_alloc_failed) {
            free(p_minimized);
            free(p_pimat);
            error("Memory allocation failed during QMC resize.");
        }

        int nonmin = 0;
        for (int i = 0; i < nimplicants; i++) {
            nonmin += !p_minimized[i];
        }

        SET_VECTOR_ELT(usage, 5, copymat = allocVector(INTSXP, (found + nonmin) * nconds));
        p_copymat = INTEGER(copymat);

        int foundlent = found * nconds;
        for (int i = 0; i < foundlent; i++) {
            p_copymat[i] = p_pimat[i];
        }

        if (nonmin > 0) {
            for (int i = 0; i < nimplicants; i++) {
                if (!p_minimized[i]) {
                    for (int r = 0; r < nconds; r++) {
                        p_copymat[foundlent + r] = p_tempmat[i * nconds + r];
                    }
                    foundlent += nconds;
                }
            }
        }

        SET_VECTOR_ELT(usage, 2, tempmat = allocMatrix(INTSXP, nconds, found + nonmin));
        p_tempmat = INTEGER(tempmat);

        for (int i = 0; i < foundlent; i++) {
            p_tempmat[i] = p_copymat[i];
        }
        free(p_minimized);
        free(p_pimat);

        //Rf_PrintValue(tempmat);

        nimplicants = ncols(tempmat);
    }



    SET_VECTOR_ELT(usage, 6, order = allocVector(INTSXP, nimplicants));
    p_order = INTEGER(order);

    SET_VECTOR_ELT(usage, 7, cl = allocVector(INTSXP, nimplicants));
    p_cl = INTEGER(cl);

    for (int c = 0; c < nimplicants; c++) {
        p_cl[c] = nconds;
        for (int r = 0; r < nconds; r++) {
            p_cl[c] -= p_tempmat[c * nconds + r] == 0 ? 1 : 0;
        }
    }

    sort_matrix(p_tempmat, p_order, p_cl, nconds, nimplicants);


    SET_VECTOR_ELT(usage, 5, copymat = allocMatrix(INTSXP, nimplicants, nconds));
    p_copymat = INTEGER(copymat);

    for (int c = 0; c < nconds; c++) {
        for (int r = 0; r < nimplicants; r++) {
            p_copymat[c * nimplicants + r] = p_tempmat[p_order[r] * nconds + c]; // transpose back to a normal shape
        }
    }


    if (!Rf_isNull(dimnames)) {
        SET_VECTOR_ELT(usage, 8, dimnames = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 1, colnms);

        setAttrib(copymat, R_DimNamesSymbol, dimnames);
    }

    UNPROTECT(1);
    return(copymat);
}

SEXP C_removeRedundants(SEXP rowno, SEXP noflevels, SEXP mbase) {

    int *pointer_next, *pointer_final, *pointer_temp1, *pointer_temp2, *pointer_rowno, *pointer_noflevels, *pointer_mbase;
    int previous, lmbase, ltemp2, lrowno, lmbasei, i, j, k, rn, finalength, lungime, flag2, flag1, templung;
    SEXP next, final, temp1, temp2;

    SEXP usage = PROTECT(allocVector(VECSXP, 7));

    SET_VECTOR_ELT(usage, 0, rowno = coerceVector(rowno, INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));

    pointer_rowno = INTEGER(rowno);
    pointer_noflevels = INTEGER(noflevels);
    pointer_mbase = INTEGER(mbase);
    lmbase = length(mbase);
    lrowno = length(rowno);


    // create the vector of next positions (similar to linked lists)
	SET_VECTOR_ELT(usage, 3, next = allocVector(INTSXP, lrowno));
	pointer_next = INTEGER(next);

	for (i = 0; i < lrowno; i++) {
	    pointer_next[i] = i + 1;
	}


	rn = 0;
	flag1 = 0;

	while (rn < lrowno) {
	    templung = 1;
	    previous = rn;

        SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, 1));
        pointer_temp1 = INTEGER(temp1);
        pointer_temp1[0] = pointer_rowno[rn];


        flag2 = 0;

        for (i = 0; i < lmbase; i++) {
            lmbasei = lmbase - i - 1;
            if (div(div(pointer_rowno[rn] - 1, pointer_mbase[lmbasei]).quot, pointer_noflevels[lmbasei] + 1).rem == 0) {
                flag2 = 1;
                lungime = templung * (pointer_noflevels[lmbasei] + 1);

                SET_VECTOR_ELT(usage, 5, temp2 = allocVector(INTSXP, lungime));
                pointer_temp2 = INTEGER(temp2);

                for (j = 0; j < length(temp1); j++) {
                    pointer_temp2[j] = pointer_temp1[j];
                    for (k = 0; k < pointer_noflevels[lmbasei]; k++) {
                        pointer_temp2[j + length(temp1)*(k + 1)] = pointer_temp1[j] + (k + 1)*pointer_mbase[lmbasei];
                    }
                }

                if (i < lmbase) {
                    SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, lungime));
                    pointer_temp1 = INTEGER(temp1);

                    for (j = 0; j < lungime; j++) {
                        pointer_temp1[j] = pointer_temp2[j];
                    }
                    templung = lungime;
                }
            }
        }

        if (flag2 == 1) { // the current row number has subsets, stored in temp2
            ltemp2 = length(temp2);

            i = pointer_next[previous];
            j = 0;
            while (i < lrowno && j < ltemp2) {
                if (pointer_rowno[i] < pointer_temp2[j]) {
                    previous = i;
                    i = pointer_next[i];
                }
                else if (pointer_rowno[i] > pointer_temp2[j]) {
                    j++;
                }
                else { // if (pointer_rowno[i] == pointer_temp2[j])
                    flag1 = 1;
                    pointer_next[previous] = pointer_next[i];
                    i = pointer_next[i];
                    j++;
                }
            }
        }

		rn = pointer_next[rn];
    }


    if (flag1 == 0) { // no single row has been found redundant
        UNPROTECT(1);
        return(rowno);
    }
    else {

        finalength = 0;
        i = 0;
        while (i < lrowno) {
            i = pointer_next[i];
            finalength++;
        }

        SET_VECTOR_ELT(usage, 6, final = allocVector(INTSXP, finalength));
        pointer_final = INTEGER(final);

        i = 0;
        j = 0;
        while (i < lrowno) {
            pointer_final[j] = pointer_rowno[i];
            i = pointer_next[i];
            j += 1;
        }

        UNPROTECT(1);
        return(final);
    }
}

SEXP C_findSubsets(SEXP rowno, SEXP noflevels, SEXP mbase, SEXP max) {
    int *prowno, *pnoflevels, *pmbase, *pmax, lmbase, lmbasei, i, j, k, lungime, flag, templung, *ptemp1, *ptemp2;
    SEXP temp1, temp2;



    SEXP usage = PROTECT(allocVector(VECSXP, 6));

    SET_VECTOR_ELT(usage, 0, rowno = coerceVector(rowno, INTSXP));
    SET_VECTOR_ELT(usage, 1, noflevels = coerceVector(noflevels, INTSXP));
    SET_VECTOR_ELT(usage, 2, mbase = coerceVector(mbase, INTSXP));

    prowno = INTEGER(rowno);
    pnoflevels = INTEGER(noflevels);
    pmbase = INTEGER(mbase);

    if (max == R_NilValue) {
        SET_VECTOR_ELT(usage, 3, max = allocVector(INTSXP, 1));
        pmax = INTEGER(max);
        pmax[0] = prowno[length(rowno) - 1];
    }
    else {
        SET_VECTOR_ELT(usage, 3, max = coerceVector(max, INTSXP));
        pmax = INTEGER(max);
    }

    SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, 1));
    ptemp1 = INTEGER(temp1);

    ptemp1[0] = prowno[0];
    flag = 0;
    lmbase = length(mbase);

    templung = 1;


    //Rprintf("length(mbase): %d\n", length(mbase));

    for (i = 0; i < lmbase; i++) {
        lmbasei = lmbase - i - 1;
        //Rprintf("rowno: %d; mbase.val: %d; noflevels.i: %d\n", prowno[0], pmbase[i], pnoflevels[i] + 1);
        //Rprintf("%d d/d %d = %d\n", prowno[0] - 1, pmbase[i], div(prowno[0] - 1, pmbase[i]).quot);
        //Rprintf("%d dd %d = %d\n", div(prowno[0] - 1, pmbase[i]).quot, pmbase[i], div(div(prowno[0] - 1, pmbase[i]).quot, pnoflevels[i] + 1).rem);
        if (div(div(prowno[0] - 1, pmbase[lmbasei]).quot, pnoflevels[lmbasei] + 1).rem == 0) {
            flag = 1;
            lungime = templung * (pnoflevels[lmbasei] + 1);
            //Rprintf("lungime: %d\n", lungime);
            SET_VECTOR_ELT(usage, 5, temp2 = allocVector(INTSXP, lungime));
            ptemp2 = INTEGER(temp2);

            for (j = 0; j < length(temp1); j++) {
                ptemp2[j] = ptemp1[j];
                for (k = 0; k < pnoflevels[lmbasei]; k++) {
                    ptemp2[j + length(temp1)*(k + 1)] = ptemp1[j] + (k + 1)*pmbase[lmbasei];
                }
            }

            if (i < length(mbase)) {
                SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, lungime));
                ptemp1 = INTEGER(temp1);

                for (j = 0; j < lungime; j++) {
                    ptemp1[j] = ptemp2[j];
                }
                templung = lungime;
            }
        }
    }

    //Rprintf("lungime: %d; templung: %d\n\n", lungime, templung);

    if (flag == 1) {

        templung = 0;

        for (i = 0; i < lungime; i++) {
            //Rprintf("%d ", ptemp2[i]);
            if (ptemp2[i] < (pmax[0] + 1)) {
                templung += 1;
            }
        }

        //Rprintf("lungime: %d; templung: %d; i: %d\n", lungime, templung, i);

        SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, templung - 1)); //"- 1" because the first element needs to be eliminated as well
        ptemp1 = INTEGER(temp1);

        j = 0;
        for (i = 1; i < lungime; i++) {
            if (ptemp2[i] < pmax[0] + 1) {
                ptemp1[j] = ptemp2[i];
                j += 1;
            }
        }
    }
    else {
        UNPROTECT(1);
        return(R_NilValue);
    }

    UNPROTECT(1);

    return(temp1);
}

SEXP C_pof(SEXP x, SEXP y, SEXP nec) {

    // x is a matrix of membership scores
    // y is the outcome, a vector of length equal to the nrows(x)
    SEXP pmin, max_ec;
    SEXP usage = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, y = coerceVector(y, REALSXP));

    double *p_x = REAL(x);
    double *p_y = REAL(y);

    int nrows_x = nrows(x);
    int ncols_x = ncols(x);

    // intersection (min) between each condition and the outcome
    SET_VECTOR_ELT(usage, 2, pmin = allocMatrix(REALSXP, nrows_x, ncols_x));
    double *p_pmin = REAL(pmin);

    // union (max) of all columns _e_xcept the _c_urrent one
    SET_VECTOR_ELT(usage, 3, max_ec = allocMatrix(REALSXP, nrows_x, ncols_x));
    double *p_max_ec = REAL(max_ec);

    // sum of the outcome
    double sum_y = 0;

    // sum of the conditions
    double sum_x[ncols_x];

    // sum of the negated conditions
    double sum_neg_x[ncols_x];

    // sum of the intersection between each condition and the outcome
    double sum_pmin[ncols_x];

    // sum of the minimum between each condition and the negation of the outcome
    double sum_pmin_negy[ncols_x];

    // sum of the negation of the minimum between each condition and the outcome
    double sum_neg_pmin[ncols_x];

    // for the Runique coverage, sum of the minimum between: union of all columns except the current one, and the outcome
    double sum_min_max_ec[ncols_x];

    for (int c = 0; c < ncols_x; c++) { // initialize values, especially to prevent valgrind errors
        sum_x[c] = 0.0;
        sum_neg_x[c] = 0.0;
        sum_pmin[c] = 0.0;
        sum_pmin_negy[c] = 0.0;
        sum_neg_pmin[c] = 0.0;
        sum_min_max_ec[c] = 0.0;
        for (int r = 0; r < nrows_x; r++) {
            p_max_ec[c * nrows_x + r] = 0.0;
        }
    }

    for (int r = 0; r < nrows_x; r++) {
        sum_y += p_y[r];
        for (int c = 0; c < ncols_x; c++) {
            sum_x[c] += p_x[c * nrows_x + r];
            sum_neg_x[c] += 1 - p_x[c * nrows_x + r];
            p_pmin[c * nrows_x + r] = (p_x[c * nrows_x + r] < p_y[r]) ? p_x[c * nrows_x + r] : p_y[r];
            sum_pmin[c] += p_pmin[c * nrows_x + r];
            sum_neg_pmin[c] += 1 - p_pmin[c * nrows_x + r];
            sum_pmin_negy[c] += (p_pmin[c * nrows_x + r] < (1 - p_y[r])) ? p_pmin[c * nrows_x + r] : (1 - p_y[r]);
        }
    }


    // pairwise max of all other columns except the current one (and also except the last column)
    for (int r = 0; r < nrows_x; r++) {
        for (int c = 0; c < ncols_x - 1; c++) { // ncols_x - 1 because the last column is the union of all columns (entire expression)
            for (int cu = 0; cu < ncols_x - 1; cu++) {
                if (cu != c) {
                    if (p_max_ec[c * nrows_x + r] < p_pmin[cu * nrows_x + r]) {
                        p_max_ec[c * nrows_x + r] = p_pmin[cu * nrows_x + r];
                    }
                }
            }
        }
    }

    // sum of the pairwise min between the current column, the outcome and the previously computed pairwise max
    // note: p_pmin already contains the pairwise min with the outcome
    for (int r = 0; r < nrows_x; r++) {
        for (int c = 0; c < ncols_x - 1; c++) {
            sum_min_max_ec[c] += (p_pmin[c * nrows_x + r] < p_max_ec[c * nrows_x + r]) ? p_pmin[c * nrows_x + r]: p_max_ec[c * nrows_x + r];
        }
    }


    SEXP inclcov;
    SET_VECTOR_ELT(usage, 4, inclcov = allocMatrix(REALSXP, ncols_x, 4));
    double *p_inclcov = REAL(inclcov);


    for (int c = 0; c < ncols_x; c++) {
        if (LOGICAL(nec)[0]) {
            // inclN
            p_inclcov[c] = sum_pmin[c] / sum_y;
            // RoN
            p_inclcov[ncols_x + c] = sum_neg_x[c] / sum_neg_pmin[c];
            // covN
            p_inclcov[2 * ncols_x + c] = sum_pmin[c] / sum_x[c];
            p_inclcov[3 * ncols_x + c] = 0; // irrelevant column
        }
        else {
            // inclS
            p_inclcov[c] = sum_pmin[c] / sum_x[c];
            // PRI
            p_inclcov[ncols_x + c] = (sum_pmin[c] - sum_pmin_negy[c]) / (sum_x[c] - sum_pmin_negy[c]);
            // covS
            p_inclcov[2 * ncols_x + c] = sum_pmin[c] / sum_y;
            // covU
            p_inclcov[3 * ncols_x + c] = p_inclcov[2 * ncols_x + c] - (sum_min_max_ec[c] / sum_y);
        }
    }


    UNPROTECT(1);
    return(inclcov);

}

typedef struct {
    int nconds;
    int *p_k;
    int *p_noflevels;
    int *p_result;
} QCAComplexityContext;

static int qca_complexity_one(
    int nconds,
    int k,
    const int *p_noflevels
) {
    int resum = 0;
    unsigned long long int maxtasks = nchoosek(nconds, k);

    for (unsigned long long int task = 0; task < maxtasks; task++) {
        int tempk[k];
        unsigned long long int combination = task;
        int x = 0;

        for (int i = 0; i < k; i++) {
            while (1) {
                unsigned long long int cval = nchoosek(nconds - (x + 1), k - (i + 1));
                if (cval == 0 || cval > combination) {
                    break;
                }
                combination -= cval;
                x++;
            }

            if (x < 0) {
                x = 0;
            }
            if (x >= nconds) {
                x = nconds - 1;
            }

            tempk[i] = x;
            x++;
        }

        int prod = 1;
        for (int i = 0; i < k; i++) {
            prod *= p_noflevels[tempk[i]];
        }

        resum += prod;
    }

    return resum;
}

static void qca_complexity_range_worker(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    QCAComplexityContext *ctx = (QCAComplexityContext *) data;
    (void) worker_id;

    for (unsigned long long ck = start; ck < end; ck++) {
        ctx->p_result[ck] = qca_complexity_one(
            ctx->nconds,
            ctx->p_k[ck],
            ctx->p_noflevels
        );
    }
}

SEXP C_omplexity(SEXP list) {
    /*
    nconds is k in the R code and the k is c.
    I've deliberately opted for this choice because the code below
    is a modified version of C_combnk which uses choose(n, k).
    */
    int nconds = INTEGER(VECTOR_ELT(list, 0))[0];
    int lk = length(VECTOR_ELT(list, 1));
    int *p_k = INTEGER(VECTOR_ELT(list, 1));
    int *p_noflevels = INTEGER(VECTOR_ELT(list, 2));

    int *p_result = (int *) calloc((size_t) lk, sizeof(int));

    if (p_result == NULL) {
        error("Memory allocation failed during complexity calculation.");
    }

    QCAComplexityContext ctx = {
        .nconds = nconds,
        .p_k = p_k,
        .p_noflevels = p_noflevels,
        .p_result = p_result
    };
    if (!qca_parallel_for((unsigned long long) lk, 0, qca_complexity_range_worker, &ctx)) {
        free(p_result);
        error("Failed to start pthread workers for complexity calculation.");
    }

    R_CheckUserInterrupt();

    SEXP result = PROTECT(allocVector(INTSXP, lk));
    for (int ck = 0; ck < lk; ck++) {
        INTEGER(result)[ck] = p_result[ck];
    }
    free(p_result);

    UNPROTECT(1);
    return(result);

}

SEXP C_expand(SEXP mat, SEXP noflevels, SEXP partial) {

    int nconds = ncols(mat); // number of conditions
    int nimp = nrows(mat); // number of implicants
    int *p_imp = INTEGER(mat);
    int *p_noflevels = INTEGER(noflevels);

    // Rprintf("nconds: %d; nrows: %d\n", nconds, nimp);
    // Rf_PrintValue(VECTOR_ELT(list, 0));

    // for (int i = 0; i < nconds; i++) {
    //     Rprintf("%d ", p_noflevels[i]);
    // }
    // Rprintf("\n");

    // check if any column has all values equal to zero
    // which means there is no point to expand there

    Rboolean zeroc[nconds];
    int nzero = 0;

    for (int c = 0; c < nconds; c++) {
        zeroc[c] = TRUE;
        int r = 0;
        while(zeroc[c] && r < nimp) {
            zeroc[c] = p_imp[c * nimp + r] == 0;
            r++;
        }

        if (zeroc[c]) {
            nzero++;
        }
    }

    int ncolsx = nconds - nzero;

    SEXP usage = PROTECT(allocVector(VECSXP, 7));
    SEXP x, rest, resmat, result, temp; //
    int rfilled = 0; // rows filled in the result matrix

    int estnrows = 1000; // estimated number of rows... TRANSPOSED!
    SET_VECTOR_ELT(usage, 0, result = allocMatrix(INTSXP, nconds, estnrows));
    int *p_result = INTEGER(result);
    Memzero(p_result, nconds * estnrows);

    SET_VECTOR_ELT(usage, 1, x = allocMatrix(INTSXP, nimp, ncolsx));
    int *p_x = INTEGER(x);
    int pos = 0;
    int poszero = 0;

    // this is the equivalent of the R code
    // x <- x[, -zeroc, drop = FALSE]
    int cnotzero[ncolsx];
    for (int c = 0; c < nconds; c++) {
        if (!zeroc[c]) {
            cnotzero[poszero] = c;
            poszero += 1;
            for (int r = 0; r < nimp; r++) {
                p_x[pos] = p_imp[c * nimp + r];
                pos++;
            }
        }
    }

    /*
    rmin <- min(apply(x, 1, function(x) sum(x == 0)))
    and to simulate:
    which(xi == 0)
    (that is "wxi"), I need a "m"atrix storing those positions
    */
    int rmin = ncolsx;
    int wxim[nimp * ncolsx]; // which(xi == 0)
    int rxi[nimp]; // rxi <- sum(xi > 0)

    for (int r = 0; r < nimp; r++) {
        rxi[r] = 0;
        for (int c = 0; c < ncolsx; c++) {
            wxim[c * nimp + r] = -1; // just to initialize
            if (p_x[c * nimp + r] == 0) {
                wxim[rxi[r] * nimp + r] = c;
                rxi[r] += 1;
            }
        }

        if (rxi[r] < rmin) {
            rmin = rxi[r];
        }
    }

    /*
    Rprintf("rmin: %d\n", rmin);

    for (int r = 0; r < nimp; r++) {
        for (int c = 0; c < ncolsx; c++) {
            Rprintf("%d ", wxim[c * nimp + r]);
        }
        Rprintf("\n");
    }
    Rprintf("\n");
    */

    for (int ri = 0; ri < nimp; ri++) {

        /*
        Rprintf("\n **** r%d: ", ri + 1);
        for (int c = 0; c < ncolsx; c++) {
            Rprintf("%d ", p_x[c * nimp + ri]);
        }
        Rprintf("\n\n");
        */

        int k = rxi[ri];
        // rmin is used here for partial expansion
        if (LOGICAL(partial)[0]) {
            k = k - rmin;
        }


        // Rprintf("rxi: %d; k: %d\n", rxi[ri], k);

        if (k > 0) {

            // now simulate
            // combs <- combnk(length(wxi), rmin - rxi)
            unsigned long long int maxtasks = nchoosek(rxi[ri], k);
            for (unsigned long long int task = 0; task < maxtasks; task++) {
                if (task > 0 && task % 1024 == 0) {
                    R_CheckUserInterrupt();
                }

                int tempk[k];
                unsigned long long int combination = task;
                int x = 0;

                for (int i = 0; i < k; i++) {
                    while (1) {
                        unsigned long long int cval = nchoosek(rxi[ri] - (x + 1), k - (i + 1));
                        if (cval == 0 || cval > combination) {
                            break;
                        }
                        combination -= cval;
                        x++;
                    }

                    if (x < 0) {
                        x = 0;
                    }
                    if (x >= rxi[ri]) {
                        x = rxi[ri] - 1;
                    }

                    tempk[i] = x;
                    x++;
                }

                // now we need to simulate
                // rest <- getMatrix(noflevels[wxic]) + 1

                //Rprintf("tempk: ");
                int nofl[k];
                for (int i = 0; i < k; i++) {
                    /*
                    Rprintf("%d; %d, %d, %d, %d;; ", tempk[i],
                                                    tempk[i] * nimp,
                                                    tempk[i] * nimp + ri,
                                                    wxim[tempk[i] * nimp + ri],
                                                    p_noflevels[wxim[tempk[i] * nimp + ri]]);
                    */
                   nofl[i] = p_noflevels[wxim[tempk[i] * nimp + ri]];
                }

                /*
                Rprintf("nofl: ");
                for (int i = 0; i < k; i++) {
                    Rprintf("%d ", nofl[i]);
                }
                Rprintf("\n");
                */

                int nrows;
                calculate_rows(&nrows, k, nofl, 0, k);

                // Rprintf("nrows: %d; ncols: %d\n", nrows, k);

                SET_VECTOR_ELT(usage, 2, rest = allocMatrix(INTSXP, nrows, k));
                int *p_rest = INTEGER(rest);

                generate_matrix(nrows, k, nofl, 0, k, p_rest);

                SET_VECTOR_ELT(usage, 3, resmat = allocMatrix(INTSXP, nrows, ncolsx));
                int *p_resmat = INTEGER(resmat);

                Rboolean cfilled[ncolsx];
                for (int c = 0; c < ncolsx; c++) {
                    cfilled[c] = FALSE;
                }

                for (int i = 0; i < k; i++) {
                    int c = wxim[tempk[i] * nimp + ri];
                    cfilled[c] = TRUE;
                    for (int r = 0; r < nrows; r++) {
                        p_resmat[c * nrows + r] = p_rest[i * nrows + r] + 1;
                    }
                }

                for (int c = 0; c < ncolsx; c++) {
                    if (!cfilled[c]) {
                        for (int r = 0; r < nrows; r++) {
                            p_resmat[c * nrows + r] = p_x[c * nimp + ri];
                        }
                    }
                }

                // Rprintf("rfilled: %d; nrows: %d; estnrows: %d\n", rfilled, nrows, estnrows);

                if ((rfilled + nrows) > estnrows) {
                    // Rprintf("%d, %d\n", ncolsx, estnrows);
                    estnrows *= 2;
                    // double the result matrix
                    // Rf_PrintValue(result);
                    // Rprintf("%d, %d\n", ncolsx, estnrows);
                    SET_VECTOR_ELT(usage, 3, result = Rresize(result, nconds * estnrows));
                    p_result = INTEGER(result);
                    // Rf_PrintValue(result);
                }

                // Rprintf("rfilled: %d; nrows: %d; estnrows: %d\n", rfilled, nrows, estnrows);

                // for (int c = 0; c < ncolsx; c++) {
                //     Rprintf("%d ", cnotzero[c]);
                // }
                // Rprintf("\n");

                for (int c = 0; c < ncolsx; c++) {
                    for (int r = 0; r < nrows; r++) {
                        // transposed...
                        // p_result[cnotzero[c] * estnrows + rfilled + r] = p_resmat[c * nrows + r];
                        p_result[(rfilled + r) * nconds + cnotzero[c]] = p_resmat[c * nrows + r];
                    }
                }

                // Rprintf("resmat:\n");
                // Rf_PrintValue(Rtranspose(resmat));
                rfilled += nrows;

                // Rprintf("result:\n");
                // Rf_PrintValue(result);

            }
        }
        else {

            if ((rfilled + 1) > estnrows) {
                estnrows *= 2;

                SET_VECTOR_ELT(usage, 3, result = Rresize(result, nconds * estnrows));
                p_result = INTEGER(result);
            }

            for (int c = 0; c < ncolsx; c++) {
                p_result[rfilled * nconds + cnotzero[c]] = p_x[c * nimp + ri];
            }

            rfilled += 1;
            // Rprintf("result:\n");
            // Rf_PrintValue(result);
        }
    }

    SET_VECTOR_ELT(usage, 4, temp = allocMatrix(INTSXP, nconds, rfilled));
    int *p_temp = INTEGER(temp);
    Memcpy(p_temp, p_result, nconds * rfilled);
    // for (int c = 0; c < nconds; c++) {
    //     for (int r = 0; r < rfilled; r++) {
    //         // transposing, because Runique() expects it transposed
    //         p_temp[r * nconds + c] = p_result[c * estnrows + r];
    //     }
    // }

    // Rf_PrintValue(result);

    SEXP unq, trsp;

    SET_VECTOR_ELT(usage, 5, unq = Runique(temp));
    SET_VECTOR_ELT(usage, 6, trsp = Rtranspose(unq));

    UNPROTECT(1);
    return(trsp);
}

SEXP C_simplify(SEXP mat, SEXP noflevels, SEXP partial) {

    SEXP umat = PROTECT(C_expand(mat, noflevels, partial));
    SEXP simplified = PROTECT(C_QMC(umat, noflevels));
    UNPROTECT(2);
    return(simplified);

}

SEXP C_Cubes(SEXP list) {

    // the tt outcome column should always be the last, and it should only contain 0s and 1s
    // the truth table is first, which is position // 0 in the zero-based C
    int posdata =      getpos(list, "data");       // 1
    int posallsol =    getpos(list, "all.sol");    // 2
    int posrowdom =    getpos(list, "row.dom");    // 3
    int pospicons =    getpos(list, "pi.cons");    // 4
    int posdepth =     getpos(list, "depth");      // 5
    int posolcons =    getpos(list, "sol.cons");   // 6
    int posolcov  =    getpos(list, "sol.cov");    // 7
    int posfs =        getpos(list, "fs");         // 8
    int posmaxcomb =   getpos(list, "max.comb");   // 9
    int pos1stmin =    getpos(list, "first.min");  // 10
    int poskeeptry =   getpos(list, "keep.trying");// 11
    int posgurobi =    getpos(list, "gurobi");     // 12
    int posolind =     getpos(list, "solind");     // 13
    int poslagrangian = getpos(list, "lagrangian");// 14


    SEXP usage = PROTECT(allocVector(VECSXP, 7));
    /*
        0:  tt
        1:  data
        2:  fsconds
        3:  pichart
        4:  dimnames
        5:  ttcolnms
        6:  colnms
    */

    SEXP   tt, data,    fsconds;

    SET_VECTOR_ELT(usage, 0, tt = coerceVector(VECTOR_ELT(list, 0), INTSXP));
    int *p_tt = INTEGER(tt);

    if (posdata > 0) {
        SET_VECTOR_ELT(usage, 1, data = coerceVector(VECTOR_ELT(list, posdata), REALSXP));
    }
    else {
        SET_VECTOR_ELT(usage, 1, data = allocMatrix(REALSXP, 2, 2)); // just initialise
        Memzero(REAL(data), 4);
    }

    double *p_data = REAL(data);

    int nrdata = nrows(data);
    // int ncdata = ncols(data); // not needed, this would be nconds + 1

    int ttrows = nrows(tt); // number of rows in the list data
    int nconds = ncols(tt) - 1; // number of conditions (columns - 1 because the last one is the outcome)


    Rboolean allsol = (posallsol >= 0) ? (LOGICAL(VECTOR_ELT(list, posallsol))[0]) : FALSE;
    Rboolean rowdom = (posrowdom >= 0) ? (LOGICAL(VECTOR_ELT(list, posrowdom))[0]) : FALSE;
    Rboolean keeptrying = (poskeeptry >= 0) ? (LOGICAL(VECTOR_ELT(list, poskeeptry))[0]) : FALSE;
    Rboolean firstmin = (pos1stmin >= 0) ? (LOGICAL(VECTOR_ELT(list, pos1stmin))[0]) : FALSE;
    Rboolean gurobi = (posgurobi >= 0) ? (LOGICAL(VECTOR_ELT(list, posgurobi))[0]) : TRUE;
    Rboolean solind = (posolind >= 0) ? (LOGICAL(VECTOR_ELT(list, posolind))[0]) : TRUE;
    Rboolean lagrangian = (poslagrangian >= 0) ? (LOGICAL(VECTOR_ELT(list, poslagrangian))[0]) : FALSE;
    double picons = (pospicons >= 0) ? (REAL(VECTOR_ELT(list, pospicons))[0]) : 0;

    int pidepth = 0;
    int soldepth = 5; // a default arbitrary value, but this is set by minimize() anyways

    if (posdepth >= 0) {
        pidepth = INTEGER(coerceVector(VECTOR_ELT(list, posdepth), INTSXP))[0];
        soldepth = INTEGER(coerceVector(VECTOR_ELT(list, posdepth), INTSXP))[1];
    }

    if (pidepth == 0 || nconds < pidepth) {
        pidepth = nconds;
    }

    double solcons = (posolcons >= 0) ? (REAL(VECTOR_ELT(list, posolcons))[0]) : 0;
    double solcov = (posolcov >= 0) ? (REAL(VECTOR_ELT(list, posolcov))[0]) : 0;

    int *p_fsconds;
    if (posfs > 0) {
        SET_VECTOR_ELT(usage, 2, fsconds = coerceVector(VECTOR_ELT(list, posfs), INTSXP));
        p_fsconds = INTEGER(fsconds);
    }
    else {
        SET_VECTOR_ELT(usage, 2, fsconds = allocVector(INTSXP, nconds)); // just initialise
        p_fsconds = INTEGER(fsconds);
        Memzero(p_fsconds, nconds);
    }



    int *p_pichart = R_Calloc(1, int); // just to initialize
    int *p_impmat = R_Calloc(1, int);  // the PIs, into an implicant matrix
    int *p_models = R_Calloc(1, int);  // p_solutions in fact, just distinguishing between "solution" (types) and "model" (the actual disjunction of PIs)
    unsigned int foundPI = 0;
    int solrows = 0;
    int solcols = 0;
    Rboolean complexpic = FALSE;

    double maxcomb = 0;
    if (posmaxcomb > 0) {
        maxcomb = REAL(VECTOR_ELT(list, posmaxcomb))[0]; // in millions
    }



    CCubes(
        p_tt, ttrows, nconds, p_data, nrdata, allsol, rowdom, picons, pidepth, p_fsconds, soldepth, solcons, solcov, maxcomb, keeptrying,

        &p_pichart, &p_impmat, &p_models, &foundPI, &solrows, &solcols, &complexpic,

        firstmin,
        lagrangian,
        gurobi,
        solind
    );


    // calculate the number of positive rows
    int posrows = 0;
    for (int r = 0; r < ttrows; r++) {
        posrows += p_tt[nconds * ttrows + r];
    }

    // for (int r = 0; r < foundPI; r++) {
    //     for (int c = 0; c < posrows; c++) {
    //         Rprintf("%d ", p_pichart[c * foundPI + r]);
    //     }
    //     Rprintf("  ");
    //     for (int c = 0; c < nconds; c++) {
    //         Rprintf("%d ", p_impmat[c * foundPI + r]);
    //     }
    //     Rprintf("\n");
    // }


    SEXP implicants, pic, solmat, complex;
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    int i, j, l_1, len;

    if (foundPI > 0) {
        SET_VECTOR_ELT(out, 0, implicants = allocMatrix(INTSXP, foundPI, nconds));
        int *p_implicants = INTEGER(implicants);
        Memzero(p_implicants, foundPI * nconds); // set everything to zero

        len = foundPI * nconds;
        l_1 = len - 1;
        for (i = 0, j = 0; i < len; i++, j += nconds) {
            if (j > l_1) j -= l_1;
            p_implicants[i] = p_impmat[j];
        }

        // transposed to a row-major with PIs on the rows
        SET_VECTOR_ELT(out, 1, pic = allocMatrix(LGLSXP, foundPI, posrows));
        int *p_pic = LOGICAL(pic);

        len = foundPI * posrows;
        l_1 = len - 1;
        for (i = 0, j = 0; i < len; i++, j += posrows) {
            if (j > l_1) j -= l_1;
            p_pic[i] = p_pichart[j];
        }

        if (hasColnames(tt)) {
            SEXP dimnames, ttcolnms,  colnms;
            SET_VECTOR_ELT(usage, 4, dimnames = allocVector(VECSXP, 2));
            SET_VECTOR_ELT(usage, 5, ttcolnms = VECTOR_ELT(getAttrib(tt, R_DimNamesSymbol), 1));
            SET_VECTOR_ELT(usage, 6, colnms = allocVector(STRSXP, nconds));

            for (int i = 0; i < nconds; i++) {
                SET_STRING_ELT(colnms, i, STRING_ELT(ttcolnms, i));
            }

            SET_VECTOR_ELT(dimnames, 1, colnms); // 1 for the columns
            setAttrib(implicants, R_DimNamesSymbol, dimnames);
        }
    }

    // Rprintf("solrows: %d; solcols: %d\n", solrows, solcols);

    if (solrows > 0 && solcols > 0) { // else it stays NULL
        SET_VECTOR_ELT(out, 2, solmat = allocMatrix(INTSXP, solrows, solcols));

        int *p_solmat = INTEGER(solmat);
        if (firstmin) {
            for (int i = 0; i < solrows; i++) {
                // p_solmat[i] = i + 1;
                // Rprintf("%d ", p_models[i]);
                p_solmat[i] = p_models[i] + 1;
            }
            // Rprintf("\n");
        }
        else {
            Memcpy(p_solmat, p_models, solrows * solcols);
        }
    }

    SET_VECTOR_ELT(out, 3, complex = allocVector(LGLSXP, 1));
    LOGICAL(complex)[0] = complexpic;

    // Rprintf("over\n");

    R_Free(p_pichart);
    R_Free(p_impmat);
    R_Free(p_models);

    UNPROTECT(2);
    return(out);
}

// no other function after this line
// C_getEC has to always be the LAST function for bwc()
SEXP C_getEC(SEXP dem, SEXP cexpr, SEXP csolm, SEXP pexpr, SEXP psolm, SEXP SA, SEXP noflevels) {

    /*
    dem = directional expectations matrix
    cexpr = conservative solution expressions
    csolm = conservative solution matrix
    pexpr = parsimonious solution expressions
    psolm = parsimonious solution matrix
    SA = simplifying assumptions
    */

    SEXP usage = PROTECT(allocVector(VECSXP, 8));

    int *p_dem   = INTEGER(dem);
    int *p_cexpr = INTEGER(cexpr);
    int *p_csolm = INTEGER(csolm);
    int *p_pexpr = INTEGER(pexpr);
    int *p_psolm = INTEGER(psolm);

    int ncols_dem   = ncols(dem);
    int nrows_dem   = nrows(dem);
    int nrows_cexpr = nrows(cexpr);
    int ncols_csolm = ncols(csolm);
    int nrows_csolm = nrows(csolm);
    int nrows_pexpr = nrows(pexpr);
    int ncols_psolm = ncols(psolm);
    int nrows_psolm = nrows(psolm);

    int rowsums[nrows_dem];
    int idecount = 0;
    int cdecount = 0;


    /*
    This code makes a count on how many DEs are on each row of the DE matrix.
    If a single DE, it is an individual IDE, and if multiple DEs on a row, it is a conjunctural CDE
    */
    for (int r = 0; r < nrows_dem; r++) {
        rowsums[r] = 0;
        for (int c = 0; c < ncols_dem; c++) {
            // Rprintf("%d ", 1*(p_dem[c * nrows_dem + r] > 0));
            if (p_dem[c * nrows_dem + r] > 0) {
                rowsums[r] += 1;
            }
        }

        if (rowsums[r] == 1) {
            idecount += 1;
        } else if (rowsums[r] > 1) {
            cdecount += 1;
        }
        // Rprintf("\n");
    }


    /*
    Using the information above, create and store the row indices for the IDEs and CDEs
    */
    SEXP iderows, cderows;
    SET_VECTOR_ELT(usage, 1, iderows = allocVector(INTSXP, idecount));
    int *p_iderows = INTEGER(iderows);
    SET_VECTOR_ELT(usage, 2, cderows = allocVector(INTSXP, cdecount));
    int *p_cderows = INTEGER(cderows);

    for (int i = 0; i < idecount; i++) {
        p_iderows[i] = 0; // initializing just to avoid a segmentation fault
    }

    for (int i = 0; i < cdecount; i++) {
        p_cderows[i] = 0; // initializing just to avoid a segmentation fault
    }

    int ipos = 0, cpos = 0;
    for (int r = 0; r < nrows_dem; r++) {
        if (rowsums[r] == 1) {
            p_iderows[ipos] = r;
            ipos++;
        } else if (rowsums[r] > 1) {
            p_cderows[cpos] = r;
            cpos++;
        }
    }

    /*
    When searching directional expectations to form the dir.exp.matrix as in R,
    this information should be enough. In R, I had a list, but here I will check
    which values are greater than zero on all idecount rows
    (the DE matrix is specified in implicants form)

    Same for the cdecount rows, where the CDE row indices are stored in p_cderows,
    but here we do _not_ need to check for multiple values, as simplify() and
    translate() already put all combinations on different rows
    Something like: simplify("A{1} + B{0}C{1,2}") is "A{1} + B{0}*C{1} + B{0}*C{2}"
    */
    
    SEXP tempmat;
    int *p_tempmat;
    int *p_intseltemp;

    SEXP output = PROTECT(allocVector(VECSXP, 3));
    SEXP EClist = PROTECT(allocVector(VECSXP, ncols_csolm * ncols_psolm));
    SEXP primes = PROTECT(allocVector(VECSXP, ncols_csolm * ncols_psolm));
    SEXP intsel = PROTECT(allocVector(VECSXP, ncols_csolm * ncols_psolm));

    int comp[ncols_dem];
    int pars[ncols_dem];
    int res[ncols_dem];

    SEXP sap; // Simplifying Assumptions from the Parsimonious solution
    int *p_sap;

    int tempcount = 0;
    int ECpos = 0;

    for (int csol = 0; csol < ncols_csolm; csol++) { //for (c.s in seq(ncol(c.sol$sol.matrix))) {
        /*
        Since the solution matrix can have different solution lengths (due to all.sol = TRUE)
        we have to check how many PIs are there in the current solution (current column from the csolm)
        */
        int csolcount = 0;
        for (int r = 0; r < nrows_csolm; r++) {
            if (p_csolm[csol * nrows_csolm + r] > 0) {
                csolcount++;
            }
        }


        for (int psol = 0; psol < ncols_psolm; psol++) { //for (p.s in seq(ncol(p.sol$sol.matrix))) {

            // same as above
            int psolcount = 0;
            for (int r = 0; r < nrows_psolm; r++) {
                if (p_psolm[psol * nrows_psolm + r] > 0) {
                    psolcount++;
                }
            }

            SET_VECTOR_ELT(usage, 6, tempmat = allocMatrix(INTSXP, ncols_dem, csolcount * psolcount));
            p_tempmat = INTEGER(tempmat);
            tempcount = 0;

            SET_VECTOR_ELT(usage, 3, sap = VECTOR_ELT(SA, psol));
            p_sap = INTEGER(sap);
            int nrows_sap = nrows(sap);

            Rboolean ecs[nrows_sap];
            for (int r = 0; r < nrows_sap; r++) {
                ecs[r] = FALSE;
            }

            Rboolean subset[csolcount];

            for (int ics = 0; ics < csolcount; ics++) {
                subset[ics] = FALSE; // assume it is not a subset

                // Rprintf("comp:\n");
                for (int c = 0; c < ncols_dem; c++) {
                    comp[c] = p_cexpr[c * nrows_cexpr + p_csolm[csol * nrows_csolm + ics] - 1];
                    // Rprintf("%d ", comp[c]);
                }
                // Rprintf("\n");

                // Rprintf("pars:\n");
                for (int ips = 0; ips < psolcount; ips++) {
                    Rboolean allequal = TRUE;
                    Rboolean parsgz[ncols_dem]; // PS terms greater than zero
                    Rboolean compgz[ncols_dem]; // CS terms greater than zero, plays the role of <notequals> in R

                    int c = 0;
                    while (allequal && c < ncols_dem) {
                        pars[c] = p_pexpr[c * nrows_pexpr + p_psolm[psol * nrows_psolm + ips] - 1];
                        // Rprintf("%d ", pars[c]);
                        parsgz[c] = pars[c] > 0;
                        compgz[c] = comp[c] > 0 && !parsgz[c]; // when pars[c] is zero, comp[c] can be a <notequal in R>
                        if (parsgz[c]) {
                            allequal = pars[c] == comp[c];
                        }
                        c++;
                    }

                    // Rprintf("\nallequal: %d\n", allequal * 1);

                    if (allequal) { // if (all(comp[pars > 0] == pars[pars > 0]))

                        subset[ics] = TRUE;

                        if (idecount > 0) { // Individual Directional Expectations

                            for (int c = 0; c < ncols_dem; c++) {
                                res[c] = comp[c]; // assume no EC will be found
                                if (compgz[c]) {
                                    // only for those conditions which are <not> part of the subset relation
                                    Rboolean dontcare = TRUE;
                                    Rboolean defound = FALSE;

                                    for (int r = 0; r < idecount; r++) {
                                        int inde = p_dem[c * nrows_dem + p_iderows[r]];
                                        if (inde > 0) {
                                            dontcare = FALSE;
                                            if (inde == comp[c]) {
                                                defound = TRUE;
                                            }
                                        }
                                    }

                                    if (!dontcare && !defound) {
                                        // there <is> a DE for that condition, and it is <not> equal to the value from the complex term
                                        res[c] = 0;
                                    }
                                }
                            }
                        }

                        if (cdecount > 0) { // Conjunctural Directional Expectations

                            for (int r = 0; r < cdecount; r++) {
                                Rboolean allcde = TRUE;

                                for (int c = 0; c < ncols_dem; c++) {
                                    res[c] = comp[c];
                                    if (compgz[c]) { // which means pars[c] is equal to 0
                                        int inde = p_dem[c * nrows_dem + p_cderows[r]];
                                        if (inde > 0) {
                                            allcde = allcde && inde == comp[c];
                                        }
                                    }
                                }

                                if (!allcde) {
                                    for (int c = 0; c < ncols_dem; c++) {
                                        if (compgz[c] && p_dem[c * nrows_dem + p_cderows[r]] > 0) {
                                            res[c] = 0;
                                        }
                                    }
                                }
                            }
                        }


                        for (int r = 0; r < ncols_dem; r++) { // populate the intermediate selection terms (res)
                            p_tempmat[tempcount * ncols_dem + r] = res[r];
                        }

                        tempcount++;

                        for (int r = 0; r < nrows_sap; r++) {
                            if (!ecs[r]) {
                                Rboolean equal = TRUE;
                                int c = 0;
                                while (equal && c < ncols_dem) {
                                    equal = (res[c] > 0) ? p_sap[c * nrows_sap + r] + 1 == res[c] : TRUE;
                                    c++;
                                }
                                ecs[r] = equal;
                            }
                        }
                    }
                } // end comparing the CS terms with all PS terms
            } // end comparing the CS terms


            // Rprintf("EC: \n");
            // for (int r = 0; r < nrows_sap; r++) {
            //     if (ecs[r]) {
            //         for (int c = 0; c < ncols_dem; c++) {
            //             Rprintf("%d ", p_sap[c * nrows_sap + r]);
            //         }
            //         Rprintf("\n");
            //     }
            // }
            // Rprintf("\n");


            int ECrows = 0;
            for (int r = 0; r < nrows_sap; r++) {
                ECrows += ecs[r] * 1;
            }


            // Rprintf("ECrows: %d\n", ECrows);
            int countnot = 0; // count the CS terms that are not subsets of any PS terms
            for (int c = 0; c < csolcount; c++) {
                if (!subset[c]) {
                    countnot++;
                }
            }

            SEXP intseltemp;
            int totalrows = tempcount + countnot;
            // Rprintf("totalrows (%d) =  tempcount(%d) + countnot (%d)\n", totalrows, tempcount, countnot);
            SET_VECTOR_ELT(usage, 7, intseltemp = allocMatrix(INTSXP, totalrows, ncols_dem));
            // SET_VECTOR_ELT(intsel, ECpos, intseltemp = allocMatrix(INTSXP, totalrows, ncols_dem));
            p_intseltemp = INTEGER(intseltemp);
            Memzero(p_intseltemp, totalrows * ncols_dem);

            for (int c = 0; c < tempcount; c++) {
                for (int r = 0; r < ncols_dem; r++) {
                    // transpose the tempmat matrix
                    p_intseltemp[r * totalrows + c] = p_tempmat[c * ncols_dem + r];
                }
            }

            // for (int c = 0; c < csolcount; c++) {
            //     Rprintf("%d ", subset[c]*1);
            // }
            // Rprintf("\n");
            // Rf_PrintValue(csolm);
            // Rf_PrintValue(cexpr);

            // comp[c] = p_cexpr[c * nrows_cexpr + p_csolm[csol * nrows_csolm + ics] - 1];

            // Rf_PrintValue(intseltemp);

            if (countnot > 0) { // add those expressions from the conservative solutions that no PS term is a superset for
                for (int r = 0; r < csolcount; r++) {
                    if (!subset[r]) {
                        for (int c = 0; c < ncols_dem; c++) {
                            p_intseltemp[c * totalrows + tempcount] = p_cexpr[c * nrows_cexpr + p_csolm[csol * nrows_csolm + r] - 1];
                        }
                        tempcount++;
                    }
                }
            }

            SEXP dim_names;

            // if (csol == ncols_csolm - 1 && psol == ncols_psolm - 1) {
            //     Rf_PrintValue(intseltemp);
            // }


            SEXP pi;
            SEXP intselmat;
            SEXP partial;
            SET_VECTOR_ELT(usage, 4, partial = allocVector(LGLSXP, 1));
            LOGICAL(partial)[0] = FALSE;
            SET_VECTOR_ELT(intsel, ECpos, intselmat = intseltemp);
            SET_VECTOR_ELT(primes, ECpos, pi = C_simplify(intseltemp, noflevels, partial));

            /*
            ---------------------------------------------------------------------------------
            // this was supposed to eliminate any PI that is not a subset of the
            // intermediate solution terms... but it doesn't work because the
            // later PI chart (against the minterms) cannot always be solved

            int *p_pi = INTEGER(pi);

            int pin  = nrows(pi);
            Rboolean subsetpi[pin];
            for (int i = 0; i < pin; i++){
                subsetpi[i] = FALSE;
            }
            int sumset = 0;

            int i = 0;
            while (i < pin && !subsetpi[i]) {

                for (int j = 0; j < totalrows; j++) {
                    Rboolean allequal = TRUE;

                    int c = 0;
                    while (allequal && c < ncols_dem) {

                        // here, comp will play the role of the pi
                        // and pars will play the role of the intermediate selection terms

                        pars[c] = p_intseltemp[c * totalrows + j];

                        if (pars[c] > 0) {
                            allequal = pars[c] == p_pi[c * pin + i];
                        }
                        c++;
                    }

                    if (allequal) {
                        subsetpi[i] = TRUE;
                        sumset++;
                    }
                }

                i++;
            }

            if (sumset < pin) {
                SEXP tempcopy;
                SET_VECTOR_ELT(usage, 4, tempcopy = allocMatrix(INTSXP, sumset, ncols_dem));
                int *p_tempcopy = INTEGER(tempcopy);
                int pos = 0;
                for (int r = 0; r < pin; r++) {
                    if (subsetpi[r]) {
                        for (int c = 0; c < ncols_dem; c++) {
                            p_tempcopy[c * sumset + pos] = p_pi[c * pin + r];
                        }
                        pos++;
                    }
                }
                SET_VECTOR_ELT(primes, ECpos, pi = tempcopy);
            }
            ---------------------------------------------------------------------------------
            */


            // SET_VECTOR_ELT(intsel, ECpos, pi = intseltemp);
            SET_VECTOR_ELT(usage, 0, dim_names = allocVector(VECSXP, 2));
            if (hasColnames(cexpr)) {
                SET_VECTOR_ELT(dim_names, 1, VECTOR_ELT(getAttrib(cexpr, R_DimNamesSymbol), 1));
            }
            setAttrib(pi, R_DimNamesSymbol, dim_names);
            setAttrib(intselmat, R_DimNamesSymbol, dim_names);

            SEXP ec;
            SET_VECTOR_ELT(EClist, ECpos, ec = allocMatrix(INTSXP, ECrows, ncols_dem));
            int *p_ec = INTEGER(ec);

            SEXP sapnms;
            SET_VECTOR_ELT(usage, 4, sapnms = VECTOR_ELT(getAttrib(sap, R_DimNamesSymbol), 0));


            SEXP rownms;
            SET_VECTOR_ELT(usage, 5, rownms = allocVector(STRSXP, ECrows));

            SET_VECTOR_ELT(usage, 0, dim_names = allocVector(VECSXP, 2));

            int ecr = 0;
            // Rprintf("r: ");
            for (int r = 0; r < nrows_sap; r++) {
                if (ecs[r]) {
                    // Rprintf("%d, ", r + 1);
                    for (int c = 0; c < ncols_dem; c++) {
                        p_ec[c * ECrows + ecr] = p_sap[c * nrows_sap + r];
                    }
                    SET_STRING_ELT(rownms, ecr, STRING_ELT(sapnms, r));
                    ecr++;
                }
            }


            // Rprintf("\n");

            if (ecr > 0) {
                SET_VECTOR_ELT(dim_names, 0, duplicate(rownms));
                if (hasColnames(cexpr)) {
                    SET_VECTOR_ELT(dim_names, 1, VECTOR_ELT(getAttrib(cexpr, R_DimNamesSymbol), 1));
                }
                setAttrib(ec, R_DimNamesSymbol, dim_names);
            }

            ECpos++;

        } // end of psol (p.s)
    } // end of csol (c.s)

    SET_VECTOR_ELT(output, 0, EClist);
    SET_VECTOR_ELT(output, 1, primes);
    SET_VECTOR_ELT(output, 2, intsel);

    UNPROTECT(5);
    return(output);
    // return(R_NilValue);
}

// no other function after this line
