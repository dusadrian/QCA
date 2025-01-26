/*
Copyright (c) 2016 - 2024, Adrian Dusa
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

#include <R_ext/RS.h>
#include <R_ext/Boolean.h>
#include <float.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <R_ext/Rdynload.h>
#include "utils.h"
#include "find_models.h"
#include "generate_matrix.h"
#include "sort_matrix.h"
#include "CCubes.h"
static R_INLINE SEXP Rtranspose(SEXP matrix) {
    SEXPTYPE type = TYPEOF(matrix);
    int nr = nrows(matrix);
    int nc = ncols(matrix);
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
    }
    else if (type == LGLSXP) {
        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            LOGICAL(out)[i] = LOGICAL(matrix)[j];
        }
    }
    else if (type == REALSXP) {
        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            REAL(out)[i] = REAL(matrix)[j];
        }
    }
    UNPROTECT(3);
    return(out);
}
static R_INLINE Rboolean hasColnames(SEXP matrix) {
    if (Rf_isNull(getAttrib(matrix, R_DimNamesSymbol))) {
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
    Rboolean objlogical = isLogical(obj); 
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
static R_INLINE SEXP Runique(SEXP mat) {
    int nc = ncols(mat);
    int nr = nrows(mat);
    Rboolean logmat = isLogical(mat);
    int *p_mat = logmat ? LOGICAL(mat) : INTEGER(mat);
    SEXP umat;
    SEXP usage = PROTECT(allocVector(VECSXP, 1));
    Rboolean survived[nc];
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
        Memcpy(INTEGER(models), p_solutions, nr * nc);
        SEXP toocomplex;
        SET_VECTOR_ELT(out, 1, toocomplex = allocVector(LGLSXP, 1));
        LOGICAL(toocomplex)[0] = too_complex(foundPI, INTEGER(k)[0], REAL(maxcomb)[0]);
    }
    R_Free(p_solutions);
    UNPROTECT(2);
    return(out);
}
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
    for (int r = 0; r < nrows; r++) {
        for (int c = 0; c < ncols; c++) {
            p_matrix[c * nrows + r] = (p_rowno[r] / p_mbase[c]) % p_noflevels[c];
        }
    }
    UNPROTECT(2);
    return(matrix);
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
    if (length(input) > 3) { 
        SEXP dimnames;
        SET_VECTOR_ELT(usage, 3, dimnames = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 1, VECTOR_ELT(input, 3)); 
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
    int estim2 = 1000; 
    int xrows = nrows(x);
    int xcols = ncols(x);
    int nconds = ncols(x); 
    if (p_depth[0] == 0) {
        p_depth[0] = nconds;
    }
    SEXP tmconj, tmdisj, ticpr_conj, ticpr_disj, combkl, tcoms_conj, tcoms_disj,
         indx_conj, indx_disj, ck_conj, ck_disj;
    SET_VECTOR_ELT(usage,  8, tmconj     = allocVector(INTSXP, nconds * estim1));
    SET_VECTOR_ELT(usage,  9, tmdisj     = allocVector(INTSXP, nconds * estim2));
    SET_VECTOR_ELT(usage, 10, ticpr_conj = allocVector(REALSXP, 3 * estim1));    
    SET_VECTOR_ELT(usage, 11, ticpr_disj = allocVector(REALSXP, 3 * estim2));    
    SET_VECTOR_ELT(usage, 12, tcoms_conj = allocMatrix(REALSXP, xrows, estim1)); 
    SET_VECTOR_ELT(usage, 13, tcoms_disj = allocMatrix(REALSXP, xrows, estim1)); 
    SET_VECTOR_ELT(usage, 14, indx_conj  = allocVector(INTSXP, p_depth[0] * estim1)); 
    SET_VECTOR_ELT(usage, 15, indx_disj  = allocVector(INTSXP, p_depth[0] * estim2)); 
    SET_VECTOR_ELT(usage, 16, ck_conj    = allocVector(INTSXP, estim1)); 
    SET_VECTOR_ELT(usage, 17, ck_disj    = allocVector(INTSXP, estim2)); 
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
    double copyline[nconds], minx[xrows], maxx[xrows];
    double incovpron[6];
    double so = 0.0,
           min, max,
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
    for (int i = 0; i < length(vo); i++) {
        so += p_vo[i];
    }
    int chkred[nconds], inclcov;
    int k = 1;
    int foundk = 1;
    while (k <= p_depth[0] && foundk) {
        if (found1 + found > 0 && k > 3) {
            foundk = 0;
        }
        int klnofl[k];
        int tempk[k];
        for (int i = 0; i < k; i++) {
            tempk[i] = i;
        }
        int e = 0;
        int h = k;
        int kmatcol = 0;
        while (tempk[0] != nconds - k || !kmatcol) {
            if (kmatcol) {
                increment(k, &e, &h, nconds, tempk, 0);
            }
            kmatcol = 1;
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
                    chkred[tempk[j]] = p_combkl[j * klrows + kli] + 1; 
                }
                sum_minx = 0;         
                sum_maxx = 0;         
                sum_min_y_minx = 0;   
                sum_min_y_maxx = 0;   
                prisum_minx = 0;      
                prisum_maxx = 0;      
                sum_1_minx = 0;       
                sum_1_min_y_minx = 0; 
                sum_1_maxx = 0;       
                sum_1_min_y_maxx = 0; 
                for (int r = 0; r < xrows; r++) { 
                    min = 1000000;        
                    max = 0;
                    for (int c = 0; c < xcols; c++) { 
                        copyline[c] = p_x[c * xrows + r];
                        if (p_fuz[c]) { 
                            if (chkred[c] == 1) {
                                copyline[c] = 1 - copyline[c];
                            }
                        }
                        else {
                            if (chkred[c] == (copyline[c] + 1)) {
                                copyline[c] = 1; 
                            }
                            else {
                                copyline[c] = 0; 
                            }
                        }
                        if (chkred[c] != 0) {
                            if (copyline[c] < min) {
                                min = copyline[c]; 
                            }
                            if (copyline[c] > max) {
                                max = copyline[c]; 
                            }
                        }
                    } 
                    minx[r] = min;   
                    maxx[r] = max;   
                    sum_minx += min; 
                    sum_maxx += max; 
                    sum_min_y_minx += (min < p_vo[r])?min:p_vo[r];
                    sum_min_y_maxx += (max < p_vo[r])?max:p_vo[r];
                    if (p_nec[0]) {  
                        sum_1_minx += 1 - min;                                  
                        sum_1_maxx += 1 - max;                                  
                        sum_1_min_y_minx += 1 - ((min < p_vo[r])?min:p_vo[r]);  
                        sum_1_min_y_maxx += 1 - ((max < p_vo[r])?max:p_vo[r]);  
                    }
                    else {           
                        tmpv11 = (min < p_vo[r])?min:p_vo[r];
                        tmpv12 = p_nec[0]?(1 - min):(1 - p_vo[r]);
                        prisum_minx += (tmpv11 < tmpv12)?tmpv11:tmpv12;
                        tmpv21 = (max < p_vo[r])?max:p_vo[r];
                        tmpv22 = 1 - max;
                        prisum_maxx += (tmpv21 < tmpv22)?tmpv21:tmpv22;
                    }
                } 
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
                inclcov = incovpron[p_nec[0]] >= p_inclcut[0] && incovpron[1 - p_nec[0]] >= p_covcut[0];
                int redundant = 0;
                if (inclcov) {
                    if (foundk1 > 0 && !p_nec[0]) { 
                        int i = 0;
                        while (i < foundk1 && !redundant) {
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
                            i += 1;
                        }
                    }
                    if (!redundant) { 
                        for (int c = 0; c < nconds; c++) {
                            p_tmconj[found1 * nconds + c] = chkred[c];
                        }
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
                        inclcov = incovpron[2] >= p_inclcut[0] && incovpron[3] >= p_covcut[0];
                        redundant = 0;
                        if (inclcov && foundk1 > 0) {
                            int i = 0;
                            while (i < foundk1 && !redundant) {
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
                                i += 1;
                            }
                        }
                        if (inclcov && foundk2 > 0 && !redundant) {
                            int i = 0;
                            while (i < foundk2 && !redundant) {
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
                                i += 1;
                            }
                        }
                        if (inclcov && !redundant) {
                            for (int c = 0; c < nconds; c++) {
                                p_tmdisj[found * nconds + c] = chkred[c];
                            }
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
                } 
            } 
        }
        foundk1 = found1;
        foundk2 = found;
        k++;
    }
    SEXP icpr_conj, icpr_disj, mconj, mdisj, coms_conj, coms_disj;
    SEXP result = PROTECT(allocVector(VECSXP, 6));
    SET_VECTOR_ELT(result, 0, icpr_conj = allocMatrix(REALSXP, found1, 3)); 
    SET_VECTOR_ELT(result, 1, icpr_disj = allocMatrix(REALSXP, found, 3)); 
    SET_VECTOR_ELT(result, 2, mconj = allocMatrix(INTSXP, found1, nconds)); 
    SET_VECTOR_ELT(result, 3, mdisj = allocMatrix(INTSXP, found, nconds)); 
    SET_VECTOR_ELT(result, 4, coms_conj = allocMatrix(REALSXP, xrows, found1)); 
    SET_VECTOR_ELT(result, 5, coms_disj = allocMatrix(REALSXP, xrows, found)); 
    double *p_icpr_conj = REAL(icpr_conj);
    double *p_icpr_disj = REAL(icpr_disj);
    int *p_mconj = INTEGER(mconj);
    int *p_mdisj = INTEGER(mdisj);
    double *p_coms_conj = REAL(coms_conj);
    double *p_coms_disj = REAL(coms_disj);
    for (int r = 0; r < found1; r++) { 
        for (int c = 0; c < 3; c++) {
            p_icpr_conj[c * found1 + r] = p_ticpr_conj[r * 3 + c];
        }
    }
    for (int r = 0; r < found1; r++) { 
        for (int c = 0; c < nconds; c++) {
            p_mconj[c * found1 + r] = p_tmconj[r * nconds + c];
        }
    }
    for (int r = 0; r < found; r++) { 
        for (int c = 0; c < 3; c++) {
            p_icpr_disj[c * found + r] = p_ticpr_disj[r * 3 + c];
        }
    }
    for (int r = 0; r < found; r++) { 
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
    SEXP pimat, tempmat, minimized, copymat, order, cl; 
    int *p_tt, *p_noflevels, *p_pimat, *p_tempmat, *p_minimized, *p_copymat, *p_order,  *p_cl;
    SEXP usage = PROTECT(allocVector(VECSXP, 10));
    SEXP dimnames, colnms;
    SET_VECTOR_ELT(usage, 8, dimnames = getAttrib(tt, R_DimNamesSymbol));
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
            p_tempmat[r * nconds + c] = p_tt[c * nimplicants + r];  
        }
    }
    int found = 1;
    int temp[nconds];
    int combs[2];
    int e, h;
    while (found > 0 && nimplicants > 1) {
        found = 0;
        SET_VECTOR_ELT(usage, 3, minimized = allocVector(LGLSXP, nimplicants));
        p_minimized = LOGICAL(minimized);
        memset(p_minimized, FALSE, nimplicants * sizeof(int));
        int estimpi = 10000;
        SET_VECTOR_ELT(usage, 4, pimat = allocMatrix(INTSXP, nconds, estimpi));
        p_pimat = INTEGER(pimat);
        combs[0] = 0;
        combs[1] = 0;
        e = 0;
        h = 2; 
        Rboolean last = (nimplicants == 2);
        while (combs[0] != nimplicants - 2 || last) { 
            if (nimplicants == 2) {
                combs[1] = 1;
            }
            else {
                increment(2, &e, &h, nimplicants, combs, 0);
            }
            last = FALSE;
            int r = 0;
            int diffs = 0;
            int which = 0;
            Rboolean comparable = TRUE;
            while (diffs < 2 && r < nconds && comparable) {
                temp[r] = p_tempmat[combs[0] * nconds + r];
                if (temp[r] != p_tempmat[combs[1] * nconds + r]) {
                    comparable = temp[r] > 0 && p_tempmat[combs[1] * nconds + r] > 0;
                    diffs++;
                    which = r;
                    temp[r] = 0;
                }
                r++;
            }
            if (diffs == 1 && comparable) {
                int minrows[p_noflevels[which]];
                minrows[0] = combs[0];
                minrows[1] = combs[1];
                int tominimize = 2;
                int c = 0;
                while (c < nimplicants && tominimize < p_noflevels[which]) {
                    if (c != combs[0] && c != combs[1]) {
                        Rboolean equal = TRUE;
                        int r = 0;
                        while (r < nconds && equal) {
                            if (r != which) {
                                equal = temp[r] == p_tempmat[c * nconds + r];
                            }
                            r++;
                        }
                        if (equal) {
                            minrows[tominimize] = c;
                            tominimize++;
                        }
                    }
                    c++;
                }
                if (tominimize == p_noflevels[which]) {
                    for (int i = 0; i < tominimize; i++) {
                        p_minimized[minrows[i]] = TRUE;
                    }
                    int f = 0;
                    Rboolean Runique = TRUE;
                    while (f < found && Runique) {
                        Rboolean equal = TRUE;
                        int r = 0;
                        while (r < nconds && equal) {
                            equal = temp[r] == p_pimat[f * nconds + r];
                            r++;
                        }
                        Runique = !equal;
                        f++;
                    }
                    if (Runique) {
                        for (int r = 0; r < nconds; r++) {
                            p_pimat[found * nconds + r] = temp[r];
                        }
                        found++;
                    }
                }
            }
            if (found == estimpi) {
                estimpi *= 2;
                int totlent = found * nconds;
                SET_VECTOR_ELT(usage, 5, copymat = allocVector(INTSXP, totlent));
                p_copymat = INTEGER(copymat);
                for (int i = 0; i < totlent; i++) {
                    p_copymat[i] = p_pimat[i];
                }
                SET_VECTOR_ELT(usage, 4, pimat = allocMatrix(INTSXP, nconds, estimpi));
                p_pimat = INTEGER(pimat);
                for (int i = 0; i < totlent; i++) {
                    p_pimat[i] = p_copymat[i];
                }
            }
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
            p_copymat[c * nimplicants + r] = p_tempmat[p_order[r] * nconds + c]; 
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
        if (flag2 == 1) { 
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
                else { 
                    flag1 = 1;
                    pointer_next[previous] = pointer_next[i];
                    i = pointer_next[i];
                    j++;
                }
            }
        }
		rn = pointer_next[rn];
    }
    if (flag1 == 0) { 
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
    for (i = 0; i < lmbase; i++) {
        lmbasei = lmbase - i - 1;
        if (div(div(prowno[0] - 1, pmbase[lmbasei]).quot, pnoflevels[lmbasei] + 1).rem == 0) {
            flag = 1;
            lungime = templung * (pnoflevels[lmbasei] + 1);
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
    if (flag == 1) {
        templung = 0;
        for (i = 0; i < lungime; i++) {
            if (ptemp2[i] < (pmax[0] + 1)) {
                templung += 1;
            }
        }
        SET_VECTOR_ELT(usage, 4, temp1 = allocVector(INTSXP, templung - 1)); 
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
    SEXP pmin, max_ec;
    SEXP usage = PROTECT(allocVector(VECSXP, 5));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, y = coerceVector(y, REALSXP));
    double *p_x = REAL(x);
    double *p_y = REAL(y);
    int nrows_x = nrows(x);
    int ncols_x = ncols(x);
    SET_VECTOR_ELT(usage, 2, pmin = allocMatrix(REALSXP, nrows_x, ncols_x));
    double *p_pmin = REAL(pmin);
    SET_VECTOR_ELT(usage, 3, max_ec = allocMatrix(REALSXP, nrows_x, ncols_x));
    double *p_max_ec = REAL(max_ec);
    double sum_y = 0;
    double sum_x[ncols_x];
    double sum_neg_x[ncols_x];
    double sum_pmin[ncols_x];
    double sum_pmin_negy[ncols_x];
    double sum_neg_pmin[ncols_x];
    double sum_min_max_ec[ncols_x];
    for (int c = 0; c < ncols_x; c++) { 
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
    for (int r = 0; r < nrows_x; r++) {
        for (int c = 0; c < ncols_x - 1; c++) { 
            for (int cu = 0; cu < ncols_x - 1; cu++) {
                if (cu != c) {
                    if (p_max_ec[c * nrows_x + r] < p_pmin[cu * nrows_x + r]) {
                        p_max_ec[c * nrows_x + r] = p_pmin[cu * nrows_x + r];
                    }
                }
            }
        }
    }
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
            p_inclcov[c] = sum_pmin[c] / sum_y;
            p_inclcov[ncols_x + c] = sum_neg_x[c] / sum_neg_pmin[c];
            p_inclcov[2 * ncols_x + c] = sum_pmin[c] / sum_x[c];
            p_inclcov[3 * ncols_x + c] = 0; 
        }
        else {
            p_inclcov[c] = sum_pmin[c] / sum_x[c];
            p_inclcov[ncols_x + c] = (sum_pmin[c] - sum_pmin_negy[c]) / (sum_x[c] - sum_pmin_negy[c]);
            p_inclcov[2 * ncols_x + c] = sum_pmin[c] / sum_y;
            p_inclcov[3 * ncols_x + c] = p_inclcov[2 * ncols_x + c] - (sum_min_max_ec[c] / sum_y);
        }
    }
    UNPROTECT(1);
    return(inclcov);
}
SEXP C_ombnk(SEXP list) {
    int nconds, k, ogte, zerobased;
    nconds = INTEGER(VECTOR_ELT(list, 0))[0];
    k = INTEGER(VECTOR_ELT(list, 1))[0];
    ogte = INTEGER(VECTOR_ELT(list, 2))[0] - 1;
    zerobased = INTEGER(VECTOR_ELT(list, 3))[0];
    int nck = 1;
    for (int i = 1; i <= k; i++) {
        nck *= nconds - (k - i);
        nck /=  i;
    }
    SEXP usage = PROTECT(allocVector(VECSXP, 2));
    SEXP out;
    SET_VECTOR_ELT(usage, 0, out = allocMatrix(INTSXP, k, nck));
    int *p_out = INTEGER(out);
    int tempk[k];
    for (int i = 0; i < k; i++) {
        tempk[i] = i;
    }
    tempk[k - 1] -= 1;
    int e = 0;
    int h = k;
    Rboolean last = (k == nconds);
    int found = 0;
    while ((tempk[0] != nconds - k) || last) {
        increment(k, &e, &h, nconds + last, tempk, ogte);
        last = FALSE;
        for (int i = 0; i < k; i++) {
            p_out[found * k + i] = tempk[i] + 1 - zerobased;
        }
        found += 1;
    }
    if (found < nck) {
        SEXP copy;
        SET_VECTOR_ELT(usage, 1, copy = duplicate(out));
        int *p_copy = INTEGER(copy);
        SET_VECTOR_ELT(usage, 0, out = allocMatrix(INTSXP, k, found));
        p_out = INTEGER(out);
        for (int i = 0; i < found * k; i++) {
            p_out[i] = p_copy[i];
        }
    }
    UNPROTECT(1);
    return(out);
}
SEXP C_omplexity(SEXP list) {
    int nconds = INTEGER(VECTOR_ELT(list, 0))[0];
    int lk = length(VECTOR_ELT(list, 1));
    int *p_k = INTEGER(VECTOR_ELT(list, 1));
    int *p_noflevels = INTEGER(VECTOR_ELT(list, 2));
    SEXP result = PROTECT(allocVector(INTSXP, lk));
    int resum = 0;
    int prod = 1;
    for (int ck = 0; ck < lk; ck++) {
        resum = 0;
        int k = p_k[ck];
        int tempk[k];
        for (int i = 0; i < k; i++) {
            tempk[i] = i;
        }
        tempk[k - 1] -= 1;
        int e = 0;
        int h = k;
        Rboolean last = (nconds == k);
        while ((tempk[0] != nconds - k) || last) {
            increment(k, &e, &h, nconds + last, tempk, 0);
            last = FALSE;
            prod = 1;
            for (int i = 0; i < k; i++) {
                prod *= p_noflevels[tempk[i]];
            }
            resum += prod;
        }
        INTEGER(result)[ck] = resum;
    }
    UNPROTECT(1);
    return(result);
}
SEXP C_expand(SEXP mat, SEXP noflevels, SEXP partial) {
    int nconds = ncols(mat); 
    int nimp = nrows(mat); 
    int *p_imp = INTEGER(mat);
    int *p_noflevels = INTEGER(noflevels);
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
    SEXP x, rest, resmat, result, temp; 
    int rfilled = 0; 
    int estnrows = 1000; 
    SET_VECTOR_ELT(usage, 0, result = allocMatrix(INTSXP, nconds, estnrows));
    int *p_result = INTEGER(result);
    Memzero(p_result, nconds * estnrows);
    SET_VECTOR_ELT(usage, 1, x = allocMatrix(INTSXP, nimp, ncolsx));
    int *p_x = INTEGER(x);
    int pos = 0;
    int poszero = 0;
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
    int rmin = ncolsx;
    int wxim[nimp * ncolsx]; 
    int rxi[nimp]; 
    for (int r = 0; r < nimp; r++) {
        rxi[r] = 0;
        for (int c = 0; c < ncolsx; c++) {
            wxim[c * nimp + r] = -1; 
            if (p_x[c * nimp + r] == 0) {
                wxim[rxi[r] * nimp + r] = c;
                rxi[r] += 1;
            }
        }
        if (rxi[r] < rmin) {
            rmin = rxi[r];
        }
    }
    for (int ri = 0; ri < nimp; ri++) {
        int k = rxi[ri];
        if (LOGICAL(partial)[0]) {
            k = k - rmin;
        }
        if (k > 0) {
            int tempk[k];
            for (int i = 0; i < k; i++) {
                tempk[i] = i;
            }
            tempk[k - 1] -= 1;
            int e = 0;
            int h = k;
            Rboolean last = (k == rxi[ri]);
            while ((tempk[0] != rxi[ri] - k) || last) {
                increment(k, &e, &h, rxi[ri] + last, tempk, 0);
                last = FALSE;
                int nofl[k];
                for (int i = 0; i < k; i++) {
                   nofl[i] = p_noflevels[wxim[tempk[i] * nimp + ri]];
                }
                int nrows;
                calculate_rows(&nrows, k, nofl, 0, k);
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
                if ((rfilled + nrows) > estnrows) {
                    estnrows *= 2;
                    SET_VECTOR_ELT(usage, 3, result = Rresize(result, nconds * estnrows));
                    p_result = INTEGER(result);
                }
                for (int c = 0; c < ncolsx; c++) {
                    for (int r = 0; r < nrows; r++) {
                        p_result[(rfilled + r) * nconds + cnotzero[c]] = p_resmat[c * nrows + r];
                    }
                }
                rfilled += nrows;
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
        }
    }
    SET_VECTOR_ELT(usage, 4, temp = allocMatrix(INTSXP, nconds, rfilled));
    int *p_temp = INTEGER(temp);
    Memcpy(p_temp, p_result, nconds * rfilled);
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
    int posdata =      getpos(list, "data");       
    int posallsol =    getpos(list, "all.sol");    
    int posrowdom =    getpos(list, "row.dom");    
    int pospicons =    getpos(list, "pi.cons");    
    int posdepth =     getpos(list, "depth");      
    int posolcons =    getpos(list, "sol.cons");   
    int posolcov  =    getpos(list, "sol.cov");    
    int posfs =        getpos(list, "fs");         
    int posmaxcomb =   getpos(list, "max.comb");   
    int pos1stmin =    getpos(list, "first.min");  
    int poskeeptry =   getpos(list, "keep.trying");
    int posgurobi =    getpos(list, "gurobi");     
    int posolind =     getpos(list, "solind");     
    SEXP usage = PROTECT(allocVector(VECSXP, 7));
    SEXP   tt, data,    fsconds;
    SET_VECTOR_ELT(usage, 0, tt = coerceVector(VECTOR_ELT(list, 0), INTSXP));
    int *p_tt = INTEGER(tt);
    if (posdata > 0) {
        SET_VECTOR_ELT(usage, 1, data = coerceVector(VECTOR_ELT(list, posdata), REALSXP));
    }
    else {
        SET_VECTOR_ELT(usage, 1, data = allocMatrix(REALSXP, 2, 2)); 
        Memzero(REAL(data), 4);
    }
    double *p_data = REAL(data);
    int nrdata = nrows(data);
    int ttrows = nrows(tt); 
    int nconds = ncols(tt) - 1; 
    Rboolean allsol = (posallsol >= 0) ? (LOGICAL(VECTOR_ELT(list, posallsol))[0]) : FALSE;
    Rboolean rowdom = (posrowdom >= 0) ? (LOGICAL(VECTOR_ELT(list, posrowdom))[0]) : FALSE;
    Rboolean keeptrying = (poskeeptry >= 0) ? (LOGICAL(VECTOR_ELT(list, poskeeptry))[0]) : FALSE;
    Rboolean firstmin = (pos1stmin >= 0) ? (LOGICAL(VECTOR_ELT(list, pos1stmin))[0]) : FALSE;
    Rboolean gurobi = (posgurobi >= 0) ? (LOGICAL(VECTOR_ELT(list, posgurobi))[0]) : TRUE;
    Rboolean solind = (posolind >= 0) ? (LOGICAL(VECTOR_ELT(list, posolind))[0]) : TRUE;
    double picons = (pospicons >= 0) ? (REAL(VECTOR_ELT(list, pospicons))[0]) : 0;
    int pidepth = 0;
    int soldepth = 5; 
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
        SET_VECTOR_ELT(usage, 2, fsconds = allocVector(INTSXP, nconds)); 
        p_fsconds = INTEGER(fsconds);
        Memzero(p_fsconds, nconds);
    }
    int *p_pichart = R_Calloc(1, int); 
    int *p_impmat = R_Calloc(1, int);  
    int *p_models = R_Calloc(1, int);  
    unsigned int foundPI = 0;
    int solrows = 0;
    int solcols = 0;
    Rboolean complexpic = FALSE;
    double maxcomb = 0;
    if (posmaxcomb > 0) {
        maxcomb = REAL(VECTOR_ELT(list, posmaxcomb))[0]; 
    }
    CCubes(
        p_tt, ttrows, nconds, p_data, nrdata, allsol, rowdom, picons, pidepth, p_fsconds, soldepth, solcons, solcov, maxcomb, keeptrying,
        &p_pichart, &p_impmat, &p_models, &foundPI, &solrows, &solcols, &complexpic,
        firstmin,
        gurobi,
        solind
    );
    int posrows = 0;
    for (int r = 0; r < ttrows; r++) {
        posrows += p_tt[nconds * ttrows + r];
    }
    SEXP implicants, pic, solmat, complex;
    SEXP out = PROTECT(allocVector(VECSXP, 4));
    int i, j, l_1, len;
    if (foundPI > 0) {
        SET_VECTOR_ELT(out, 0, implicants = allocMatrix(INTSXP, foundPI, nconds));
        int *p_implicants = INTEGER(implicants);
        Memzero(p_implicants, foundPI * nconds); 
        len = foundPI * nconds;
        l_1 = len - 1;
        for (i = 0, j = 0; i < len; i++, j += nconds) {
            if (j > l_1) j -= l_1;
            p_implicants[i] = p_impmat[j];
        }
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
            SET_VECTOR_ELT(dimnames, 1, colnms); 
            setAttrib(implicants, R_DimNamesSymbol, dimnames);
        }
    }
    if (solrows > 0 && solcols > 0) { 
        SET_VECTOR_ELT(out, 2, solmat = allocMatrix(INTSXP, solrows, solcols));
        int *p_solmat = INTEGER(solmat);
        if (firstmin) {
            for (int i = 0; i < solrows; i++) {
                p_solmat[i] = p_models[i] + 1;
            }
        }
        else {
            Memcpy(p_solmat, p_models, solrows * solcols);
        }
    }
    SET_VECTOR_ELT(out, 3, complex = allocVector(LGLSXP, 1));
    LOGICAL(complex)[0] = complexpic;
    R_Free(p_pichart);
    R_Free(p_impmat);
    R_Free(p_models);
    UNPROTECT(2);
    return(out);
}
SEXP C_getEC(SEXP dem, SEXP cexpr, SEXP csolm, SEXP pexpr, SEXP psolm, SEXP SA, SEXP noflevels) {
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
    for (int r = 0; r < nrows_dem; r++) {
        rowsums[r] = 0;
        for (int c = 0; c < ncols_dem; c++) {
            if (p_dem[c * nrows_dem + r] > 0) {
                rowsums[r] += 1;
            }
        }
        if (rowsums[r] == 1) {
            idecount += 1;
        } else if (rowsums[r] > 1) {
            cdecount += 1;
        }
    }
    SEXP iderows, cderows;
    SET_VECTOR_ELT(usage, 1, iderows = allocVector(INTSXP, idecount));
    int *p_iderows = INTEGER(iderows);
    SET_VECTOR_ELT(usage, 2, cderows = allocVector(INTSXP, cdecount));
    int *p_cderows = INTEGER(cderows);
    for (int i = 0; i < idecount; i++) {
        p_iderows[i] = 0; 
    }
    for (int i = 0; i < cdecount; i++) {
        p_cderows[i] = 0; 
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
    SEXP sap; 
    int *p_sap;
    int tempcount = 0;
    int ECpos = 0;
    for (int csol = 0; csol < ncols_csolm; csol++) { 
        int csolcount = 0;
        for (int r = 0; r < nrows_csolm; r++) {
            if (p_csolm[csol * nrows_csolm + r] > 0) {
                csolcount++;
            }
        }
        for (int psol = 0; psol < ncols_psolm; psol++) { 
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
                subset[ics] = FALSE; 
                for (int c = 0; c < ncols_dem; c++) {
                    comp[c] = p_cexpr[c * nrows_cexpr + p_csolm[csol * nrows_csolm + ics] - 1];
                }
                for (int ips = 0; ips < psolcount; ips++) {
                    Rboolean allequal = TRUE;
                    Rboolean parsgz[ncols_dem]; 
                    Rboolean compgz[ncols_dem]; 
                    int c = 0;
                    while (allequal && c < ncols_dem) {
                        pars[c] = p_pexpr[c * nrows_pexpr + p_psolm[psol * nrows_psolm + ips] - 1];
                        parsgz[c] = pars[c] > 0;
                        compgz[c] = comp[c] > 0 && !parsgz[c]; 
                        if (parsgz[c]) {
                            allequal = pars[c] == comp[c];
                        }
                        c++;
                    }
                    if (allequal) { 
                        subset[ics] = TRUE;
                        if (idecount > 0) { 
                            for (int c = 0; c < ncols_dem; c++) {
                                res[c] = comp[c]; 
                                if (compgz[c]) {
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
                                        res[c] = 0;
                                    }
                                }
                            }
                        }
                        if (cdecount > 0) { 
                            for (int r = 0; r < cdecount; r++) {
                                Rboolean allcde = TRUE;
                                for (int c = 0; c < ncols_dem; c++) {
                                    res[c] = comp[c];
                                    if (compgz[c]) { 
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
                        for (int r = 0; r < ncols_dem; r++) { 
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
                } 
            } 
            int ECrows = 0;
            for (int r = 0; r < nrows_sap; r++) {
                ECrows += ecs[r] * 1;
            }
            int countnot = 0; 
            for (int c = 0; c < csolcount; c++) {
                if (!subset[c]) {
                    countnot++;
                }
            }
            SEXP intseltemp;
            int totalrows = tempcount + countnot;
            SET_VECTOR_ELT(usage, 7, intseltemp = allocMatrix(INTSXP, totalrows, ncols_dem));
            p_intseltemp = INTEGER(intseltemp);
            Memzero(p_intseltemp, totalrows * ncols_dem);
            for (int c = 0; c < tempcount; c++) {
                for (int r = 0; r < ncols_dem; r++) {
                    p_intseltemp[r * totalrows + c] = p_tempmat[c * ncols_dem + r];
                }
            }
            if (countnot > 0) { 
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
            SEXP pi;
            SEXP intselmat;
            SEXP partial;
            SET_VECTOR_ELT(usage, 4, partial = allocVector(LGLSXP, 1));
            LOGICAL(partial)[0] = FALSE;
            SET_VECTOR_ELT(intsel, ECpos, intselmat = intseltemp);
            SET_VECTOR_ELT(primes, ECpos, pi = C_simplify(intseltemp, noflevels, partial));
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
            for (int r = 0; r < nrows_sap; r++) {
                if (ecs[r]) {
                    for (int c = 0; c < ncols_dem; c++) {
                        p_ec[c * ECrows + ecr] = p_sap[c * nrows_sap + r];
                    }
                    SET_STRING_ELT(rownms, ecr, STRING_ELT(sapnms, r));
                    ecr++;
                }
            }
            if (ecr > 0) {
                SET_VECTOR_ELT(dim_names, 0, duplicate(rownms));
                if (hasColnames(cexpr)) {
                    SET_VECTOR_ELT(dim_names, 1, VECTOR_ELT(getAttrib(cexpr, R_DimNamesSymbol), 1));
                }
                setAttrib(ec, R_DimNamesSymbol, dim_names);
            }
            ECpos++;
        } 
    } 
    SET_VECTOR_ELT(output, 0, EClist);
    SET_VECTOR_ELT(output, 1, primes);
    SET_VECTOR_ELT(output, 2, intsel);
    UNPROTECT(5);
    return(output);
}
