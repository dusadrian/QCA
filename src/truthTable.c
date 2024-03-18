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

# include <R.h>
# include <Rinternals.h>
# include <R_ext/Rdynload.h>
SEXP C_truthTable(SEXP x, SEXP vo, SEXP tt, SEXP fuz) {
    int i, j, k, index;
    double *p_x, *p_inclpri, *p_vo, min, sumx, sumpmin, prisum, temp1, temp2; 
    int xrows, xcols, ttrows, ncut, *p_tt, *p_fuz;
    SEXP usage = PROTECT(allocVector(VECSXP, 4));
    SET_VECTOR_ELT(usage, 0, x = coerceVector(x, REALSXP));
    SET_VECTOR_ELT(usage, 1, vo = coerceVector(vo, REALSXP));
    SET_VECTOR_ELT(usage, 2, tt = coerceVector(tt, INTSXP));
    SET_VECTOR_ELT(usage, 3, fuz = coerceVector(fuz, INTSXP));
    xrows = nrows(x);
    xcols = ncols(x);
    ttrows = nrows(tt);
    double copyline[xcols];
    p_x = REAL(x);
    p_vo = REAL(vo);
    p_tt = INTEGER(tt);
    p_fuz = INTEGER(fuz);
    SEXP inclpri = PROTECT(allocMatrix(REALSXP, 3 + xrows, ttrows));
    p_inclpri = REAL(inclpri);
    for (k = 0; k < ttrows; k++) { 
        sumx = 0;
        sumpmin = 0;
        prisum = 0;  
        ncut = 0;
        for (i = 0; i < xrows; i++) { 
            min = 1000;
            for (j = 0; j < xcols; j++) { 
                copyline[j] = p_x[j * xrows + i];
                index = j * ttrows + k;
                if (p_fuz[j] == 1) { 
                    if (p_tt[index] == 0) {
                        copyline[j] = 1 - copyline[j];
                    }
                }
                else {
                    if (p_tt[index] != (copyline[j])) {
                        copyline[j] = 0;
                    }
                    else {
                        copyline[j] = 1;
                    }
                }
                if (copyline[j] < min) {
                    min = copyline[j];
                }
            } 
            sumx += min;
            sumpmin += (min < p_vo[i])?min:p_vo[i];
            temp1 = (min < p_vo[i])?min:p_vo[i];
            temp2 = 1 - p_vo[i];
            prisum += (temp1 < temp2)?temp1:temp2;
            ncut += (min > 0.5)?1:0;
            p_inclpri[k*(3 + xrows) + 3 + i] = min;
        } 
        p_inclpri[k*(3 + xrows) + 0] = ncut;
        p_inclpri[k*(3 + xrows) + 1] = sumpmin/sumx;
        p_inclpri[k*(3 + xrows) + 2] = (sumpmin - prisum)/(sumx - prisum);
    } 
    UNPROTECT(2);
    return(inclpri);
}
