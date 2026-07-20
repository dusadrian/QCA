#include "qca_rinternals.h"
#include "qca_rinternals.h"
# include <R_ext/Rdynload.h>


SEXP C_truthTable(SEXP x, SEXP vo, SEXP tt, SEXP fuz) {
    // x is the initial data conditions (might be fuzzy as well)
    // y is the truth table (crisp)
    // fuz is signaling which conditions are fuzzy
    // vo is the initial outcome column (might be fuzzy as well)
    int i, j, k, index;
    double *p_x, *p_inclpri, *p_vo, min, sumx, sumpmin, prisum, temp1, temp2; // , so
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
    
    
    // so = 0;
    
    // // sum of the outcome variable
    // for (i = 0; i < length(vo); i++) {
    //     so += p_vo[i];
    // }
    
    
    for (k = 0; k < ttrows; k++) { // loop for every line of the truth table matrix
        
        sumx = 0;
        sumpmin = 0;
        prisum = 0;  
        ncut = 0;
        
        for (i = 0; i < xrows; i++) { // loop over every line of the data matrix
            
            min = 1000;
            for (j = 0; j < xcols; j++) { // loop over each column of the data matrix
                copyline[j] = p_x[j * xrows + i];
                
                index = j * ttrows + k;
                
                if (p_fuz[j] == 1) { // for the fuzzy variables, invert those who have the 3k value equal to 1 ("onex3k" in R)
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
                
            } // end of j loop, over columns
            
            sumx += min;
            sumpmin += (min < p_vo[i])?min:p_vo[i];
            temp1 = (min < p_vo[i])?min:p_vo[i];
            temp2 = 1 - p_vo[i];
            prisum += (temp1 < temp2)?temp1:temp2;
            ncut += (min > 0.5)?1:0;
            
            p_inclpri[k*(3 + xrows) + 3 + i] = min;
            
        } // end of i loop
        
        
        p_inclpri[k*(3 + xrows) + 0] = ncut;
        p_inclpri[k*(3 + xrows) + 1] = sumpmin/sumx;
        p_inclpri[k*(3 + xrows) + 2] = (sumpmin - prisum)/(sumx - prisum);
        
    } // end of k loop
    
    UNPROTECT(2);
    
    return(inclpri);
}
