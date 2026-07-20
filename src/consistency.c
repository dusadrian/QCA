#include <R_ext/RS.h> // for R_Calloc, R_free, Memset etc.
#include "qca_r.h"
#include "utils.h"


double consistency(
    const double p_x[],
    const int nrowsx,
    const int nconds,
    int k,
    int tempk[],
    int val[],
    int fuzzy[]
) {
    
    double *p_y = (double *) R_Calloc(nrowsx * k, double);
    // memset(p_y, 0, nrowsx * k * sizeof(double));
    
    
    for (int c = 0; c < k; c++) {
        
        if (fuzzy[c]) {
            Rboolean negation = val[c] == 0;
            for (int r = 0; r < nrowsx; r++) {
                p_y[c * nrowsx + r] = negation ? (1 - p_x[tempk[c] * nrowsx + r]) : p_x[tempk[c] * nrowsx + r];
            }
        }
        else {
            for (int r = 0; r < nrowsx; r++) {
                p_y[c * nrowsx + r] = (p_x[tempk[c] * nrowsx + r] == val[c]) ? 1 : 0;
            }
        }
    }
    
    
    double pminx;
    double sumx = 0, sumxy = 0;
    

    for (int r = 0; r < nrowsx; r++) {
        pminx = 1;
        for (int c = 0; c < k; c++) {
            if (p_y[c * nrowsx + r] < pminx) {
                pminx = p_y[c * nrowsx + r];
            }
        }
        
        sumx += pminx;
        sumxy += ((pminx < p_x[nconds * nrowsx + r]) ? pminx : p_x[nconds * nrowsx + r]);
        
    }
    
    R_Free(p_y);
    return(sumxy / sumx);
}
