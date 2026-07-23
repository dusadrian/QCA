#include <stdlib.h>
#include <math.h>
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
    
    double *p_y = (double *) calloc((size_t)nrowsx * (size_t)k, sizeof(double));
    if (p_y == NULL) {
        return NAN;
    }
    
    
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
    
    free(p_y);
    return(sumxy / sumx);
}
