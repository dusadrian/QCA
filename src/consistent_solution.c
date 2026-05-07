
#include "qca_r.h"
#include <R_ext/RS.h>
#include "qca_r.h"
#include <math.h>
#include <stdlib.h>
#include "utils.h"
#include "consistent_solution.h"



Rboolean consistent_solution(
    const double p_data[],
    const int nconds,
    const int nrdata,
    const int k,
    const int tempk[],
    const unsigned int foundPI,
    const int p_implicants[],
    const int ck[],
    const int indx[],
    const int p_fsconds[],
    const double solcons,
    const double solcov) {
    
    
    int cindx, val;
    double *p_y = (double *) calloc(1, sizeof(double));
    double *ymat = (double *) calloc((size_t) nrdata * (size_t) k, sizeof(double));
    if (p_y == NULL || ymat == NULL) {
        free(p_y);
        free(ymat);
        return FALSE;
    }
    
    double sumy = 0;
    for (int r = 0; r < nrdata; r++) {
        sumy += p_data[nconds * nrdata + r];
    }
    
    // Rboolean aici = false;
    // if (k == 2) {
    //     aici = (tempk[0] == 0 && tempk[1] == 2) || (tempk[0] == 1 && tempk[1] == 3);
    // }

    for (int i = 0; i < k; i++) { // for each conjunction
        
        int k2 = ck[tempk[i]];
        if (k2 <= 0) {
            free(p_y);
            free(ymat);
            return FALSE;
        }
        free(p_y);
        p_y = (double *) calloc((size_t) nrdata * (size_t) k2, sizeof(double));
        if (p_y == NULL) {
            free(ymat);
            return FALSE;
        }
        
        
        for (int c = 0; c < k2; c++) {
            // cindx = indx[c * foundPI + tempk[i]];
            cindx = indx[tempk[i] * nconds + c] - 1;
            if (cindx < 0 || cindx >= nconds) {
                free(p_y);
                free(ymat);
                return FALSE;
            }
            
            // val = p_implicants[cindx * foundPI + tempk[i]] - 1;
            val = p_implicants[tempk[i] * nconds + cindx] - 1;
            
            // if (aici) {
            //     printf("%d%d; cindx: %d, val: %d\n", tempk[0],tempk[1], cindx, val);
            // }

            if (p_fsconds[cindx]) {
                
                Rboolean negation = val == 0;
                for (int r = 0; r < nrdata; r++) {
                    p_y[c * nrdata + r] = negation ? (1 - p_data[cindx * nrdata + r]) : p_data[cindx * nrdata + r];
                }
            }
            else {
                for (int r = 0; r < nrdata; r++) {
                    p_y[c * nrdata + r] = (p_data[cindx * nrdata + r] == val) ? 1 : 0;
                }
            }
        }
        
        double pminx;
        
        for (int r = 0; r < nrdata; r++) {
            pminx = 1;
            
            for (int c = 0; c < k2; c++) {
                if (p_y[c * nrdata + r] < pminx) {
                    pminx = p_y[c * nrdata + r];
                }
            }
            
            ymat[i * nrdata + r] = pminx;
        }
        
    }

    // if (aici) {
    //     for (int r = 0; r < nrdata; r++) {
    //         for (int c = 0; c < k; c++) {
    //             printf("%1.0f ", ymat[c * nrdata + r]);
    //         }
    //         printf("\n");
    //     }
    //     printf("\n");
    // }
    
    double pmaxx;
    double sumx = 0, sumxy = 0;
    
    for (int r = 0; r < nrdata; r++) {
        pmaxx = 0;
        for (int c = 0; c < k; c++) {
            if (ymat[c * nrdata + r] > pmaxx) {
                pmaxx = ymat[c * nrdata + r];
            }
        }
        
        sumx += pmaxx;
        sumxy += ((pmaxx < p_data[nconds * nrdata + r]) ? pmaxx: p_data[nconds * nrdata + r]);
        
    }
    
    // Guard against divisions by zero (e.g., unattainable solutions with solcov = 1)
    if (sumx <= 0 || sumy <= 0) {
        free(p_y);
        free(ymat);
        return FALSE;
    }
    
    // if (aici) {
        // Rprintf("sumxy: %5.3f; sumx: %5.3f; sumy: %5.3f, solcons: %5.3f, solcov: %5.3f\n", sumxy, sumx, sumy, solcons, solcov);
        // Rprintf(
        //     "incl: %5.3f; cov: %5.3f; decision: %d\n",
        //     sumxy / sumx,
        //     sumxy / sumy,
        //     agteb(sumxy / sumx, solcons) && agteb(sumxy / sumy, solcov)
        // );
    // }
    free(p_y);
    free(ymat);
    return(agteb(sumxy / sumx, solcons) && agteb(sumxy / sumy, solcov));
}
