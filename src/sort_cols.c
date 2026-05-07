#include "qca_r.h"
#include "sort_cols.h"


void sort_cols(
    int *p_matrix,
    int *sortcols,
    int *p_ck,
    const int nconds,
    unsigned int foundPI
) {
    
    /*
    // simulate this:
    sortExpressions <- function(mat) {
        
        for (i in rev(seq(ncol(mat)))) {
            mat <- mat[order(mat[, i], decreasing = TRUE), , drop = FALSE]
            if (length(wx <- which(mat[, i] > 0)) > 0) {
                rest <- if (max(wx) == nrow(mat)) NULL else seq(max(wx) + 1, nrow(mat))
                mat <- mat[c(order(mat[wx, i]), rest), , drop = FALSE]
            }
        }
        
        return(mat[order(apply(mat, 1, function(x) sum(x > 0))), , drop = FALSE])
    }
    */
    
    int temp;
    
    for (int i = nconds - 1; i >= 0; i--) {
        for (unsigned int c1 = 0; c1 < foundPI; c1++) {
            for (unsigned int c2 = c1 + 1; c2 < foundPI; c2++) {
                if (p_matrix[sortcols[c1] * nconds + i] < p_matrix[sortcols[c2] * nconds + i]) {
                    
                    temp = sortcols[c2];
                    
                    for (int c3 = c2; c3 > c1; c3--) {
                        sortcols[c3] = sortcols[c3 - 1];
                    }
                    
                    sortcols[c1] = temp;
                }
            }
        }
        
        Rboolean nonzero = true;
        unsigned int zeroidx = 0;
        while(zeroidx < foundPI && nonzero) {
            
            nonzero = p_matrix[sortcols[zeroidx] * nconds + i];
            zeroidx++;
            
        }
        
        zeroidx--;
        
        for (unsigned int c1 = 0; c1 < zeroidx; c1++) {
            for (unsigned int c2 = c1 + 1; c2 < zeroidx; c2++) {
                if (p_matrix[sortcols[c1] * nconds + i] > p_matrix[sortcols[c2] * nconds + i]) {
                    
                    temp = sortcols[c2];
                    
                    for (int c3 = c2; c3 > c1; c3--) {
                        sortcols[c3] = sortcols[c3 - 1];
                    }
                    
                    sortcols[c1] = temp;
                }
            }
        }
    }
    
    for (unsigned int c1 = 0; c1 < foundPI; c1++) {
        for (unsigned int c2 = c1 + 1; c2 < foundPI; c2++) {
            if (p_ck[sortcols[c1]] > p_ck[sortcols[c2]]) {
                
                temp = sortcols[c2];
                    
                for (int c3 = c2; c3 > c1; c3--) {
                    sortcols[c3] = sortcols[c3 - 1];
                }
                
                sortcols[c1] = temp;
            }
        }
    }
}
