
#include "generate_matrix.h"

void generate_matrix(int nrows, int ncols, int nofl[], int arrange, int maxprod, int *p_matrix) {
    
    int e, h, k, prod;
    
    if (arrange == 0) {
        int cols[ncols];
        for (int c = 0; c < ncols; c++) {
            cols[c] = c;
        }
        
        fill_matrix(nrows, ncols, nofl, p_matrix, 0, cols, 0);
    }
    else { // p_arrange[0] == 1
        
        for (int i = 0; i < nrows * ncols; i++) {
            p_matrix[i] = 0;
        }
        
        
        // start filling the matrix, one section of rows after another
        int startrow = 0;
        
        for (k = 1; k <= maxprod; k++) {
            
            int tempk[k];
            
            int nck;
            nck = 1;
            for (int i = 1; i <= k; i++) {
                nck *= ncols - (k - i);
                nck /=  i;
            }
            
            for (int i = 0; i < k; i++) {
                tempk[i] = i;
            }
            
            e = 0;
            h = k;
            
            for (int count = 0; count < nck; count++) {
                
                if (count > 0) {
                    increment(k, &e, &h, ncols, tempk, 0);
                }
                
                prod = 1;
                int colsk[k];
                int noflk[k];
                for (int c = 0; c < k; c++) {
                    prod *= (nofl[tempk[c]] - 1);
                    colsk[c] = tempk[c];
                    noflk[c] = nofl[tempk[c]] - 1;
                }
                
                fill_matrix(nrows, k, noflk, p_matrix, startrow, colsk, 1);
                
                startrow += prod;
            }
        }
    }
}
