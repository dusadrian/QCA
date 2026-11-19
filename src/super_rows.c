#include "qca_r.h"
#include "super_rows.h"

void super_rows(
    int *p_matrix,
    int rows,
    unsigned int *survcols,
    int *p_cols,
    int depth
) {
    
    unsigned int cols = *survcols;
    int colsums[cols];
    int sortcol[cols];
    int temp;
    unsigned int eligible = 0;
    
    for (unsigned int c = 0; c < cols; c++) {
        colsums[c] = 0;
        
        for (int r = 0; r < rows; r++) {
            colsums[c] += p_matrix[c * rows + r];
        }
        
        // Adding PIs for later configurations can only increase this size.
        if (depth > 0 && colsums[c] > depth) {
            p_cols[c] = false;
            --(*survcols);
        } else {
            sortcol[eligible++] = c;
        }
    }
    
    
    for (unsigned int c1 = 0; c1 + 1 < eligible; c1++) {
        for (unsigned int c2 = c1 + 1; c2 < eligible; c2++) {
            if (colsums[sortcol[c1]] > colsums[sortcol[c2]]) {
                temp = sortcol[c1];
                sortcol[c1] = sortcol[c2];
                sortcol[c2] = temp;
            }
        }
    }
    
    for (unsigned int c1 = 0; c1 + 1 < eligible; c1++) {
        if (p_cols[sortcol[c1]]) {
            for (unsigned int c2 = c1 + 1; c2 < eligible; c2++) {
                if (p_cols[sortcol[c2]]) {
                    if (colsums[sortcol[c1]] <= colsums[sortcol[c2]]) {
                        
                        Rboolean cover = true; // assume it covers
                        int r = 0;
                        
                        while (r < rows && cover) {
                            if (p_matrix[sortcol[c1] * rows + r]) {
                                cover = p_matrix[sortcol[c2] * rows + r];
                            }
                            r++;
                        }
                        
                        
                        if (cover) {
                            p_cols[sortcol[c2]] = false;
                            --(*survcols);
                        }
                    }
                }
            }
        }
    }
    
    if (*survcols < cols) {
        unsigned int s = 0;
        for (unsigned int c = 0; c < cols; c++) {
            if (p_cols[c]) {
                for (int r = 0; r < rows; r++) {
                    p_matrix[s * rows + r] = p_matrix[c * rows + r];
                }
                s++;
            }
        }
    }
}
