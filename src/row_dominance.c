
#include "qca_r.h"
#include "row_dominance.h"


void row_dominance(
    int p_pichart[],
    int p_implicants[],
    int *p_ck,
    int pirows,
    int *foundPI,
    int nconds
) {
    
    unsigned int picols = *foundPI;


    // int* survcols = (int *) R_Calloc (picols, int);
    
    // for (int i = 0; i < picols; i++) {
    //     survcols[i] = true;
    // }
    
    Rboolean survcols[picols];
    int colsums[picols];
    int sortcol[picols];
    int temp;
    
    for (unsigned int c = 0; c < picols; c++) {
        colsums[c] = 0;
        
        for (int r = 0; r < pirows; r++) {
            colsums[c] += p_pichart[c * pirows + r];
        }
        
        sortcol[c] = c;
        survcols[c] = true;
    }
    
    for (unsigned int c1 = 0; c1 < picols; c1++) {
        for (unsigned int c2 = c1 + 1; c2 < picols; c2++) {
            if (colsums[sortcol[c1]] < colsums[sortcol[c2]]) {
                temp = sortcol[c1];
                sortcol[c1] = sortcol[c2];
                sortcol[c2] = temp;
            }
        }
    }

    for (unsigned int c1 = 0; c1 < picols; c1++) {
        if (survcols[sortcol[c1]]) {
            for (unsigned int c2 = c1 + 1; c2 < picols; c2++) {
                if (survcols[sortcol[c2]]) {
                    if (colsums[sortcol[c1]] > colsums[sortcol[c2]]) {
                        
                        Rboolean itcovers = true; // assume it's covered
                        int r = 0;
                        
                        while (r < pirows && itcovers) {
                            if (p_pichart[sortcol[c2] * pirows + r]) {
                                itcovers = p_pichart[sortcol[c1] * pirows + r];
                            }
                            r++;
                        }
                        
                        if (itcovers) {
                            survcols[sortcol[c2]] = false;
                            --(*foundPI);
                        }
                    }
                }
            }
        }
    }
        
    
    if (*foundPI < picols) {

        // move (overwrite) all surviving columns towards the beginning
        int s = 0;
        for (unsigned int c = 0; c < picols; c++) {
            if (survcols[c]) {
                for (int r = 0; r < pirows; r++) {
                    p_pichart[s * pirows + r] = p_pichart[c * pirows + r];
                }

                for (int r = 0; r < nconds; r++) {
                    p_implicants[s * nconds + r] = p_implicants[c * nconds + r];
                }
                
                // same with the vector positions
                p_ck[s] = p_ck[c];

                s++;
            }
        }
    }

    return;
}
