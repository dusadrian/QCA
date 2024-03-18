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
    else { 
        for (int i = 0; i < nrows * ncols; i++) {
            p_matrix[i] = 0;
        }
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
