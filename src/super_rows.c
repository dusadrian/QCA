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

#include <R_ext/Boolean.h>
#include "super_rows.h"
void super_rows(
    int *p_matrix,
    int rows,
    unsigned int *survcols,
    int *p_cols
) {
    unsigned int cols = *survcols;
    int colsums[cols];
    int sortcol[cols];
    int temp;
    for (unsigned int c = 0; c < cols; c++) {
        colsums[c] = 0;
        for (int r = 0; r < rows; r++) {
            colsums[c] += p_matrix[c * rows + r];
        }
        sortcol[c] = c;
    }
    for (unsigned int c1 = 0; c1 < cols; c1++) {
        for (unsigned int c2 = c1 + 1; c2 < cols; c2++) {
            if (colsums[sortcol[c1]] > colsums[sortcol[c2]]) {
                temp = sortcol[c1];
                sortcol[c1] = sortcol[c2];
                sortcol[c2] = temp;
            }
        }
    }
    for (unsigned int c1 = 0; c1 < cols - 1; c1++) {
        if (p_cols[sortcol[c1]]) {
            for (unsigned int c2 = c1 + 1; c2 < cols; c2++) {
                if (p_cols[sortcol[c2]]) {
                    if (colsums[sortcol[c1]] <= colsums[sortcol[c2]]) {
                        Rboolean cover = true; 
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
