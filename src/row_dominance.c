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
#include "row_dominance.h"
void row_dominance(
    int p_pichart[],
    int p_implicants[],
    int *p_ck,
    int pirows,
    unsigned int *foundPI,
    int nconds
) {
    unsigned int picols = *foundPI;
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
                        Rboolean itcovers = true; 
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
        int s = 0;
        for (unsigned int c = 0; c < picols; c++) {
            if (survcols[c]) {
                for (int r = 0; r < pirows; r++) {
                    p_pichart[s * pirows + r] = p_pichart[c * pirows + r];
                }
                for (int r = 0; r < nconds; r++) {
                    p_implicants[s * nconds + r] = p_implicants[c * nconds + r];
                }
                p_ck[s] = p_ck[c];
                s++;
            }
        }
    }
    return;
}
