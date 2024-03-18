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

#include <R_ext/RS.h> 
#include <R_ext/Boolean.h>
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
