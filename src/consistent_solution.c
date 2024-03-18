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
#include <math.h>
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
    double *p_y = R_Calloc(1, double);
    double ymat[nrdata * k];
    double sumy = 0;
    for (int r = 0; r < nrdata; r++) {
        sumy += p_data[nconds * nrdata + r];
    }
    for (int i = 0; i < k; i++) { 
        int k2 = ck[tempk[i]];
        R_Free(p_y);
        p_y = R_Calloc(nrdata * k2, double);
        for (int c = 0; c < k2; c++) {
            cindx = indx[tempk[i] * nconds + c] - 1;
            val = p_implicants[tempk[i] * nconds + cindx] - 1;
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
    R_Free(p_y);
    return(agteb(sumxy / sumx, solcons) && agteb(sumxy / sumy, solcov));
}
