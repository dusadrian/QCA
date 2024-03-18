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
#include "find_consistent_models.h"
#include "utils.h"
void find_consistent_models(
    const int p_implicants[],
    const int indx[],
    const int ck[],
    const double p_data[],
    const int p_fuzzy[],
    const int nconds,
    const int nrdata,
    const int posrows,
    const double solcons,
    const double solcov,
    const Rboolean allsol,
    const int soldepth,
    const unsigned int foundPI,
    const double maxcomb,
    int **solutions,
    int *nr,
    int *nc
) {
    unsigned int estimsol = 1000;
    int maxk = posrows;
    if (foundPI < maxk) {
        maxk = foundPI;
    }
    if (soldepth < maxk && soldepth > 0) {
        maxk = soldepth;
    }
    int *p_sol = R_Calloc(maxk * estimsol, int);
    int *cksol = R_Calloc(estimsol, int);
    unsigned int solfound = 0;
    unsigned int prevfound = 0;
    Rboolean keep_searching = true;
    int k = 1;
    double counter = 1;
    while (keep_searching && k <= maxk) {
            int tempk[k];
            for (int i = 0; i < k; i++) {
                tempk[i] = i; 
            }
            tempk[k - 1] -= 1; 
            int e = 0;
            int h = k;
            Rboolean last = (foundPI == k);
            while (keep_searching && ((tempk[0] != foundPI - k) || last)) {
                increment(k, &e, &h, foundPI + last, tempk, 0);
                last = false;
                Rboolean nonred = true;
                int i = 0;
                while (i < prevfound && nonred) {
                int sumeq = 0;
                int v = 0;
                while (sumeq == v && v < cksol[i]) {
                    for (int c = 0; c < k; c++) {
                        if (p_sol[i * maxk + v] == tempk[c]) {
                            sumeq++;
                        }
                    }
                    v++;
                }
                if (sumeq == v) { 
                    nonred = false; 
                }
                i++;
            }
            if (nonred) {
                if (consistent_solution(p_data, nconds, nrdata, k, tempk, foundPI, p_implicants, ck, indx, p_fuzzy, solcons, solcov)) {
                    for (int c = 0; c < k; c++) {
                        p_sol[solfound * maxk + c] = tempk[c];
                    }
                    cksol[solfound] = k;
                    solfound++;
                    if (solfound == estimsol) {
                        estimsol += 1000;
                        resize(&p_sol,  maxk,  estimsol, solfound);
                        resize(&cksol,   1,  estimsol, solfound);
                    }
                }
            }
            if (maxcomb > 0) {
                counter++;
                if ((counter / 1000000000) >= maxcomb) {
                    keep_searching = false;
                }
            }
        }
        prevfound = solfound;
        k += 1;
    }
    int *p_tempmat = R_Calloc(1, int);
    if (solfound > 0) {
        int finalrows = cksol[solfound - 1];
        R_Free(p_tempmat);
        p_tempmat = R_Calloc(finalrows * solfound, int);
        for (int c = 0; c < solfound; c++) {
            for (int r = 0; r < cksol[c]; r++) {
                p_tempmat[c * finalrows + r] = p_sol[c * maxk + r] + 1; 
            }
        }
        *nr = finalrows;
        *nc = solfound;
    }
    R_Free(p_sol);
    R_Free(cksol);
    R_Free(*solutions);
    *solutions = p_tempmat;
}
