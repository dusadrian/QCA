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
#include <float.h> 
#include <math.h>
#include "utils.h"
void resize(
    int **array,
    const int rows,
    const unsigned int newcols,
    const unsigned int oldcols
) {
    int* tmp = R_Calloc(rows * newcols, int);
    unsigned int cols = (newcols > oldcols) ? oldcols : newcols;
    for (unsigned int i = 0; i < rows * cols; i++) {
        *(tmp + i) = *(*array + i);
    }
    R_Free(*array);
    *array = tmp;
}
Rboolean too_complex(
    const unsigned int foundPI,
    const int solmin,
    const double maxcomb
) {
    double result = 1;
    unsigned int n = foundPI;
    int k = solmin;
    for (int i = 1; i <= k; i++) {
        result *= n - (k - i);
        result /= i;
    }
    if ((result / 1000000000) > maxcomb && maxcomb > 0) {
        return(true);
    }
    return(false);
}
void over_transpose(
    int matrix[],
    const int nr,
    const int nc,
    const int type
) {
    int len = nr * nc;
    int i, j, l_1 = len - 1;
    if (type == 0) {
        Rboolean tmp[nr * nc];
        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            tmp[i] = matrix[j];
        }
        for (int i = 0; i < len; i++) {
            matrix[i] = tmp[i];
        }
    }
    else if (type == 1) {
        int tmp[nr * nc];
        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            tmp[i] = matrix[j];
        }
        for (int i = 0; i < len; i++) {
            matrix[i] = tmp[i];
        }
    }
    else if (type == 2) {
        double tmp[nr * nc];
        for (i = 0, j = 0; i < len; i++, j += nr) {
            if (j > l_1) j -= l_1;
            tmp[i] = matrix[j];
        }
        for (int i = 0; i < len; i++) {
            matrix[i] = tmp[i];
        }
    }
}
Rboolean altb(double a, double b) {
    return (b - a) > ( (fabs(a) < fabs(b) ? fabs(b) : fabs(a)) * DBL_EPSILON);
}
Rboolean agteb(double a, double b) {
    return((a > b) || (fabs(a - b) <= DBL_EPSILON));
}
Rboolean nonredundant(
    const int p_implicants[],
    const int p_indx[],
    const int p_ck[],
    const int tempk[],
    const int tempc[],
    const int nconds,
    const int k,
    const unsigned int prevfoundPI
) {
    Rboolean nonred = true;
    if (prevfoundPI > 0) {
        unsigned int i = 0;
        while (i < prevfoundPI && nonred) {
            int sumeq = 0;
            int v = 0;
            while (sumeq == v && v < p_ck[i]) {
                for (int c = 0; c < k; c++) {
                    if (p_indx[i * nconds + v] == tempk[c] + 1) { 
                        sumeq += (p_implicants[i * nconds + p_indx[i * nconds + v] - 1] == tempc[c]);
                    }
                }
                v += 1;
            }
            if (sumeq == v) { 
                nonred = false; 
            }
            i += 1;
        }
    }
    return(nonred);
}
void push_PI(
    int p_implicants[],
    int p_indx[],
    int p_ck[],
    int p_pichart[],
    const int tempk[],
    const int tempc[],
    const int nconds,
    const int k,
    const int f,
    const int decpos[],
    const int frows[],
    const int posrows,
    const unsigned int foundPI
) {
    for (int c = 0; c < k; c++) {
        p_implicants[nconds * foundPI + tempk[c]] = tempc[c];
    }
    for (int c = 0; c < k; c++) {
        p_indx[nconds * foundPI + c] = tempk[c] + 1; 
    }
    p_ck[foundPI] = k;
    for (int r = 0; r < posrows; r++) {
        p_pichart[posrows * foundPI + r] = decpos[r] == decpos[frows[f]];
    }
}
void increment(
    int k,
    int *e,
    int *h,
    int nconds,
    int *tempk,
    int minval
) {
    if (k == 1) {
        tempk[0] += 1;
    }
    else {
        if (*e < nconds - *h) {
            *h = 1;
            tempk[k - 1] += 1;
            *e = tempk[k - 1];
            if (tempk[k - 1] < minval) {
                tempk[k - 1] = minval;
                *e = minval;
            }
        }
        else {
            *e = tempk[k - *h - 1] + 1;
            ++*h;
            Rboolean under = true;
            for (int j = 0; j < *h; j++) {
                under = under && (*e + j < minval);
                tempk[k - *h + j] = *e + j;
            }
            if (under) {
                *h = 1;
                tempk[k - *h] = minval;
                *e = minval;
            }
        }
    }
}
void populate_posneg(
    int *rowpos,
    int *rowneg,
    int nconds,
    int ttrows,
    int posrows,
    const int p_tt[],
    int posmat[],
    int negmat[]
) {
    int negrows = ttrows - posrows;
    for (int r = 0; r < ttrows; r++) {
        if (p_tt[nconds * ttrows + r] == 1) { 
            for (int c = 0; c < nconds; c++) {
                posmat[c * posrows + *rowpos] = p_tt[c * ttrows + r];
            }
            *rowpos += 1; 
        }
        else { 
            for (int c = 0; c < nconds; c++) {
                negmat[c * negrows + *rowneg] = p_tt[c * ttrows + r];
            }
            *rowneg += 1; 
        }
    }
    return;
}
void get_noflevels(
    int noflevels[],
    const int p_tt[],
    int nconds,
    int ttrows
) {
    for (int c = 0; c < nconds; c++) {
        noflevels[c] = 0; 
    }
    for (int c = 0; c < nconds; c++) {
        for (int r = 0; r < ttrows; r++) {
            if (noflevels[c] < p_tt[c * ttrows + r]) {
                noflevels[c] = p_tt[c * ttrows + r];
            }
        }
        noflevels[c] += 1; 
    }
    return;
}
void fill_mbase(
    int mbase[],
    const int tempk[],
    const int noflevels[],
    const int k
) {
    for (int c = 1; c < k; c++) {
        mbase[c] = mbase[c - 1] * noflevels[tempk[c - 1]];
    }
}
void get_decimals(
    const int posrows,
    const int negrows,
    const int k,
    int decpos[],
    int decneg[],
    const int posmat[],
    const int negmat[],
    const int tempk[],
    const int mbase[]
) {
    for (int r = 0; r < posrows; r++) {
        decpos[r] = 0;
        for (int c = 0; c < k; c++) {
            decpos[r] += posmat[tempk[c] * posrows + r] * mbase[c];
        }
    }
    for (int r = 0; r < negrows; r++) {
        decneg[r] = 0;
        for (int c = 0; c < k; c++) {
            decneg[r] += negmat[tempk[c] * negrows + r] * mbase[c];
        }
    }
}
void get_uniques(
    const int posrows,
    int *found,
    int decpos[],
    Rboolean possiblePI[],
    int possiblePIrows[]
) {
    for (int r = 1; r < posrows; r++) {
        int prev = 0;
        Rboolean unique = true; 
        while (prev < *found && unique) {
            unique = decpos[possiblePIrows[prev]] != decpos[r];
            prev += 1;
        }
        if (unique) {
            possiblePIrows[*found] = r;
            possiblePI[*found] = true;
            *found += 1;
        }
    }
}
void verify_possible_PI(
    const int compare,
    const int negrows,
    int *found,
    Rboolean possiblePI[],
    const int possiblePIrows[],
    const int decpos[],
    const int decneg[]
) {
    for (int i = 0; i < compare; i++) {
        int j = 0;
        while (j < negrows && possiblePI[i]) {
            if (decpos[possiblePIrows[i]] == decneg[j]) {
                possiblePI[i] = false;
                *found -= 1;
            }
            j += 1;
        }
    }
}
void get_frows(
    int frows[],
    const Rboolean possiblePI[],
    const int possiblePIrows[],
    const int compare
) {
    int pos = 0;
    for (int i = 0; i < compare; i++) {
        if (possiblePI[i]) {
            frows[pos] = possiblePIrows[i];
            pos += 1;
        }
    }
}
void fill_matrix(
    int nrows,
    int ncols,
    int nofl[],
    int *matrix,
    int startrow,
    int cols[],
    int plus1
) {
    int mbase[ncols];
    int orep[ncols];
    for (int c = 0; c < ncols; c++) {
        if (c == 0) {
            mbase[ncols - c - 1] = 1;
            orep[c] = 1;
        }
        else {
            mbase[ncols - c - 1] = mbase[ncols - c] * nofl[ncols - c];
            orep[c] = orep[c - 1] * nofl[c - 1];
        }
    }
    for (int c = 0; c < ncols; c++) {
        int lt = mbase[c] * nofl[c];
        for (int o = 0; o < orep[c]; o++) {
            for (int l = 0; l < nofl[c]; l++) {
                for (int i = 0; i < mbase[c]; i++) {
                    matrix[startrow + nrows * cols[c] + lt * o + mbase[c] * l + i] = l + plus1;
                }
            }
        }
    }
}
void calculate_rows(
    int *nrows,
    int ncols,
    int nofl[],
    int arrange,
    int maxprod
) {
    *nrows = 0;
    int e, h, k, prod;
    if (arrange == 0) {
        *nrows = 1;
        for (int c = 0; c < ncols; c++) {
            *nrows *= nofl[c]; 
        }
    }
    else {
        for (k = 1; k <= maxprod; k++) {
            int tempk[k];
            int nck = 1;
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
                for (int c = 0; c < k; c++) {
                    prod *= (nofl[tempk[c]] - 1);
                }
                *nrows += prod;
            }
        }
    }
}
int fill_combination_tasks(
    const int n, 
    int k, 
    int (*comb)[3], 
    const int arrayMaxSize, 
    int *numPrefixCombinations
) {
    *numPrefixCombinations = 0;
    if (k <= 1) {
        for (int x = 0; x <= n-k; x++)
        {
            comb[*numPrefixCombinations][0] = x;
            (*numPrefixCombinations)++;
        }
        return 1;
    }
    else if (k <= 2) {
        for (int x = 0; x <= n - k; x++)
            for (int y = x + 1; y <= n - k + 1; y++)
            {
                const int currIndex = *numPrefixCombinations;
                comb[currIndex][0] = x;
                comb[currIndex][1] = y;
                (*numPrefixCombinations)++;
            }
        return 2;
    }
    else { 
        for (int X = 0; X <= n-k; X++)
            for (int Y = X + 1; Y <= n-k+1; Y++)
                for (int Z = Y + 1; Z <=  n-k+2; Z++)
                {
                    const int currIndex = *numPrefixCombinations;
                    comb[currIndex][0] = X;
                    comb[currIndex][1] = Y;
                    comb[currIndex][2] = Z;
                    (*numPrefixCombinations)++;
                }
        return 3;
    }
}
int get_next_combination(
    int *arr,       
    int k,
    int n,          
    int indexLimit
) {
    int i;
    for (i = k - 1; i >= indexLimit; --i) {
        arr[i]++;
        if (arr[i] + (k - i - 1) != n)
            break;
    }
    if (i == indexLimit-1)
        return 0;
    for (int j = i + 1; j < k; j++)
        arr[j] = arr[j-1] + 1;
    return 1;
}
Rboolean all_covered(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols
) {
    Rboolean allrows = true;
    int r = 0;
    while (r < pirows && allrows) {
        Rboolean covered = false;
        unsigned int c = 0;
        while (c < picols && !covered) {
            covered = p_pichart[c * pirows + r];
            c++;
        }
        allrows = covered;
        r++;
    }
    return(allrows);
}
