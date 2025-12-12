/*
Copyright (c) 2016 - 2025, Adrian Dusa
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

#include <R.h>
#include <R_ext/RS.h>
#include <R_ext/Boolean.h>
#include <float.h> 
#include <math.h>
#include <limits.h>
#include "utils.h"
int compute_value_bit_width(
    int max_value
) {
    if (max_value < 2) {
        return 1;  
    }
    int bits_needed = ceil(log2(max_value + 1));  
    int power_of_2 = 1;
    while (power_of_2 < bits_needed) {
        power_of_2 *= 2;  
    }
    return power_of_2;
}
void resize(
    void **array,
    int type,
    int increase,
    int size,
    int nrows
) {
    if (type != 1 && type != 2) {
        error("Invalid type for resizing.");
    }
    void *tmp = NULL;
    if (type == 1) { 
        tmp = (int *)R_Calloc((size + increase) * nrows, int);
    } else if (type == 2) { 
        tmp = (unsigned int *)R_Calloc((size + increase) * nrows, unsigned int);
    }
    if (tmp == NULL) {
        error("Memory allocation failed during resize.");
    }
    if (type == 1) { 
        memcpy(tmp, *array, size * nrows * sizeof(int));
    } else if (type == 2) { 
        memcpy(tmp, *array, size * nrows * sizeof(unsigned int));
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
Rboolean redundant(
    unsigned int p_implicants_pos[],
    unsigned int p_implicants_val[],
    int implicant_words,
    unsigned int fixed_bits[],
    unsigned int value_bits[],
    unsigned int prevfoundPI,
    bool debug
) {
    Rboolean redundant = false;
    if (debug) {
        Rprintf("implicant_words: %d\n", implicant_words);
    }
    unsigned int i = 0;
    while (i < prevfoundPI && !redundant) {
        Rboolean is_subset = true; 
        for (int w = 0; w < implicant_words; w++) {
            unsigned int pos_mask = p_implicants_pos[i * implicant_words + w];
            if ((fixed_bits[w] & pos_mask) != pos_mask) {
                is_subset = false;
                break;
            }
            if ((value_bits[w] & pos_mask) != (p_implicants_val[i * implicant_words + w] & pos_mask)) {
                is_subset = false;
                break;
            }
        }
        redundant = is_subset;
        i++;
    }
    return(redundant);
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
    int negmat[],
    int *max_value
) {
    int negrows = ttrows - posrows;
    for (int r = 0; r < ttrows; r++) {
        if (p_tt[nconds * ttrows + r] == 1) { 
            for (int c = 0; c < nconds; c++) {
                int value = p_tt[c * ttrows + r];
                posmat[c * posrows + *rowpos] = p_tt[c * ttrows + r];
                if (value > *max_value) {
                    *max_value = value;
                }
            }
            *rowpos += 1; 
        }
        else { 
            for (int c = 0; c < nconds; c++) {
                int value = p_tt[c * ttrows + r];
                negmat[c * negrows + *rowneg] = p_tt[c * ttrows + r];
                if (value > *max_value) {
                    *max_value = value;
                }
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
        for (int r = 0; r < ttrows; r++) {
            if (noflevels[c] < p_tt[c * ttrows + r]) {
                noflevels[c] = p_tt[c * ttrows + r];
            }
        }
        noflevels[c] += 1; 
        if (noflevels[c] == 1) {
            noflevels[c] = 2;
        }
    }
    return;
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
            prev++;
        }
        if (unique) {
            possiblePIrows[*found] = r;
            possiblePI[*found] = true;
            (*found)++;
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
                (*found)--;
            }
            j++;
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
            pos++;
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
unsigned long long int nchoosek(int n, int k) {
    if (k > n) return 0;
    if (k == 0 || k == n) return 1;
    unsigned long long int result = 1;
    if (k > n - k) {
        k = n - k;
    }
    for (int i = 0; i < k; i++) {
        if (result > ULLONG_MAX / (n - i)) {
            return 0; 
        }
        result *= (n - i);
        if (result % (i + 1) != 0) {
            return 0; 
        }
        result /= (i + 1);
    }
    return result;
}
