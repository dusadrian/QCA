
#include "qca_r.h"
#include <R_ext/RS.h>
#include "qca_r.h"
#include <float.h> // for DBL_EPSILON
#include <math.h>
#include <limits.h>
#include "utils.h"


// #define BITS_PER_WORD 32
// #define BITS_PER_VALUE 4


// void printbits(unsigned int x, int bits_per_word, int value_bit_width)
// {
//     for (int i = bits_per_word - value_bit_width; i >= 0; i -= value_bit_width)
//     {
//         unsigned int chunk = (x >> i) & ((1U << value_bit_width) - 1);

//         // Print each chunk in binary
//         for (int j = value_bit_width - 1; j >= 0; j--)
//         {
//             printf("%d", (chunk >> j) & 1);
//         }
//         printf(" "); // Separate chunks for readability
//     }
//     printf("\n");
// }

int compute_value_bit_width(
    int max_value
) {
    if (max_value < 2) {
        return 1;  // Minimum width should be at least 1 bit
    }
    int bits_needed = ceil(log2(max_value + 1));  // Compute the necessary bits
    int power_of_2 = 1;
    while (power_of_2 < bits_needed) {
        power_of_2 *= 2;  // Round up to the next power of 2
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
    if (type == 1) { // type 1: int
        tmp = (int *)R_Calloc((size + increase) * nrows, int);
    } else if (type == 2) { // type 2: unsigned int
        tmp = (unsigned int *)R_Calloc((size + increase) * nrows, unsigned int);
    }

    if (tmp == NULL) {
        error("Memory allocation failed during resize.");
    }

    if (type == 1) { // type 1: int
        memcpy(tmp, *array, size * nrows * sizeof(int));
    } else if (type == 2) { // type 2: unsigned int
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

    // simulate n choose k
    for (int i = 1; i <= k; i++) {
        result *= n - (k - i);
        result /= i;
    }
    // printf("complexity: %5.3f\n", result / 1000000000);
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
        // for (int r = 0; r < nr; r++) {
        //     for (int c = 0; c < nc; c++) {
        //         tmp[r * nc + c] = matrix[c * nr + r];
        //     }
        // }
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
    const int covered[],
    unsigned int covered_limit,
    bool debug
) {
    Rboolean redundant = false;

    if (debug) {
        Rprintf("implicant_words: %d\n", implicant_words);
    }

    // check if the current possible PI is redundant
    unsigned int i = 0;
    unsigned int limit = (covered != NULL) ? covered_limit : prevfoundPI;
    while (i < limit && !redundant) {
        unsigned int pi_index = (covered != NULL) ? (unsigned int) covered[i] : i;
        // /*
        // a redundant PI is one for which all values from a previous PI are exactly the same:
        // 0 0 1 2 0, let's say previously found PI
        // then
        // 0 0 1 2 1 is redundant because on both columns 3 and 4 the values are equal
        // */

        Rboolean is_subset = true; // Assume it's a subset unless proven otherwise
        
        for (int w = 0; w < implicant_words; w++) {
            // if (debug && w == 0) {
            //     Rprintf("old_pos_%d: %u\n", i, p_implicants_pos[i * implicant_words + w]);
            //     printbits(p_implicants_pos[i * implicant_words]);
                
            //     Rprintf("old_val_%d: %u\n", i, p_implicants_val[i * implicant_words + w]);
            //     printbits(p_implicants_val[i * implicant_words]);
    
            //     Rprintf("new_pos: %u\n", fixed_bits[w]);
            //     printbits(fixed_bits[w]);
    
            //     Rprintf("new_val: %u\n", value_bits[w]);
            //     printbits(value_bits[w]);
                
            //     Rprintf("new_val & old_pos_%d: \n", i);
            //     printbits(value_bits[w] & p_implicants_pos[i * implicant_words + w]);
    
            //     Rprintf("old_val_%d & old_pos_%d: \n", i, i);
            //     printbits(p_implicants_val[i * implicant_words + w] & p_implicants_pos[i * implicant_words + w]);
            // }

            unsigned int pos_mask = p_implicants_pos[pi_index * implicant_words + w];
    
            if ((fixed_bits[w] & pos_mask) != pos_mask) {
                is_subset = false;
                break;
            }
    
            // Ensure P2 has the same values as P1 at those fixed positions
            if ((value_bits[w] & pos_mask) != (p_implicants_val[pi_index * implicant_words + w] & pos_mask)) {
                is_subset = false;
                break;
            }
        }

        redundant = is_subset;

        // if (debug) {
        //     Rprintf("redundant: %d\n", redundant);
        // }
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
        if (p_tt[nconds * ttrows + r] == 1) { // positive
            for (int c = 0; c < nconds; c++) {
                int value = p_tt[c * ttrows + r];
                posmat[c * posrows + *rowpos] = p_tt[c * ttrows + r];
                if (value > *max_value) {
                    *max_value = value;
                }
            }
            *rowpos += 1; // (*rowpos)++;
        }
        else { // negative
            for (int c = 0; c < nconds; c++) {
                int value = p_tt[c * ttrows + r];
                negmat[c * negrows + *rowneg] = p_tt[c * ttrows + r];
                if (value > *max_value) {
                    *max_value = value;
                }
            }
            *rowneg += 1; // (*rowneg)++;
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
        noflevels[c] = 0; // initiate with 0
        
        for (int r = 0; r < ttrows; r++) {
            if (noflevels[c] < p_tt[c * ttrows + r]) {
                noflevels[c] = p_tt[c * ttrows + r];
            }
        }

        noflevels[c] += 1; // add 1 because if the biggest number is 2 then it has three levels: 0, 1 and 2

        if (noflevels[c] == 1) {
            // no conditions ever has less than two levels
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
        Rboolean unique = true; // boolean flag, assume the row is unique
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

    // startrow is important, to append the matrix to another matrix
    // cols is there for the same reason, to indicate on which columns
    // from the previous matrix should this matrix be appended

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
        // number of rows is maximal
        *nrows = 1;
        // int cols[ncols];

        for (int c = 0; c < ncols; c++) {
            *nrows *= nofl[c]; // number of rows is maximal
            // cols[c] = c;
        }

    }
    else {

        // start counting the total number of rows
        // number of rows is not maximal
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
        // end counting the total number of rows
    }
}

Rboolean all_covered(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols
) {

    Rboolean allrows = true;

    for (int r = 0; r < pirows && allrows; r++) {

        Rboolean covered = false;

        for (unsigned int c = 0; c < picols && !covered; c++) {
            covered = p_pichart[c * pirows + r];
        }

        allrows = covered;
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
        // result = result * (n - i) / (i + 1);
        
        // Check for potential overflow before multiplication
        if (result > ULLONG_MAX / (n - i)) {
            return 0; // Indicate overflow
        }

        result *= (n - i);

        // Check for potential overflow before division
        if (result % (i + 1) != 0) {
            return 0; // Indicate overflow
        }
        
        result /= (i + 1);
    }

    return result;
}
