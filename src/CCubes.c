
#include "qca_rinternals.h"
#include <R_ext/RS.h> // for R_Calloc, R_free etc.
#include "qca_rinternals.h"
#include <R_ext/Utils.h>
#include "qca_rinternals.h"
#include <Rmath.h>
#include <stdlib.h>
#include <string.h>
#include "CCubes.h"
#include "solvePIchart_gurobi.h"
#include "solvePIchart_lagrangian.h"
#include "findmin_lpsolve.h"

static Rboolean resize_worker_buffer(
    void **array,
    int type,
    int increase,
    int size,
    int nrows
) {
    size_t elem_size = type == 1 ? sizeof(int) : sizeof(unsigned int);
    size_t old_count = (size_t) size * (size_t) nrows;
    size_t new_count = (size_t) (size + increase) * (size_t) nrows;
    void *tmp = calloc(new_count, elem_size);

    if (tmp == NULL) {
        return FALSE;
    }

    memcpy(tmp, *array, old_count * elem_size);
    free(*array);
    *array = tmp;
    return TRUE;
}

// void printfarray(int* arr, int rows, int cols) {
//     for (int c = 0; c < cols; c++) {
//         for (int r = 0; r < rows; r++) {
//             printf("%d ", arr[rows * c + r]);
//         }
//         printf("\n");
//     }
// }


#define BITS_PER_WORD 32
#define INTERRUPT_EVERY 1024

typedef struct {
    int *implicants;
    int *indx;
    int *ck;
    int *covsum;
    unsigned int *implicants_pos;
    unsigned int *implicants_val;
    unsigned int *pichart_pos;
    int found;
} ThreadBuffer;

void CCubes(const int p_tt[],          // the truth table
            const int ttrows,          // number of rows in the truth table
            const int nconds,          // number of conditions in the data
            const double p_data[],     // calibrated dataset
            const int nrdata,          // number of rows in the raw data
            const Rboolean allsol,     // all solution
            const Rboolean rowdom,     // row dominance
            const double picons,       // PI consistency
            const int pidepth,         // depth in number of conditions for each PI
            const int p_fsconds[],     // are conditions fuzzy?
            const int soldepth,        // depth in number of PIs for each solution
            const double solcons,
            const double solcov,
            const double maxcomb,
            const Rboolean keeptrying, // try hard to find at least one solution if everything else fails

            // pointers to save the results
            int **pichart,             // final PI chart
            int **implicants,          // final implicants in matrix form
            int **models,              // final solution models
            unsigned int *foundPI_,    // final number of PIs found
            int *solrows,              // number of rows for the solution matrix
            int *solcols,              // number of columns for the solution matrix
            Rboolean *complex,         // signal if the returned solution is incomplete due to a too complex PI chart

            const Rboolean firstmin,   // IEEE switch
            const Rboolean lagrangian, // whether using the lagrangian backend via findmin()
            const Rboolean gurobi,     // whether to use Gurobi (if installed), default is TRUE
            const Rboolean solind      // IF using Gurobi, findmin() returns the indexes of the solutions
) {

    // calculate the number of positive rows
    int posrows = 0;
    for (int r = 0; r < ttrows; r++) {
        posrows += p_tt[nconds * ttrows + r];
    }

    // calculate the number of negative rows
    int negrows = ttrows - posrows;

    // split the list data into two matrices:
    // one for the positive and the other for the negative rows
    int posmat[posrows * nconds];
    int negmat[nconds * negrows];
    int rowpos = 0, rowneg = 0;
    int max_value = 0;

    populate_posneg(&rowpos, &rowneg, nconds, ttrows, posrows, p_tt, posmat, negmat, &max_value);

    int value_bit_width = compute_value_bit_width(max_value + 1); // + 1 because shifting PI values uses 0 as a don't care

    // just to make sure we _do_ initialize negmat (... why... ?)
    // calculate the number of levels for each causal condition, finding the biggest number in each column
    int noflevels[nconds];

    get_noflevels(noflevels, p_tt, nconds, ttrows);

    // preallocating for an estimated large number of 1000 found PIs
    // this number will be iteratively increased when the found PIs reach the upper limit
    int estimPI = 25000;
    int increase = 25000; // how much to increase the size of the arrays when needed

    // p_pichart = malloc(posrows * estimPI * sizeof(int));
    // memset(p_pichart, false, posrows * estimPI * sizeof(int));
    int *p_pichart = (int *) calloc((size_t) posrows * (size_t) estimPI, sizeof(int));
    // calloc() prefills all values with 0s

    int *p_indx = (int *) calloc((size_t) nconds * (size_t) estimPI, sizeof(int));
    int *p_implicants = (int *) calloc((size_t) nconds * (size_t) estimPI, sizeof(int));

    int pichart_words = (posrows + BITS_PER_WORD - 1) / BITS_PER_WORD; // Words needed per PI chart columns
    unsigned int *p_pichart_pos = (unsigned int *) calloc((size_t) estimPI * (size_t) pichart_words, sizeof(unsigned int));
    int implicant_words = (nconds + BITS_PER_WORD - 1) / BITS_PER_WORD; // Words needed per PI representation
    unsigned int *p_implicants_pos = (unsigned int *) calloc((size_t) estimPI * (size_t) implicant_words, sizeof(unsigned int));
    unsigned int *p_implicants_val = (unsigned int *) calloc((size_t) estimPI * (size_t) implicant_words, sizeof(unsigned int));

    // vector of ints containing the complexity level for each found, non-redundant PI
    int *p_ck = (int *) calloc((size_t) estimPI, sizeof(int));
    if (
        p_pichart == NULL || p_indx == NULL || p_implicants == NULL ||
        p_pichart_pos == NULL || p_implicants_pos == NULL ||
        p_implicants_val == NULL || p_ck == NULL
    ) {
        error("Memory allocation failed during PI buffer initialization.");
    }

    Rboolean stop_searching = false;
    Rboolean solind_failed = false;
    Rboolean native_gurobi_checked = false;
    Rboolean native_gurobi_available = false;

    // prev: in the PREVIOUS level of complexity
    int prevfoundPI = 0;  // the number of previously found PIs
    int foundPI = 0;
    int prevsolmin = 0;     // the minimum number of PIs that solve the PI chart
    int solmin = 0;

    // the positions of the PIs solving the PI chart
    // a vector which can never be lengthier than the number of minterms (posrows)
    int previndices[posrows];
    int indices[posrows];
    int *covered = (int *) calloc((size_t) estimPI, sizeof(int));
    if (covered == NULL) {
        error("Memory allocation failed during PI coverage initialization.");
    }
    int *last_index = (int *) R_Calloc(posrows, int);
    int *k_last_index = (int *) R_Calloc(posrows, int);

    int nthreads = 1;
    ThreadBuffer *buffers = (ThreadBuffer *) R_Calloc(nthreads, ThreadBuffer);

    for (int t = 0; t < nthreads; t++) {
        buffers[t].implicants = (int *) R_Calloc(posrows * nconds, int);
        buffers[t].indx = (int *) R_Calloc(posrows * nconds, int);
        buffers[t].ck = (int *) R_Calloc(posrows, int);
        buffers[t].covsum = (int *) R_Calloc(posrows, int);
        buffers[t].implicants_pos = (unsigned int *) R_Calloc(posrows * implicant_words, unsigned int);
        buffers[t].implicants_val = (unsigned int *) R_Calloc(posrows * implicant_words, unsigned int);
        buffers[t].pichart_pos = (unsigned int *) R_Calloc(posrows * pichart_words, unsigned int);
        buffers[t].found = 0;
    }


    for (int i = 0; i < posrows; i++) {
        previndices[i] = 0;
        indices[i] = 0;
    }


    Rboolean solution_exists = false;
    int counter = 0; // to stop if two consecutive levels of complexity yield no more PIs
    Rboolean pi_resize_failed = false;
    int k;
    for (k = 1; k <= pidepth; k++) {
        R_CheckUserInterrupt();
        // Rprintf("====================\nk: %d\n\n", k);

        Rboolean foundk = false;
        for (int i = 0; i < posrows; i++) {
            k_last_index[i] = last_index[i];
        }

        unsigned long long int maxtasks = nchoosek(nconds, k);
        if (maxtasks == 0) {
            // overflow, too many tasks
            // return(R_NilValue);
        }

        for (unsigned long long int task = 0; task < maxtasks; task++) {
            if (task > 0 && task % INTERRUPT_EVERY == 0) {
                R_CheckUserInterrupt();
            }
            int tid = 0;

            ThreadBuffer *tb = &buffers[tid];
            tb->found = 0;
            // this will contain each instance of n choose k
            int tempk[k];
            int x = 0;
            unsigned long long int combination = task;

            // fill the combination for the current task / combination number
            for (int i = 0; i < k; i++) {
                while (1) {
                    unsigned long long int cval = nchoosek(nconds - (x + 1), k - (i + 1));
                    if (cval == 0 || cval > combination) {
                        break;
                    }
                    combination -= cval;
                    x++;
                }
                if (x < 0) x = 0;
                if (x >= nconds) x = nconds - 1;
                tempk[i] = x;
                x++;
            }

            // allocate vectors of decimal row numbers for the positive and negative rows
            int decpos[posrows];
            int decneg[negrows];

            // create the vector of multiple bases, useful when calculating the decimal representation
            // of a particular combination of columns, for each row
            int mbase[k];
            mbase[0] = 1; // the first number is _always_ equal to 1, irrespective of the number of levels in a certain column

            // tempk now contains the combinations of columns from the initial data, for a given complexity level k

            // calculate the vector of multiple bases, for example if we have k = 3 (three columns) with 2, 3 and 2 levels
            // then mbase will be [1, 2, 6] 1, 1 * 2 = 2, 2 * 3 = 6
            for (int c = 1; c < k; c++) {
                mbase[c] = mbase[c - 1] * noflevels[tempk[c - 1]];
            }

            // calculate decimal numbers, using mbase
            get_decimals(posrows, negrows, k, decpos, decneg, posmat, negmat, tempk, mbase);

            int selected_word_index[k];
            int selected_bit_index[k];
            int bits_per_chunk = BITS_PER_WORD / value_bit_width;
            for (int c = 0; c < k; c++) {
                selected_word_index[c] = tempk[c] / bits_per_chunk;
                selected_bit_index[c] = (tempk[c] % bits_per_chunk) * value_bit_width;
            }


            int possiblePIrows[posrows];
            possiblePIrows[0] = 0; // first row is always a possible PI, to compare with the negative rows

            Rboolean possiblePI[posrows];
            possiblePI[0] = true; // boolean flag, to be set with 0 if found among the negative rows

            int found = 1;

            // identifies all unique decimal rows
            get_uniques(posrows, &found, decpos, possiblePI, possiblePIrows);

            int compare = found;

            if (picons > 0) {

                int val[k];
                int fuzzy[k];

                for (int i = 0; i < compare; i++) {

                    for (int c = 0; c < k; c++) {
                        val[c] = posmat[tempk[c] * posrows + possiblePIrows[i]];
                        fuzzy[c] = p_fsconds[tempk[c]] * 1;
                    }

                    if (altb(consistency(p_data, nrdata, nconds, k, tempk, val, fuzzy), picons)) {
                        possiblePI[i] = false;
                        found--;
                    }
                }
            }
            else if (negrows > 0) {
                // check if the same decimal number is found in the negative rows
                verify_possible_PI(compare, negrows, &found, possiblePI, possiblePIrows, decpos, decneg);
            }


            if (found > 0) {

                // some of the positive row numbers are possible PIs
                // (not found in the negative rows, or having a consistency greater than pi.cons)
                int frows[found];
                get_frows(frows, possiblePI, possiblePIrows, compare);

                for (int f = 0; f < found; f++) {

                    // Rprintf("k: %d; task: %llu; f: %d\n", k, task, f);

                    // create a temporary vector of length k, containing the values from the initial positive
                    // matrix plus 1 (because 0 now signals a minimization, it becomes 1, and 1 becomes 2 etc.
                    int tempc[k];

                    // using bit shifting, store the fixed bits and value bits
                    unsigned int fixed_bits[implicant_words];
                    unsigned int value_bits[implicant_words];
                    
                    for (int i = 0; i < implicant_words; i++) {
                        fixed_bits[i] = 0U;
                        value_bits[i] = 0U;
                    }

                    for (int c = 0; c < k; c++) {
                        int value = posmat[tempk[c] * posrows + frows[f]];
                        tempc[c] = value + 1;
                        // (+ 1 because 0 now signals a don't care)

                        int word_index = selected_word_index[c];
                        int bit_index = selected_bit_index[c];

                        fixed_bits[word_index] |= (((1 << value_bit_width) - 1) << bit_index);
                        value_bits[word_index] |= ((unsigned int)tempc[c] << bit_index);
                    }

                    int covsum = 0;
                    int covered_rows[posrows];
                    for (int r = 0; r < posrows; r++) {
                        if (decpos[r] == decpos[frows[f]]) {
                            covered_rows[covsum] = r;
                            covsum++;
                        }
                    }

                    // Rboolean debug = k == 2 && task == 3 && f == 3;
                    Rboolean debug = false;

                    unsigned int redundant_limit = (covsum > 0) ? (unsigned int) last_index[covsum - 1] : (unsigned int) prevfoundPI;

                    if (redundant(
                        p_implicants_pos,
                        p_implicants_val,
                        implicant_words,
                        fixed_bits,
                        value_bits,
                        prevfoundPI,
                        covered,
                        redundant_limit,
                        debug
                    )) {
                        continue;
                    }

                    unsigned int pichart_values[pichart_words];
                    Memzero(pichart_values, pichart_words);

                    for (int r = 0; r < covsum; r++) {
                        int row = covered_rows[r];
                        int word_index = row / BITS_PER_WORD;
                        int bit_index = row % BITS_PER_WORD;
                        pichart_values[word_index] |= (1U << bit_index);
                    }

                    // Rprintf("It is a prime implicant\n");
                    // This operation first gets a new index to push in the global array in a concurent way
                    // Then adds the result there.
                    // WE CAN Synchronize only the index and let the copy operation happen in parallel BUT this
                    // creates a false sharing problem and the performance is down by several factors.

                    int bf = tb->found;
                    int base_imp = bf * nconds;
                    int base_bits = bf * implicant_words;
                    int base_pic = bf * pichart_words;

                    Memzero(&tb->implicants[base_imp], nconds);
                    Memzero(&tb->indx[base_imp], nconds);

                    for (int c = 0; c < k; c++) {
                        tb->implicants[base_imp + tempk[c]] = tempc[c];
                        tb->indx[base_imp + c] = tempk[c] + 1;
                    }

                    tb->ck[bf] = k;
                    tb->covsum[bf] = covsum;

                    for (int w = 0; w < implicant_words; w++) {
                        tb->implicants_pos[base_bits + w] = fixed_bits[w];
                        tb->implicants_val[base_bits + w] = value_bits[w];
                    }

                    for (int w = 0; w < pichart_words; w++) {
                        tb->pichart_pos[base_pic + w] = pichart_values[w];
                    }

                    tb->found++;
                }
            }

            if (tb->found > 0) {
                {
                    if (!pi_resize_failed) {
                        while ((double)(foundPI + tb->found) / (double)estimPI > 0.9) {
                            if (
                                !resize_worker_buffer((void**)&p_pichart,        1, increase, estimPI, posrows) ||
                                !resize_worker_buffer((void**)&p_implicants,     1, increase, estimPI, nconds) ||
                                !resize_worker_buffer((void**)&p_indx,           1, increase, estimPI, nconds) ||
                                !resize_worker_buffer((void**)&p_implicants_val, 2, increase, estimPI, implicant_words) ||
                                !resize_worker_buffer((void**)&p_implicants_pos, 2, increase, estimPI, implicant_words) ||
                                !resize_worker_buffer((void**)&p_ck,             1, increase, estimPI, 1) ||
                                !resize_worker_buffer((void**)&p_pichart_pos,    2, increase, estimPI, pichart_words) ||
                                !resize_worker_buffer((void**)&covered,          1, increase, estimPI, 1)
                            ) {
                                pi_resize_failed = true;
                                break;
                            }
                            estimPI += increase;
                        }
                    }

                    for (int bf = 0; !pi_resize_failed && bf < tb->found; bf++) {
                        int base_imp = bf * nconds;
                        int base_bits = bf * implicant_words;
                        int base_pic = bf * pichart_words;
                        int covsum = tb->covsum[bf];
                        if (covsum < 1) covsum = 1;
                        if (covsum > posrows) covsum = posrows;

                        int samelevel_start = last_index[covsum - 1];
                        int samelevel_end = k_last_index[covsum - 1];
                        int samelevel_count = samelevel_end - samelevel_start;
                        if (samelevel_count > 0) {
                            Rboolean samelevel_redundant = redundant(
                                p_implicants_pos,
                                p_implicants_val,
                                implicant_words,
                                &tb->implicants_pos[base_bits],
                                &tb->implicants_val[base_bits],
                                foundPI,
                                covered + samelevel_start,
                                samelevel_count,
                                false
                            );
                            if (samelevel_redundant) {
                                continue;
                            }
                        }

                        Memcpy(&p_implicants[nconds * foundPI], &tb->implicants[base_imp], nconds);
                        Memcpy(&p_indx[nconds * foundPI], &tb->indx[base_imp], nconds);

                        p_ck[foundPI] = tb->ck[bf];

                        Memcpy(&p_implicants_pos[implicant_words * foundPI], &tb->implicants_pos[base_bits], implicant_words);
                        Memcpy(&p_implicants_val[implicant_words * foundPI], &tb->implicants_val[base_bits], implicant_words);

                        Memcpy(&p_pichart_pos[foundPI * pichart_words], &tb->pichart_pos[base_pic], pichart_words);

                        for (int r = 0; r < posrows; r++) {
                            int word_index = r / BITS_PER_WORD;
                            int bit_index = r % BITS_PER_WORD;
                            p_pichart[posrows * foundPI + r] =
                                (tb->pichart_pos[base_pic + word_index] & (1U << bit_index)) > 0;
                        }

                        int insert_at = k_last_index[covsum - 1];
                        if (foundPI > insert_at) {
                            memmove(&covered[insert_at + 1], &covered[insert_at], (size_t) (foundPI - insert_at) * sizeof(int));
                        }
                        covered[insert_at] = foundPI;

                        for (int i = 0; i < covsum; i++) {
                            k_last_index[i]++;
                        }

                        ++foundPI;
                    }

                    foundk = true;
                }
            }
        }
        R_CheckUserInterrupt();
        if (pi_resize_failed) {
            error("Memory allocation failed during PI buffer resize.");
        }
        if (foundPI > 0) {

            *complex = !gurobi && too_complex(foundPI, (solmin > 0 ? solmin : k), maxcomb);

            solution_exists = all_covered(p_pichart, posrows, foundPI);

            if (solution_exists) {
                stop_searching = *complex;
            }

            // printf("posrows: %d; foundPI: %d; complex: %d\n", posrows, foundPI, *complex);
            // printf("PI chart can be solved: %d\n", solution_exists);

            if (!stop_searching && solution_exists) { // not too complex

                // SEXP pic = PROTECT(allocMatrix(INTSXP, foundPI, posrows));

                // int i, j, l_1, len;
                // len = foundPI * posrows;
                // l_1 = len - 1;
                // // transpose to a row-major with PIs on the rows
                // for (i = 0, j = 0; i < len; i++, j += posrows) {
                //     if (j > l_1) j -= l_1;
                //     INTEGER(pic)[i] = p_pichart[j];
                // }

                Rboolean used_native_gurobi = false;
                Rboolean used_native_lagrangian = false;

                if (!lagrangian && gurobi) {
                    if (!native_gurobi_checked) {
                        native_gurobi_available = gurobi_runtime_available();
                        native_gurobi_checked = true;
                    }

                    if (native_gurobi_available) {
                        used_native_gurobi = solvePIchart_gurobi(
                            p_pichart,
                            foundPI,
                            posrows,
                            indices,
                            &solmin
                        );

                        if (used_native_gurobi) {
                            solind_failed = false;
                        }
                        else {
                            native_gurobi_available = false;
                        }
                    }
                }

                if (lagrangian) {
                    solvePIchart_lagrangian(
                        p_pichart,
                        foundPI,
                        posrows,
                        NULL,
                        indices,
                        &solmin,
                        NULL,
                        NULL
                    );
                    used_native_lagrangian = solmin > 0;
                    solind_failed = !used_native_lagrangian;
                    if (solmin < 0) {
                        solmin = 0;
                    }
                }

                if (!used_native_gurobi && !used_native_lagrangian) {
                    solmin = 0;
                    solind_failed = !solvePIchart_lpsolve(
                        p_pichart,
                        posrows,
                        foundPI,
                        indices,
                        &solmin
                    );
                }
                // Rprintf("solution minima: %d\n", solmin);

                // find_min(p_pichart, posrows, foundPI, &solmin, indices);

                if (solmin == prevsolmin) {
                    // the minimum did not change in the current level of complexity

                    for (int i = 0; i < solmin; i++) {
                        indices[i] = previndices[i];
                    }

                    if (firstmin) {
                        counter += 1;
                    }
                }
                else {
                    // this means solmin is in fact smaller than the previously found solmin
                    // or it is the very first time a solmin was found
                    // only here it makes sense to overwrite prevsolmin and previndices,
                    // otherwise they are just as good as the ones from the previous level

                    prevsolmin = solmin;
                    for (int i = 0; i < solmin; i++) {
                        previndices[i] = indices[i];
                    }

                    counter = 0; // this means for sure foundk is true
                }

                if (!firstmin) {
                    if (foundk) {
                        counter = 0;
                        // prevminPI = foundPI;
                    }
                    else {
                        counter += 1;
                    }
                }
            }

            prevfoundPI = foundPI;
            for (int i = 0; i < posrows; i++) {
                last_index[i] = k_last_index[i];
            }
        }


        // printf("counter: %d; stop: %d\n", counter, stop_searching * 1);
        
        if (stop_searching || counter > 1) {
            break;
        }
    }

    int *copy_implicants = R_Calloc(1, int);

    int *p_solutions = R_Calloc(1, int);
    int nr = 0, nc = 0; // for the solution matrix

    int *p_tempic = R_Calloc(1, int);

    if ((firstmin || *complex) && solcons == 0) {

        if (solmin > 0) {
            // copy (only) the PIs from the solution
            R_Free(p_solutions);
            p_solutions = R_Calloc(solmin, int);

            for (int c = 0; c < solmin; c++) {
                p_solutions[c] = indices[c];
            }

            nr = solmin;
            nc = 1;
        }

        if (foundPI > 0) { // make sure foundPI > 0 to allocate memory > 0...!!

            // too complex problem to find a solution for
            // just copy the (sorted) PIs and the PI chart

            int *p_sorted = R_Calloc(foundPI, int);
            // int sorted[foundPI];
            for (unsigned int i = 0; i < foundPI; i++) {
                p_sorted[i] = i;
            }

            // sort PIs with the least complex (negative) first
            sort_cols(p_implicants, p_sorted, p_ck, nconds, foundPI);

            R_Free(copy_implicants);
            copy_implicants = R_Calloc(nconds * foundPI, int);

            R_Free(p_tempic);
            p_tempic = R_Calloc(posrows * foundPI, int);


            for (int r = 0; r < solmin; r++) {
                for (int c = 0; c < foundPI; c++) {
                    if (p_sorted[c] == p_solutions[r]) {
                        indices[r] = c;
                    }
                }
            }

            for (int r = 0; r < solmin; r++) {
                p_solutions[r] = indices[r];
            }


            for (unsigned int c = 0; c < foundPI; c++) {
                Memcpy(&copy_implicants[c * nconds], &p_implicants[p_sorted[c] * nconds], nconds);
                Memcpy(&p_tempic[c * posrows], &p_pichart[p_sorted[c] * posrows], posrows);
            }

            R_Free(p_sorted);
            // printf("keeptrying: %d\n", keeptrying * 1);
            if (keeptrying && solind_failed) {
                find_models(
                    p_tempic,
                    posrows,
                    foundPI,
                    false, // allsol
                    k + 1,
                    maxcomb,
                    true, // firstmin
                    &p_solutions,
                    &nr,
                    &nc
                );
                // printf("nr: %d; nc: %d\n", nr, nc);
            }
        }
    }
    else if (foundPI > 0) { // make sure foundPI > 0 to allocate memory > 0...!!

        if (rowdom) {
            // "PI dominance", which changes foundPI unless there is no dominated PI
            // this automatically adjusts p_pichart, p_implicants and p_ck
            row_dominance(p_pichart, p_implicants, p_ck, posrows, &foundPI, nconds);
        }

        int *p_sorted = R_Calloc(foundPI, int);
        // int sorted[foundPI];
        for (unsigned int i = 0; i < foundPI; i++) {
            p_sorted[i] = i;
        }

        // sort PIs with the least complex (negative) first
        sort_cols(p_implicants, p_sorted, p_ck, nconds, foundPI);

        R_Free(copy_implicants);
        copy_implicants = R_Calloc(nconds * foundPI, int);

        R_Free(p_tempic);
        p_tempic = R_Calloc(posrows * foundPI, int);


        for (unsigned int c = 0; c < foundPI; c++) {
            Memcpy(&copy_implicants[c * nconds], &p_implicants[p_sorted[c] * nconds], nconds);
            Memcpy(&p_tempic[c * posrows], &p_pichart[p_sorted[c] * posrows], posrows);
        }

        if (solcons > 0) {
            // Allocate using the larger dimension to avoid overflow if future code
            // ever iterates by posrows instead of nconds.
            int max_dim = (nconds > posrows) ? nconds : posrows;
            int *p_tempindx = R_Calloc(max_dim * foundPI, int);
            int *p_tempck = R_Calloc(foundPI, int);

            for (unsigned int c = 0; c < foundPI; c++) {
                Memcpy(&p_tempindx[c * nconds], &p_indx[p_sorted[c] * nconds], nconds);
                p_tempck[c] = p_ck[p_sorted[c]];
            }

            find_consistent_models(
                copy_implicants,
                p_tempindx,
                p_tempck,
                p_data,
                p_fsconds,
                nconds,
                nrdata,
                posrows,
                solcons,
                solcov,
                allsol,
                soldepth,
                foundPI,
                maxcomb,
                &p_solutions,
                &nr,
                &nc
            );

            R_Free(p_tempindx);
            R_Free(p_tempck);

        }
        else if (solmin > 0) {
            // the solutions matrix has no less than solmin number of rows, but if allsol it can have more
            nr = solmin;

            find_models(
                p_tempic,
                posrows,
                foundPI,
                allsol,
                solmin,
                maxcomb,
                false, // firstmin
                &p_solutions,
                &nr,
                &nc
            );
        }

        R_Free(p_sorted);
    }

    free(p_pichart);
    free(p_implicants);
    free(p_indx);
    free(p_ck);
    free(p_pichart_pos);
    free(p_implicants_pos);
    free(p_implicants_val);
    for (int t = 0; t < nthreads; t++) {
        R_Free(buffers[t].implicants);
        R_Free(buffers[t].indx);
        R_Free(buffers[t].ck);
        R_Free(buffers[t].covsum);
        R_Free(buffers[t].implicants_pos);
        R_Free(buffers[t].implicants_val);
        R_Free(buffers[t].pichart_pos);
    }
    R_Free(buffers);
    free(covered);
    R_Free(last_index);
    R_Free(k_last_index);
    R_Free(*models);
    R_Free(*implicants);
    R_Free(*pichart);

    *solrows = nr;
    *solcols = nc;
    *foundPI_ = foundPI;
    *implicants = copy_implicants;
    *pichart = p_tempic;
    *models = p_solutions;

    return;
}
