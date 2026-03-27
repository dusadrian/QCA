/*
Copyright (c) 2016 - 2026, Adrian Dusa
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
#include <R_ext/Utils.h>
#include <Rinternals.h>
#include <Rmath.h>
#include <string.h>
#include "CCubes.h"
#include "solvePIchart_gurobi.h"
#include "solvePIchart_lagrangian.h"
#ifdef _OPENMP
    #undef match
    #include <omp.h>
    #define OMP_NUM_PROCS omp_get_num_procs()
    #define OMP_THREAD_LIMIT omp_get_thread_limit()
    #define OMP_MAX_THREADS omp_get_max_threads()
#else
    #define OMP_NUM_PROCS 1
    #define OMP_THREAD_LIMIT 1
    #define OMP_MAX_THREADS 1
#endif
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
void CCubes(const int p_tt[],          
            const int ttrows,          
            const int nconds,          
            const double p_data[],     
            const int nrdata,          
            const Rboolean allsol,     
            const Rboolean rowdom,     
            const double picons,       
            const int pidepth,         
            const int p_fsconds[],     
            const int soldepth,        
            const double solcons,
            const double solcov,
            const double maxcomb,
            const Rboolean keeptrying, 
            int **pichart,             
            int **implicants,          
            int **models,              
            unsigned int *foundPI_,    
            int *solrows,              
            int *solcols,              
            Rboolean *complex,         
            const Rboolean firstmin,   
            const Rboolean lagrangian, 
            const Rboolean gurobi,     
            const Rboolean solind      
) {
    int posrows = 0;
    for (int r = 0; r < ttrows; r++) {
        posrows += p_tt[nconds * ttrows + r];
    }
    int negrows = ttrows - posrows;
    int posmat[posrows * nconds];
    int negmat[nconds * negrows];
    int rowpos = 0, rowneg = 0;
    int max_value = 0;
    populate_posneg(&rowpos, &rowneg, nconds, ttrows, posrows, p_tt, posmat, negmat, &max_value);
    int value_bit_width = compute_value_bit_width(max_value + 1); 
    int noflevels[nconds];
    get_noflevels(noflevels, p_tt, nconds, ttrows);
    int estimPI = 25000;
    int increase = 25000; 
    int *p_pichart = (int *) R_Calloc(posrows * estimPI, int);
    int *p_indx = (int *) R_Calloc(nconds * estimPI, int);
    int *p_implicants = (int *) R_Calloc(nconds * estimPI, int);
    int pichart_words = (posrows + BITS_PER_WORD - 1) / BITS_PER_WORD; 
    unsigned int *p_pichart_pos = (unsigned int *) R_Calloc(estimPI * pichart_words, unsigned int);
    int implicant_words = (nconds + BITS_PER_WORD - 1) / BITS_PER_WORD; 
    unsigned int *p_implicants_pos = (unsigned int *) R_Calloc(estimPI * implicant_words, unsigned int);
    unsigned int *p_implicants_val = (unsigned int *) R_Calloc(estimPI * implicant_words, unsigned int);
    int *p_ck = (int *) R_Calloc(estimPI, int);
    Rboolean stop_searching = false;
    Rboolean solind_failed = false;
    Rboolean native_gurobi_checked = false;
    Rboolean native_gurobi_available = false;
    int prevfoundPI = 0;  
    int foundPI = 0;
    int prevsolmin = 0;     
    int solmin = 0;
    int previndices[posrows];
    int indices[posrows];
    int *covered = (int *) R_Calloc(estimPI, int);
    int *last_index = (int *) R_Calloc(posrows, int);
    int *k_last_index = (int *) R_Calloc(posrows, int);
    int nthreads = OMP_MAX_THREADS;
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
    #if defined(SHOW_DEBUG_PROFILE) && defined(_OPENMP)
        const double findingPIsStart_time = omp_get_wtime();
    #endif
    Rboolean solution_exists = false;
    int counter = 0; 
    int k;
    for (k = 1; k <= pidepth; k++) {
        R_CheckUserInterrupt();
        Rboolean foundk = false;
        for (int i = 0; i < posrows; i++) {
            k_last_index[i] = last_index[i];
        }
        unsigned long long int maxtasks = nchoosek(nconds, k);
        if (maxtasks == 0) {
        }
        #ifdef _OPENMP
            #pragma omp parallel for schedule(static, 1)
        #endif
        for (unsigned long long int task = 0; task < maxtasks; task++) {
            #ifndef _OPENMP
                if (task > 0 && task % INTERRUPT_EVERY == 0) {
                    R_CheckUserInterrupt();
                }
            #endif
            int tid = 0;
            #ifdef _OPENMP
                tid = omp_get_thread_num();
            #endif
            ThreadBuffer *tb = &buffers[tid];
            tb->found = 0;
            int tempk[k];
            int x = 0;
            unsigned long long int combination = task;
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
            int decpos[posrows];
            int decneg[negrows];
            int mbase[k];
            mbase[0] = 1; 
            #ifdef SHOW_DEBUG_OUTPUT
                #pragma omp critical
                {
                    int threadid = omp_get_thread_num();
                }
            #endif
            for (int c = 1; c < k; c++) {
                mbase[c] = mbase[c - 1] * noflevels[tempk[c - 1]];
            }
            get_decimals(posrows, negrows, k, decpos, decneg, posmat, negmat, tempk, mbase);
            int selected_word_index[k];
            int selected_bit_index[k];
            int bits_per_chunk = BITS_PER_WORD / value_bit_width;
            for (int c = 0; c < k; c++) {
                selected_word_index[c] = tempk[c] / bits_per_chunk;
                selected_bit_index[c] = (tempk[c] % bits_per_chunk) * value_bit_width;
            }
            int possiblePIrows[posrows];
            possiblePIrows[0] = 0; 
            Rboolean possiblePI[posrows];
            possiblePI[0] = true; 
            int found = 1;
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
                verify_possible_PI(compare, negrows, &found, possiblePI, possiblePIrows, decpos, decneg);
            }
            if (found > 0) {
                int frows[found];
                get_frows(frows, possiblePI, possiblePIrows, compare);
                for (int f = 0; f < found; f++) {
                    int tempc[k];
                    unsigned int fixed_bits[implicant_words];
                    unsigned int value_bits[implicant_words];
                    for (int i = 0; i < implicant_words; i++) {
                        fixed_bits[i] = 0U;
                        value_bits[i] = 0U;
                    }
                    for (int c = 0; c < k; c++) {
                        int value = posmat[tempk[c] * posrows + frows[f]];
                        tempc[c] = value + 1;
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
                #ifdef _OPENMP
                    #pragma omp critical
                #endif
                {
                    while ((double)(foundPI + tb->found) / (double)estimPI > 0.9) {
                        resize((void**)&p_pichart,        1, increase, estimPI, posrows);
                        resize((void**)&p_implicants,     1, increase, estimPI, nconds);
                        resize((void**)&p_indx,           1, increase, estimPI, nconds);
                        resize((void**)&p_implicants_val, 2, increase, estimPI, implicant_words);
                        resize((void**)&p_implicants_pos, 2, increase, estimPI, implicant_words);
                        resize((void**)&p_ck,             1, increase, estimPI, 1);
                        resize((void**)&p_pichart_pos,    2, increase, estimPI, pichart_words);
                        resize((void**)&covered,          1, increase, estimPI, 1);
                        estimPI += increase;
                    }
                    for (int bf = 0; bf < tb->found; bf++) {
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
        if (foundPI > 0) {
            *complex = !gurobi && too_complex(foundPI, (solmin > 0 ? solmin : k), maxcomb);
            solution_exists = all_covered(p_pichart, posrows, foundPI);
            if (solution_exists) {
                stop_searching = *complex;
            }
            if (!stop_searching && solution_exists) { 
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
                        &solmin
                    );
                    used_native_lagrangian = solmin > 0;
                    solind_failed = !used_native_lagrangian;
                    if (solmin < 0) {
                        solmin = 0;
                    }
                }
                if (!used_native_gurobi && !used_native_lagrangian) {
                    SEXP pic = PROTECT(allocMatrix(INTSXP, posrows, foundPI));
                    Memcpy(INTEGER(pic), p_pichart, posrows * foundPI);
                    SEXP cpi = PROTECT(allocVector(INTSXP, 1));
                    INTEGER(cpi)[0] = 0;
                    setAttrib(pic, install("C_PI"), cpi);
                    SEXP gurobipic = PROTECT(allocVector(LGLSXP, 1));
                    LOGICAL(gurobipic)[0] = gurobi;
                    setAttrib(pic, install("gurobi"), gurobipic);
                    R_ParseEvalString("library(QCA)", R_GlobalEnv);
                    SEXP pkg_env = PROTECT(R_FindNamespace(mkString("QCA")));
                    SEXP call = PROTECT(Rf_lang2(Rf_install("findmin"), pic));
                    SEXP evalinR = PROTECT(R_tryEval(call, pkg_env, NULL));
                    int leneval = length(evalinR);
                    solind_failed = leneval == 1 && TYPEOF(evalinR) == INTSXP;
                    if (solind_failed) {
                        solmin = INTEGER(evalinR)[0];
                    } else {
                        solmin = 0;
                        for (int i = 0; i < leneval; i++) {
                            if (REAL(evalinR)[i] > 0) {
                                indices[solmin] = i;
                                solmin++;
                            }
                        }
                    }
                    UNPROTECT(6);
                }
                if (solmin == prevsolmin) {
                    for (int i = 0; i < solmin; i++) {
                        indices[i] = previndices[i];
                    }
                    if (firstmin) {
                        counter += 1;
                    }
                }
                else {
                    prevsolmin = solmin;
                    for (int i = 0; i < solmin; i++) {
                        previndices[i] = indices[i];
                    }
                    counter = 0; 
                }
                if (!firstmin) {
                    if (foundk) {
                        counter = 0;
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
        if (stop_searching || counter > 1) {
            break;
        }
    }
    int *copy_implicants = R_Calloc(1, int);
    int *p_solutions = R_Calloc(1, int);
    int nr = 0, nc = 0; 
    int *p_tempic = R_Calloc(1, int);
    if ((firstmin || *complex) && solcons == 0) {
        if (solmin > 0) {
            R_Free(p_solutions);
            p_solutions = R_Calloc(solmin, int);
            for (int c = 0; c < solmin; c++) {
                p_solutions[c] = indices[c];
            }
            nr = solmin;
            nc = 1;
        }
        if (foundPI > 0) { 
            int *p_sorted = R_Calloc(foundPI, int);
            for (unsigned int i = 0; i < foundPI; i++) {
                p_sorted[i] = i;
            }
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
            if (keeptrying && solind_failed) {
                find_models(
                    p_tempic,
                    posrows,
                    foundPI,
                    false, 
                    k + 1,
                    maxcomb,
                    true, 
                    &p_solutions,
                    &nr,
                    &nc
                );
            }
        }
    }
    else if (foundPI > 0) { 
        if (rowdom) {
            row_dominance(p_pichart, p_implicants, p_ck, posrows, &foundPI, nconds);
        }
        int *p_sorted = R_Calloc(foundPI, int);
        for (unsigned int i = 0; i < foundPI; i++) {
            p_sorted[i] = i;
        }
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
            nr = solmin;
            find_models(
                p_tempic,
                posrows,
                foundPI,
                allsol,
                solmin,
                maxcomb,
                false, 
                &p_solutions,
                &nr,
                &nc
            );
        }
        R_Free(p_sorted);
    }
    R_Free(p_pichart);
    R_Free(p_implicants);
    R_Free(p_indx);
    R_Free(p_ck);
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
    R_Free(covered);
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
