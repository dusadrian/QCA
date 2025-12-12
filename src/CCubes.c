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
#include <Rinternals.h>
#include <Rmath.h>
#include "CCubes.h"
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
    int prevfoundPI = 0;  
    int foundPI = 0;
    int prevsolmin = 0;     
    int solmin = 0;
    int previndices[posrows];
    int indices[posrows];
    for (int i = 0; i < posrows; i++) {
        previndices[i] = 0;
        indices[i] = 0;
    }
    #ifdef SHOW_DEBUG_PROFILE
        const double findingPIsStart_time = omp_get_wtime();
    #endif
    Rboolean solution_exists = false;
    int counter = 0; 
    int k;
    for (k = 1; k <= pidepth; k++) {
        Rboolean foundk = false;
        unsigned long long int maxtasks = nchoosek(nconds, k);
        if (maxtasks == 0) {
        }
        #ifdef _OPENMP
            #pragma omp parallel for schedule(static, 1)
        #endif
        for (unsigned long long int task = 0; task < maxtasks; task++) {
            int tempk[k];
            int x = 0;
            int combination = task;
            for (int i = 0; i < k; i++) {
                while (nchoosek(nconds - (x + 1), k - (i + 1)) <= combination) {
                    combination -= nchoosek(nconds - (x + 1), k - (i + 1));
                    x++;
                }
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
                        int word_index = tempk[c] / (BITS_PER_WORD / value_bit_width);
                        int bit_index = (tempk[c] % (BITS_PER_WORD / value_bit_width)) * value_bit_width;
                        fixed_bits[word_index] |= (((1 << value_bit_width) - 1) << bit_index);
                        value_bits[word_index] |= ((unsigned int)tempc[c] << bit_index);
                    }
                    Rboolean debug = false;
                    if (redundant(
                        p_implicants_pos,
                        p_implicants_val,
                        implicant_words,
                        fixed_bits,
                        value_bits,
                        prevfoundPI,
                        debug
                    )) {
                        continue;
                    }
                    Rboolean coverage[posrows];
                    unsigned int pichart_values[pichart_words];
                    for (int w = 0; w < pichart_words; w++) {
                        pichart_values[w] = 0U;
                    }
                    for (int r = 0; r < posrows; r++) {
                        coverage[r] = decpos[r] == decpos[frows[f]];
                        if (coverage[r]) {
                            int word_index = r / BITS_PER_WORD;
                            int bit_index = r % BITS_PER_WORD;
                            pichart_values[word_index] |= (1U << bit_index);
                        }
                    }
                    #ifdef _OPENMP
                        #pragma omp critical
                    #endif
                    {
                        for (int c = 0; c < k; c++) {
                            p_implicants[nconds * foundPI + tempk[c]] = tempc[c];
                        }
                        p_ck[foundPI] = k;
                        for (int c = 0; c < k; c++) {
                            p_indx[nconds * foundPI + c] = tempk[c] + 1;
                        }
                        for (int w = 0; w < implicant_words; w++) {
                            p_implicants_pos[implicant_words * foundPI + w] = fixed_bits[w];
                            p_implicants_val[implicant_words * foundPI] = value_bits[w];
                        }
                        for (int r = 0; r < posrows; r++) {
                            for (int w = 0; w < pichart_words; w++) {
                                p_pichart_pos[foundPI * pichart_words + w] = pichart_values[w];
                            }
                            p_pichart[posrows * foundPI + r] = coverage[r];
                        }
                        ++foundPI;
                        foundk = true;
                        if ((double)foundPI / (double)estimPI > 0.9) {
                            resize((void**)&p_pichart,        1, increase, estimPI, posrows);
                            resize((void**)&p_implicants,     1, increase, estimPI, nconds);
                            resize((void**)&p_indx,           1, increase, estimPI, nconds);
                            resize((void**)&p_implicants_val, 2, increase, estimPI, implicant_words);
                            resize((void**)&p_implicants_pos, 2, increase, estimPI, implicant_words);
                            resize((void**)&p_ck,             1, increase, estimPI, 1);
                            resize((void**)&p_pichart_pos,    2, increase, estimPI, pichart_words);
                            estimPI += increase;
                        }
                    }
                }
            }
        }
        if (foundPI > 0) {
            *complex = !gurobi && too_complex(foundPI, (solmin > 0 ? solmin : k), maxcomb);
            solution_exists = all_covered(p_pichart, posrows, foundPI);
            if (solution_exists) {
                stop_searching = *complex;
            }
            if (!stop_searching && solution_exists) { 
                SEXP pic = PROTECT(allocMatrix(INTSXP, posrows, foundPI));
                for (unsigned int i = 0; i < posrows * foundPI; i++) {
                    INTEGER(pic)[i] = p_pichart[i];
                }
                SEXP cpi = PROTECT(allocVector(INTSXP, 1));
                INTEGER(cpi)[0] = 0;
                setAttrib(pic, install("C_PI"), cpi);
                SEXP gurobipic = PROTECT(allocVector(LGLSXP, 1));
                LOGICAL(gurobipic)[0] = gurobi;
                setAttrib(pic, install("gurobi"), gurobipic);
                SEXP solindpic = PROTECT(allocVector(LGLSXP, 1));
                LOGICAL(solindpic)[0] = solind;
                setAttrib(pic, install("solind"), solindpic);
                R_ParseEvalString("library(QCA)", R_GlobalEnv);
                SEXP pkg_env = PROTECT(R_FindNamespace(mkString("QCA")));
                SEXP findmin = PROTECT(Rf_findVarInFrame(pkg_env, Rf_install("findmin")));
                SEXP evalinR = PROTECT(R_tryEval(Rf_lang2(findmin, pic), pkg_env, NULL));
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
                UNPROTECT(7);
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
        }
        if (stop_searching || counter > 1) {
            break;
        }
    }
    #ifdef SHOW_DEBUG_PROFILE
        const double findingPIsEnd_time = omp_get_wtime();
    #endif
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
                for (int r = 0; r < nconds; r++) {
                    copy_implicants[c * nconds + r] = p_implicants[p_sorted[c] * nconds + r];
                }
                for (int r = 0; r < posrows; r++) {
                    p_tempic[c * posrows + r] = p_pichart[p_sorted[c] * posrows + r];
                }
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
            for (int r = 0; r < nconds; r++) {
                copy_implicants[c * nconds + r] = p_implicants[p_sorted[c] * nconds + r];
            }
            for (int r = 0; r < posrows; r++) {
                p_tempic[c * posrows + r] = p_pichart[p_sorted[c] * posrows + r];
            }
        }
        if (solcons > 0) {
            int max_dim = (nconds > posrows) ? nconds : posrows;
            int *p_tempindx = R_Calloc(max_dim * foundPI, int);
            int *p_tempck = R_Calloc(foundPI, int);
            for (unsigned int c = 0; c < foundPI; c++) {
                for (int r = 0; r < nconds; r++) {
                    p_tempindx[c * nconds + r] = p_indx[p_sorted[c] * nconds + r];
                }
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
