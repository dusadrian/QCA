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

#include <R.h>
#include <R_ext/RS.h> 
#include <R_ext/Boolean.h>
#include <Rinternals.h>
#include <math.h>
#include "consistency.h"
#include "find_models.h"
#include "find_consistent_models.h"
#include "CCubes.h"
void CCubes(const int p_tt[],       
            const int ttrows,       
            const int nconds,       
            const double p_data[],  
            const int nrdata,       
            const Rboolean allsol,      
            const Rboolean rowdom,      
            const Rboolean minpin,      
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
            const Rboolean firstmin     
) {
    int *p_pichart, *p_implicants, *p_indx, *p_ck;
    int posrows = 0;
    for (int r = 0; r < ttrows; r++) {
        posrows += p_tt[nconds * ttrows + r];
    }
    int negrows = ttrows - posrows;
    int posmat[posrows * nconds];
    int negmat[nconds * ((negrows > 0) ? negrows : 1)]; 
    int rowpos = 0, rowneg = 0;
    populate_posneg(&rowpos, &rowneg, nconds, ttrows, posrows, p_tt, posmat, negmat);
    int noflevels[nconds];
    for (int i = 0; i < nconds; i++) {
        noflevels[i] = 2; 
    }
    get_noflevels(noflevels, p_tt, nconds, ttrows);
    unsigned int estimPI = 10000;
    p_pichart = R_Calloc(posrows * estimPI, int);
    p_implicants = R_Calloc(nconds * estimPI, int);
    p_indx = R_Calloc(nconds * estimPI, int);
    p_ck = R_Calloc(estimPI, int);
    Rboolean stop_searching = false;
    unsigned int prevfoundPI = 0;  
    unsigned int foundPI = 0;
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
        #ifdef SHOW_DEBUG_OUTPUT
        #endif 
        const int numMaxCombinations = (nconds*(nconds-1)*(nconds-2))/6 + (nconds*(nconds-1)/2) + nconds;
        int combinationsSetFixedPrefix[numMaxCombinations][3];
        int numPrefixCombinations = 0;
	    
        const int INDEX_LIMIT = fillCombinationTasks(nconds, k, combinationsSetFixedPrefix, numMaxCombinations, &numPrefixCombinations);
        for (int taskIter = 0; taskIter < numPrefixCombinations; taskIter++) {
            int tempk[k];
            for (int i = 0; i < INDEX_LIMIT; i++) {
                tempk[i] = combinationsSetFixedPrefix[taskIter][i];
            }
            for (int i = INDEX_LIMIT; i < k; i++) {
                tempk[i] = tempk[i - 1] + 1;
            }
            if (k > INDEX_LIMIT)
                tempk[k - 1]--; 
            int decpos[posrows];
            int decneg[(negrows > 0) ? negrows : 1]; 
            int finishedAll = 0;
            int mbase[k];
            mbase[0] = 1; 
            while (!finishedAll) {
                if (k > INDEX_LIMIT) {
                    int res = GetNextComb(tempk, k, nconds, INDEX_LIMIT);
                    if (!res) {
                        finishedAll = 1;
                        break; 
                    }                    
                }
                else {
                    finishedAll = 1;
                }
                #ifdef SHOW_DEBUG_OUTPUT
                    #pragma omp critical
                    {
                        int threadid = omp_get_thread_num();
                    }
                #endif
                fill_mbase(mbase, tempk, noflevels, k);
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
                            found -= 1;
                        }
                    }
                }
                else if (negrows > 0) {
                    verify_possible_PI(compare, negrows, &found, possiblePI, possiblePIrows, decpos, decneg);
                }
                if (found) { 
                    int frows[found];
                    get_frows(frows, possiblePI, possiblePIrows, compare);
                    for (int f = 0; f < found; f++) {
                        int tempc[k];
                        for (int c = 0; c < k; c++) {
                            tempc[c] = posmat[tempk[c] * posrows + frows[f]] + 1;
                        }
                        if (nonredundant(p_implicants, p_indx, p_ck, tempk, tempc, nconds, k, prevfoundPI)) {
                            {
                                push_PI(
                                    p_implicants,
                                    p_indx,
                                    p_ck,
                                    p_pichart,
                                    tempk,
                                    tempc,
                                    nconds,
                                    k,
                                    f,
                                    decpos,
                                    frows,
                                    posrows,
                                    foundPI
                                );
                                ++foundPI;
                                foundk = true;
                                if (foundPI / estimPI > 0.9) {
                                    estimPI += 10000;
                                    resize(&p_pichart,    posrows, estimPI, foundPI);
                                    resize(&p_implicants, nconds,  estimPI, foundPI);
                                    resize(&p_indx,       nconds,  estimPI, foundPI);
                                    resize(&p_ck,         1,       estimPI, foundPI);
                                }
                            }
                        }
                    }
                }
            }
        }
        if (foundPI > 0) {
            *complex = too_complex(foundPI, (solmin > 0 ? solmin : k), maxcomb);
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
                Rf_eval(lang2(Rf_install("library"), mkString("QCA")), R_GlobalEnv);
                SEXP pkg_env = R_FindNamespace(mkString("QCA"));
                SEXP findmin = findVarInFrame(pkg_env, Rf_install("findmin"));
                solmin = INTEGER(R_tryEval(lang2(findmin, pic), pkg_env, NULL))[0];
                UNPROTECT(2);
                if (solmin == prevsolmin) {
                    if (firstmin || minpin) {
                        for (int i = 0; i < solmin; i++) {
                            indices[i] = previndices[i];
                        }
                    }
                }
                else {
                    prevsolmin = solmin;
                    for (int i = 0; i < solmin; i++) {
                        previndices[i] = indices[i];
                    }
                }
                if (foundk) {
                    counter = 0;
                }
                else {
                    counter += 1;
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
            if (keeptrying) {
                find_models(p_tempic, posrows, foundPI, false, k + 1, maxcomb, true, &p_solutions, &nr, &nc);
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
            int *p_tempindx = R_Calloc(posrows * foundPI, int);
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
            find_models(p_tempic, posrows, foundPI, allsol, solmin, maxcomb, false, &p_solutions, &nr, &nc);
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
