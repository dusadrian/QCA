#include "qca_r.h"
#include <R_ext/RS.h>
#include "qca_r.h"
#include <math.h>
#include <stdlib.h>
#include "find_consistent_models.h"
#include "qca_threads.h"
#include "utils.h"

typedef struct {
    const int *p_implicants;
    const int *indx;
    const int *ck;
    const double *p_data;
    const int *p_fuzzy;
    int nconds;
    int nrdata;
    double solcons;
    double solcov;
    unsigned int foundPI;
    int k;
    int maxk;
    const int *p_sol;
    const int *cksol;
    unsigned int prevfound;
    unsigned char *valid;
} QCAFindConsistentModelsContext;

static void qca_unrank_combination(
    unsigned long long int task,
    unsigned int nitems,
    int k,
    int tempk[]
) {
    unsigned long long int combination = task;
    int x = 0;

    for (int i = 0; i < k; i++) {
        while (1) {
            unsigned long long int cval = nchoosek((int) nitems - (x + 1), k - (i + 1));
            if (cval == 0 || cval > combination) {
                break;
            }
            combination -= cval;
            x++;
        }

        if (x < 0) {
            x = 0;
        }
        if (x >= (int) nitems) {
            x = (int) nitems - 1;
        }

        tempk[i] = x;
        x++;
    }
}

static Rboolean qca_consistent_candidate_nonredundant(
    const QCAFindConsistentModelsContext *ctx,
    const int tempk[]
) {
    Rboolean nonred = true;
    unsigned int i = 0;

    while (i < ctx->prevfound && nonred) {
        int sumeq = 0;
        int v = 0;

        while (sumeq == v && v < ctx->cksol[i]) {
            for (int c = 0; c < ctx->k; c++) {
                if (ctx->p_sol[i * ctx->maxk + v] == tempk[c]) {
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

    return nonred;
}

static void qca_find_consistent_models_range(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    QCAFindConsistentModelsContext *ctx = (QCAFindConsistentModelsContext *) data;
    (void) worker_id;

    for (unsigned long long int task = start; task < end; task++) {
        int tempk[ctx->k];
        qca_unrank_combination(task, ctx->foundPI, ctx->k, tempk);

        if (
            qca_consistent_candidate_nonredundant(ctx, tempk) &&
            consistent_solution(
                ctx->p_data,
                ctx->nconds,
                ctx->nrdata,
                ctx->k,
                tempk,
                ctx->foundPI,
                ctx->p_implicants,
                ctx->ck,
                ctx->indx,
                ctx->p_fuzzy,
                ctx->solcons,
                ctx->solcov
            )
        ) {
            ctx->valid[task] = 1;
        }
    }
}



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
    
    // printf("soldepth: %d; maxk: %d\n", soldepth, maxk);
    
    
    // printf("maxk: %d\n", maxk);
    Rboolean keep_searching = true;
    int k = 1;
    double counter = 0; // counts combinations only when maxcomb > 0
    if (maxk == 0 || foundPI == 0) {
        // Nothing to search; return an empty matrix
        R_Free(*solutions);
        *solutions = R_Calloc(1, int);
        *nr = 0;
        *nc = 0;
        R_Free(p_sol);
        R_Free(cksol);
        return;
    }
    // while ((k <= maxk && counter > 0) || solfound == 0) {
    while (keep_searching && k <= maxk) {
        if (maxcomb <= 0) {
            unsigned long long int maxtasks = nchoosek((int) foundPI, k);
            unsigned char *valid = R_Calloc(maxtasks, unsigned char);
            unsigned int found_this_k = 0;

            QCAFindConsistentModelsContext ctx = {
                .p_implicants = p_implicants,
                .indx = indx,
                .ck = ck,
                .p_data = p_data,
                .p_fuzzy = p_fuzzy,
                .nconds = nconds,
                .nrdata = nrdata,
                .solcons = solcons,
                .solcov = solcov,
                .foundPI = foundPI,
                .k = k,
                .maxk = maxk,
                .p_sol = p_sol,
                .cksol = cksol,
                .prevfound = prevfound,
                .valid = valid
            };

            if (!qca_parallel_for(maxtasks, 0, qca_find_consistent_models_range, &ctx)) {
                R_Free(valid);
                R_Free(p_sol);
                R_Free(cksol);
                error("Failed to start pthread workers while finding consistent models.");
            }

            for (unsigned long long int task = 0; task < maxtasks; task++) {
                found_this_k += valid[task] != 0;
            }

            R_CheckUserInterrupt();

            while (solfound + found_this_k >= estimsol) {
                resize((void**)&p_sol, 1, 1000, estimsol, maxk);
                resize((void**)&cksol, 1, 1000, estimsol, 1);
                estimsol += 1000;
            }

            for (unsigned long long int task = 0; task < maxtasks; task++) {
                if (task > 0 && task % 1024 == 0) {
                    R_CheckUserInterrupt();
                }
                if (!valid[task]) {
                    continue;
                }

                int tempk[k];
                qca_unrank_combination(task, foundPI, k, tempk);

                for (int c = 0; c < k; c++) {
                    p_sol[solfound * maxk + c] = tempk[c];
                }
                cksol[solfound] = k;
                solfound++;
            }

            R_Free(valid);
        }
        else {
            int *tempk = R_Calloc(k, int);
            for (int i = 0; i < k; i++) {
                tempk[i] = i; // prepopulate with 0, 1, 2, ... , k - 1
            }
            
            tempk[k - 1] -= 1; // the last value is always increased with 1 by the function increment()
            
            int e = 0;
            int h = k;
            
            Rboolean last = (foundPI == k);
            
            // while ((tempk[0] != foundPI - k) || last) {
            while (keep_searching && ((tempk[0] != foundPI - k) || last)) {
                increment(k, &e, &h, foundPI + last, tempk, 0);
                last = false;
        

                Rboolean nonred = true;
                
                int i = 0;

                while (i < prevfound && nonred) {
        // printf("i: %d\n", i);
                // prevfound activates this loop only when allsol is activated
                int sumeq = 0;
                int v = 0;
                
                while (sumeq == v && v < cksol[i]) {
                    for (int c = 0; c < k; c++) {
        // printf("%d ", p_sol[i * maxk + v]);
                        if (p_sol[i * maxk + v] == tempk[c]) {
                            sumeq++;
                        }
                    }
                    v++;
                }
                
                if (sumeq == v) { // sumeq == ck[i], same thing
                    nonred = false; // make non-redundant false, the PI _is_ redundant
                }
                
                i++;
            }
            
        // printf("nonred: %d\n", nonred);
            if (nonred) {
                // test consistency and coverage


                if (consistent_solution(p_data, nconds, nrdata, k, tempk, foundPI, p_implicants, ck, indx, p_fuzzy, solcons, solcov)) {

        // printf("tempk: ");
        // for (int kk = 0; kk < k; kk++) {
        //     printf("%d ", tempk[kk]);
        // }
        // printf("\n");
        // printf("test: %d; incl: %5.3f; cov: %5.3f\n", consistent_solution(p_data, nconds, nrdata, k, tempk, foundPI, p_implicants, ck, indx, p_fuzzy, solcons, solcov), solcons, solcov);
                    for (int c = 0; c < k; c++) {
                        p_sol[solfound * maxk + c] = tempk[c];
                    }
                    
                    
                    cksol[solfound] = k;
                    
                    solfound++;
                    // counter = 2;
                    
                    // printf("solfound: %d", solfound);
                    if (solfound == estimsol) {
                        
                        resize((void**)&p_sol, 1, 1000, estimsol, maxk);
                        resize((void**)&cksol, 1, 1000, estimsol, 1);
                        
                        estimsol += 1000;
                        
                    }
                }
            }

            if (maxcomb > 0) {
                counter += 1.0;
                if ((counter / 1000000000.0) >= maxcomb) {
                    keep_searching = false;
                }
            }
            }
        
            R_Free(tempk);
        }

        prevfound = solfound;
        
        k += 1;
    }
    
    int *p_tempmat = R_Calloc(1, int);

    if (solfound > 0) {

        // Use the maximum k observed (not the last one) to size the matrix safely
        int finalrows = cksol[0];
        for (unsigned int i = 1; i < solfound; i++) {
            if (cksol[i] > finalrows) {
                finalrows = cksol[i];
            }
        }

        R_Free(p_tempmat);
        p_tempmat = R_Calloc(finalrows * solfound, int);
        
        for (int c = 0; c < solfound; c++) {
            for (int r = 0; r < cksol[c]; r++) {
                p_tempmat[c * finalrows + r] = p_sol[c * maxk + r] + 1; // + 1 to bring it in R notation
            }
        }

        *nr = finalrows;
        *nc = solfound;
        
    }
    else {
        // No consistent solutions found; return an empty matrix without dereferencing cksol.
        *nr = 0;
        *nc = 0;
    }

    R_Free(p_sol);
    R_Free(cksol);
    R_Free(*solutions);
    *solutions = p_tempmat;
}
