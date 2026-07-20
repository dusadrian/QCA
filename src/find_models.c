#include "qca_r.h"
#include <R_ext/RS.h> // for R_Calloc, R_free etc.
#include "qca_r.h"
#include <R_ext/Utils.h>
#include <stdlib.h>
#include "find_models.h"
#include "qca_threads.h"

#define INTERRUPT_EVERY 1024

typedef struct {
    const int *p_pichart;
    int pirows;
    unsigned int picols;
    int k;
    unsigned char *valid;
} QCAFindModelsValidContext;

static void qca_find_models_valid_range(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    QCAFindModelsValidContext *ctx = (QCAFindModelsValidContext *) data;
    (void) worker_id;

    for (unsigned long long int task = start; task < end; task++) {
        int tempk[ctx->k];
        unsigned long long int combination = task;
        int x = 0;

        for (int i = 0; i < ctx->k; i++) {
            while (1) {
                unsigned long long int cval = nchoosek(ctx->picols - (x + 1), ctx->k - (i + 1));
                if (cval == 0 || cval > combination) {
                    break;
                }
                combination -= cval;
                x++;
            }

            if (x < 0) {
                x = 0;
            }
            if (x >= (int) ctx->picols) {
                x = (int) ctx->picols - 1;
            }

            tempk[i] = x;
            x++;
        }

        Rboolean allrows = true;
        int r = 0;
        while (r < ctx->pirows && allrows) {
            Rboolean covered = false;
            int c = 0;
            while (c < ctx->k && !covered) {
                covered = ctx->p_pichart[tempk[c] * ctx->pirows + r];
                c++;
            }
            allrows = covered;
            r++;
        }

        if (allrows) {
            ctx->valid[task] = 1;
        }
    }
}

void find_models(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols,
    const Rboolean allsol,
    const int k,
    const double maxcomb,
    const Rboolean firstmin,
    int **solutions,
    int *nr,
    int *nc) {
    /*
    <pichart> is the (column-major) PI chart with the initial
    minterms on the rows and the PIs on the columns.

    This function assumes that solmin (here, k) is greater than zero
    and the caller already checked if a solution exists

    If the PI chart is certainly too complex (i.e. 34000 rows and 214 columns)
    lp_solve will strugle. Perhaps another idea would be to search combinations
    of k PIs, up to maxcomb combinations. If a certain combination pops up,
    the search should be stopped immediately and the solution reported.
    If no such minimal combination pops up until maxcomb, stop the search and report
    the problem as too complex to bother.
    */


    if (pirows <= 0 || picols == 0) {
        *nr = 0;
        *nc = 0;
        return;
    }

    // when k is equal to the number of columns, it means
    // all PIs are part of the solution and simply stop
    if (k == picols) {

        int *p_temp = (int *) R_Calloc(k, int);

        for (int i = 0; i < k; i++) {
            p_temp[i] = i + 1;
        }

        R_Free(*solutions);
        *solutions = p_temp;
        *nr = k;
        *nc = 1;

        return;
    }

    int *p_temp1 = (int *) R_Calloc(1, int);
    int *p_temp2 = (int *) R_Calloc(1, int);

    if (allsol) {

        // implements Petrick's method to derive all possible solutions

        int indmat[picols * pirows];
        int *mintpis = (int *) R_Calloc(pirows, int);

        for (int r = 0; r < pirows; r++) {
            if (r > 0 && r % INTERRUPT_EVERY == 0) {
                R_CheckUserInterrupt();
            }
            for (unsigned int c = 0; c < picols; c++) {
                if (p_pichart[c * pirows + r]) {
                    indmat[r * picols + mintpis[r]] = c;
                    mintpis[r]++;
                }
            }
        }

        // printfarray3(indmat, pirows, picols);
        R_Free(p_temp2);
        p_temp2 = (int *) R_Calloc(picols * mintpis[0], int);
        int *p_cols = (int *) R_Calloc(1, int); // just to initialize

        for (int i = 0; i < mintpis[0]; i++) {
            p_temp2[i * picols + indmat[i]] = 1;
        }

        unsigned int tempcols = mintpis[0];
        // printf("pirows: %d\n", pirows);
        for (int i = 1; i < pirows; i++) {
            if (i > 0 && i % INTERRUPT_EVERY == 0) {
                R_CheckUserInterrupt();
            }

            R_Free(p_temp1);
            p_temp1 = (int *) R_Calloc(picols * tempcols * mintpis[i], int);

            for (int j = 0; j < mintpis[i]; j++) {

                Memcpy(&p_temp1[j * tempcols * picols], p_temp2, tempcols * picols);

                for (unsigned int tc = 0; tc < tempcols; tc++) {
                    p_temp1[(j * tempcols + tc) * picols + indmat[i * picols + j]] = 1;
                }
            }
            // printf("temp2cols: %d\n", tempcols * mintpis[i]);
            unsigned int temp2cols = tempcols * mintpis[i];

            R_Free(p_cols);
            p_cols = (int *) R_Calloc(temp2cols, int);

            for (unsigned int i = 0; i < temp2cols; i++) {
                p_cols[i] = true;
            }

            unsigned int survcols = temp2cols;
            super_rows(p_temp1, picols, &survcols, p_cols);


            R_Free(p_temp2);
            p_temp2 = (int *) R_Calloc(picols * survcols, int);
            Memcpy(p_temp2, p_temp1, picols * survcols);

            tempcols = survcols;
            // printfarray3(p_temp2, picols, survcols);

        }

        R_Free(mintpis);

        R_Free(p_temp1);
        p_temp1 = (int *) R_Calloc(picols * tempcols, int);

        R_Free(p_cols);
        p_cols = (int *) R_Calloc(tempcols, int);

        int maxr = 0;

        for (int c = 0; c < tempcols; c++) {

            for (unsigned int r = 0; r < picols; r++) {
                if (p_temp2[c * picols + r]) {
                    p_temp1[c * picols + p_cols[c]] = r + 1;
                    p_cols[c]++;
                }

                if (maxr < p_cols[c]) {
                    maxr = p_cols[c];
                }
            }
        }

        R_Free(p_temp2);
        p_temp2 = (int *) R_Calloc(maxr * tempcols, int);

        for (unsigned int c = 0; c < tempcols; c++) {
            for (int r = 0; r < maxr; r++) {
                p_temp2[c * maxr + r] = p_temp1[c * picols + r];
            }
        }


        int temp;
        int order[tempcols];

        for (unsigned int c = 0; c < tempcols; c++) {
            order[c] = c;
        }


        // sort the entire matrix in ascending order by rows (PIs)
        for (int r = maxr - 1; r >= 0; r--) {
            for (unsigned int c1 = 0; c1 < tempcols; c1++) {
                for (unsigned int c2 = c1 + 1; c2 < tempcols; c2++) {
                    if (p_temp2[order[c1] * maxr + r] > p_temp2[order[c2] * maxr + r]) {

                        temp = order[c2];

                        for (int i = c2; i > c1; i--) {
                            order[i] = order[i - 1];
                        }

                        order[c1] = temp;
                    }
                }
            }
        }


        // then sort by the number of PIs in each solution, ascending order
        for (unsigned int c1 = 0; c1 < tempcols; c1++) {
            for (unsigned int c2 = c1 + 1; c2 < tempcols; c2++) {
                if (p_cols[order[c1]] > p_cols[order[c2]]) {

                    temp = order[c2];

                    for (int i = c2; i > c1; i--) {
                        order[i] = order[i - 1];
                    }

                    order[c1] = temp;
                }
            }
        }


        R_Free(p_cols);
        R_Free(p_temp1);
        p_temp1 = (int *) R_Calloc(maxr * tempcols, int);

        for (unsigned int c = 0; c < tempcols; c++) {
            for (int r = 0; r < maxr; r++) {
                p_temp1[c * maxr + r] = p_temp2[order[c] * maxr + r];
            }
        }

        *nr = maxr;
        *nc = tempcols;

    }
    else {

        unsigned int solfound = 0;

        if (!firstmin && maxcomb <= 0) {
            unsigned long long int maxtasks = nchoosek(picols, k);
            unsigned char *valid = (unsigned char *) R_Calloc(maxtasks, unsigned char);

            QCAFindModelsValidContext ctx = {
                .p_pichart = p_pichart,
                .pirows = pirows,
                .picols = picols,
                .k = k,
                .valid = valid
            };
            if (!qca_parallel_for(maxtasks, 0, qca_find_models_valid_range, &ctx)) {
                R_Free(valid);
                error("Failed to start pthread workers while finding models.");
            }
            for (unsigned long long int task = 0; task < maxtasks; task++) {
                solfound += valid[task] != 0;
            }

            R_CheckUserInterrupt();

            if (solfound > 0) {
                R_Free(p_temp1);
                p_temp1 = (int *) R_Calloc(k * solfound, int);

                unsigned int sol_index = 0;
                for (unsigned long long int task = 0; task < maxtasks; task++) {
                    if (task > 0 && task % INTERRUPT_EVERY == 0) {
                        R_CheckUserInterrupt();
                    }
                    if (!valid[task]) {
                        continue;
                    }

                    unsigned long long int combination = task;
                    int x = 0;

                    for (int i = 0; i < k; i++) {
                        while (1) {
                            unsigned long long int cval = nchoosek(picols - (x + 1), k - (i + 1));
                            if (cval == 0 || cval > combination) {
                                break;
                            }
                            combination -= cval;
                            x++;
                        }

                        if (x < 0) {
                            x = 0;
                        }
                        if (x >= (int) picols) {
                            x = (int) picols - 1;
                        }

                        p_temp1[sol_index * k + i] = x + 1;
                        x++;
                    }

                    sol_index++;
                }
            }
            else {
                R_Free(p_temp1);
                p_temp1 = (int *) R_Calloc(1, int);
            }

            R_Free(valid);
        }
        else {
            unsigned int estimsol = 100;

            R_Free(p_temp1);
            p_temp1 = (int *) R_Calloc(k * estimsol, int);
            Rboolean keep_searching = true;
            unsigned long long int maxtasks = nchoosek(picols, k);
            unsigned long long int counter = 0;

            for (unsigned long long int task = 0; keep_searching && task < maxtasks; task++) {
                if (task > 0 && task % INTERRUPT_EVERY == 0) {
                    R_CheckUserInterrupt();
                }

                int tempk[k];
                unsigned long long int combination = task;
                int x = 0;

                for (int i = 0; i < k; i++) {
                    while (1) {
                        unsigned long long int cval = nchoosek(picols - (x + 1), k - (i + 1));
                        if (cval == 0 || cval > combination) {
                            break;
                        }
                        combination -= cval;
                        x++;
                    }

                    if (x < 0) {
                        x = 0;
                    }
                    if (x >= (int) picols) {
                        x = (int) picols - 1;
                    }

                    tempk[i] = x;
                    x++;
                }

                // assume all initial minterms (on rows) are covered
                Rboolean allrows = true;

                int r = 0;
                while (r < pirows && allrows) {

                    // assume none of the PIs cover this initial minterm
                    Rboolean covered = false;

                    int c = 0;
                    while (c < k && !covered) {

                        // as soon as we find a PI which covers it, exit the loop
                        covered = p_pichart[tempk[c] * pirows + r];
                        c++;

                    }

                    // as soon as we find an initial minterm which is not covered by any of the PIs, exit the loop
                    allrows = covered;

                    r++;
                }


                if (allrows) {
                    for (int c = 0; c < k; c++) {
                        p_temp1[solfound * k + c] = tempk[c] + 1; // + 1 to bring it in R notation
                    }

                    solfound++;

                    if (solfound == estimsol) {
                        resize((void**)&p_temp1, 1, 100, estimsol, k);
                        estimsol += 100;
                    }
                }

                if (firstmin && solfound > 0) {
                    keep_searching = false;
                }

                if (maxcomb > 0) {
                    counter++;

                    if (((double) counter / 1000000000.0) >= maxcomb) {
                        keep_searching = false;
                    }
                }
            }            

            if (solfound > 0) {
                p_temp1 = (int *) R_Realloc(p_temp1, k * solfound, int);
            }
            else {
                R_Free(p_temp1);
                p_temp1 = (int *) R_Calloc(1, int);
            }
        }

        *nr = k;
        *nc = solfound;
    }

    R_Free(p_temp2);
    R_Free(*solutions);
    *solutions = p_temp1;

}
