#include "qca_rinternals.h"
#include <R_ext/RS.h>
#include <stdbool.h>
#include <stdlib.h>
#include <limits.h>
#include <stdint.h>
#include <string.h>
#include "find_consistent_models.h"
#include "qca_threads.h"
#include "utils.h"

typedef struct {
    const QCAConsistentFit *fit;
    const QCAConsistentBounds *bounds;
    unsigned int foundPI;
    int k;
    int maxk;
    const int *p_sol;
    const int *cksol;
    unsigned int prevfound;
} QCAFindConsistentModelsContext;

static Rboolean qca_consistent_candidate_nonredundant(
    const QCAFindConsistentModelsContext *ctx,
    const int tempk[],
    int terms
) {
    Rboolean nonred = true;
    unsigned int i = 0;

    while (i < ctx->prevfound && nonred) {
        int sumeq = 0;
        int v = 0;

        while (sumeq == v && v < ctx->cksol[i]) {
            for (int c = 0; c < terms; c++) {
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

/* Each prefix occupies a contiguous interval in lexicographic rank order.
   Remember passes within the worker's range and jump over rejected prefixes. */
static unsigned long long qca_skip_failed_prefix(
    const QCAFindConsistentModelsContext *ctx,
    unsigned long long rank,
    const int selected[],
    const unsigned long long prefix_end[],
    unsigned long long checked_until[]
) {
    for (int d = 0; d < ctx->k - 1; ++d) {
        if (rank < checked_until[d]) continue;
        checked_until[d] = prefix_end[d];
        if (prefix_end[d] > rank + 1 &&
            !qca_consistent_prefix_possible(ctx->fit, ctx->bounds, selected, d + 1, ctx->k - d - 1))
            return prefix_end[d];
    }
    return rank;
}

/* A worker resumes a lexicographic prefix search after each bounded batch.
   Storage grows with the models found, never with choose(foundPI, k). */
typedef struct {
    int *selected, *models;
    int depth, root_start, root_end;
    size_t count, capacity;
    bool done, failed, larger_possible;
} QCAConsistentTree;

typedef struct {
    QCAFindConsistentModelsContext candidates;
    QCAConsistentTree *trees;
} QCAConsistentLevel;

static bool qca_larger_consistent_model_possible(
    const QCAFindConsistentModelsContext *ctx, const int selected[], int depth
) {
    int available = (int)ctx->foundPI - selected[depth - 1] - 1;
    if (ctx->maxk <= ctx->k || available < ctx->k + 1 - depth) return false;
    if (!ctx->bounds->suffix_memberships) return true;
    int slots = ctx->maxk - depth;
    if (slots > available) slots = available;
    return qca_consistent_prefix_possible(
        ctx->fit, ctx->bounds, selected, depth, slots);
}

static void qca_record_consistent_model(QCAConsistentTree *tree, int terms)
{
    size_t limit = (size_t)INT_MAX / terms;
    if (tree->count == limit) {
        tree->failed = true;
        return;
    }
    if (tree->count == tree->capacity) {
        size_t capacity = tree->capacity ? tree->capacity * 2 : 64;
        if (capacity > limit) capacity = limit;
        if (capacity > SIZE_MAX / sizeof(int) / terms) {
            tree->failed = true;
            return;
        }
        int *models = realloc(tree->models, capacity * terms * sizeof(int));
        if (!models) {
            tree->failed = true;
            return;
        }
        tree->models = models;
        tree->capacity = capacity;
    }
    memcpy(tree->models + tree->count * terms, tree->selected, terms * sizeof(int));
    ++tree->count;
}

static void qca_search_consistent_tree(QCAConsistentLevel *level, int tree_id)
{
    const QCAFindConsistentModelsContext *ctx = &level->candidates;
    QCAConsistentTree *tree = &level->trees[tree_id];
    if (tree->done || tree->failed) return;
    if (!tree->selected) {
        tree->selected = malloc((size_t)ctx->k * sizeof(int));
        if (!tree->selected) {
            tree->failed = true;
            return;
        }
        tree->selected[0] = tree->root_start - 1;
    }

    for (int work = 0; work < 4096 && !tree->done && !tree->failed; ++work) {
        int depth = tree->depth;
        int last = (int)ctx->foundPI - (ctx->k - depth);
        if (depth == 0 && last >= tree->root_end) last = tree->root_end - 1;
        if (++tree->selected[depth] > last) {
            if (depth == 0) tree->done = true;
            else --tree->depth;
            continue;
        }

        int terms = depth + 1;
        if (!qca_consistent_candidate_nonredundant(ctx, tree->selected, terms))
            continue;

        if (terms == ctx->k) {
            if (consistent_solution(ctx->fit, terms, tree->selected)) {
                qca_record_consistent_model(tree, terms);
            } else if (!tree->larger_possible) {
                tree->larger_possible = qca_larger_consistent_model_possible(
                    ctx, tree->selected, terms);
            }
            continue;
        }

        if (ctx->bounds->suffix_memberships &&
            !qca_consistent_prefix_possible(
                ctx->fit, ctx->bounds, tree->selected, terms, ctx->k - terms)) {
            /* Failure at this size does not exclude larger models. Check the
               whole remaining depth before declaring this branch exhausted. */
            if (!tree->larger_possible)
                tree->larger_possible = qca_larger_consistent_model_possible(
                    ctx, tree->selected, terms);
            continue;
        }
        tree->selected[terms] = tree->selected[depth];
        ++tree->depth;
    }
}

static void qca_search_consistent_prefixes(
    unsigned long long start, unsigned long long end, int worker_id, void *data
) {
    (void)worker_id;
    for (unsigned long long tree = start; tree < end; ++tree)
        qca_search_consistent_tree(data, (int)tree);
}

typedef struct {
    QCAConsistentFit fit;
    QCAConsistentBounds bounds;
    unsigned int foundPI;
    int maxk;
    double maxcomb;
    int **solutions;
    int *nr, *nc;
    int *models, *sizes;
    QCAConsistentTree *trees;
    int workers;
} QCAConsistentSearch;

static void qca_reserve_consistent_models(
    QCAConsistentSearch *search, size_t needed, unsigned int *capacity
) {
    size_t limit = (size_t)INT_MAX / search->maxk;
    if (needed > limit) error("Too many consistent models for the solution matrix.");
    if (needed <= *capacity) return;
    size_t grown = (size_t)*capacity * 2;
    if (grown < needed) grown = needed;
    if (grown > limit) grown = limit;
    search->models = R_Realloc(search->models, grown * search->maxk, int);
    search->sizes = R_Realloc(search->sizes, grown, int);
    *capacity = grown;
}

static void free_consistent_search(void *data)
{
    QCAConsistentSearch *search = data;
    if (search->trees) {
        for (int i = 0; i < search->workers; ++i) {
            free(search->trees[i].selected);
            free(search->trees[i].models);
        }
        free(search->trees);
    }
    R_Free(search->models);
    R_Free(search->sizes);
    free(search->fit.memberships);
    free(search->bounds.suffix_memberships);
    free(search->bounds.suffix_overlap);
}

static SEXP search_consistent_models(void *data)
{
    QCAConsistentSearch *search = data;
    const unsigned int foundPI = search->foundPI;
    const int maxk = search->maxk;
    const double maxcomb = search->maxcomb;
    int **solutions = search->solutions;
    int *nr = search->nr, *nc = search->nc;
    if (maxk == 0 || foundPI == 0) {
        R_Free(*solutions);
        *solutions = R_Calloc(1, int);
        *nr = 0;
        *nc = 0;
        return R_NilValue;
    }

    qca_cache_pi_memberships(&search->fit, foundPI);
    if (!search->fit.memberships)
        error("Cannot allocate memberships for consistent-model search.");
    qca_prepare_consistent_bounds(&search->fit, foundPI, maxk, &search->bounds);
    unsigned int estimsol = (unsigned int)INT_MAX / maxk;
    if (estimsol > 1000) estimsol = 1000;

    int *p_sol = search->models = R_Calloc((size_t)maxk * estimsol, int);
    int *cksol = search->sizes = R_Calloc(estimsol, int);
    
    unsigned int solfound = 0;
    unsigned int prevfound = 0;
    
    // printf("soldepth: %d; maxk: %d\n", soldepth, maxk);
    
    
    // printf("maxk: %d\n", maxk);
    Rboolean keep_searching = true;
    int k = 1;
    double counter = 0; // counts combinations only when maxcomb > 0
    // while ((k <= maxk && counter > 0) || solfound == 0) {
    while (keep_searching && k <= maxk) {
        QCAFindConsistentModelsContext ctx = {
            .fit = &search->fit, .bounds = &search->bounds,
            .foundPI = foundPI, .k = k, .maxk = maxk,
            .p_sol = p_sol, .cksol = cksol, .prevfound = prevfound
        };
        const bool possible = !search->bounds.suffix_memberships ||
            qca_consistent_prefix_possible(&search->fit, &search->bounds, NULL, 0, k);
        bool larger_possible = true;
        if (maxcomb <= 0) {
            unsigned int roots = foundPI - k + 1;
            search->workers = qca_default_thread_count();
            if ((unsigned int)search->workers > roots) search->workers = roots;
            search->trees = calloc((size_t)search->workers, sizeof(QCAConsistentTree));
            if (!search->trees) error("Cannot allocate consistent-model search workers.");
            unsigned int start = 0;
            for (int w = 0; w < search->workers; ++w) {
                unsigned int width = roots / search->workers +
                    ((unsigned int)w < roots % search->workers);
                search->trees[w].root_start = start;
                search->trees[w].root_end = start + width;
                start += width;
            }
            QCAConsistentLevel level = {.candidates = ctx, .trees = search->trees};

            bool finished = false;
            while (!finished) {
                if (!qca_parallel_for(search->workers, search->workers, qca_search_consistent_prefixes, &level))
                    error("Failed to start consistent-model search workers.");
                finished = true;
                for (int w = 0; w < search->workers; ++w) {
                    if (search->trees[w].failed)
                        error("Consistent-model output exceeds available memory or matrix limits.");
                    if (!search->trees[w].done) finished = false;
                }
                // Workers have joined, so interrupts cannot leave live threads behind.
                R_CheckUserInterrupt();
            }

            larger_possible = false;
            for (int w = 0; w < search->workers; ++w) {
                QCAConsistentTree *tree = &search->trees[w];
                if (tree->larger_possible) larger_possible = true;
                qca_reserve_consistent_models(search, solfound + tree->count, &estimsol);
                p_sol = search->models;
                cksol = search->sizes;
                for (size_t m = 0; m < tree->count; ++m) {
                    memcpy(p_sol + (size_t)solfound * maxk,
                           tree->models + m * k, (size_t)k * sizeof(int));
                    cksol[solfound++] = k;
                }
                free(tree->selected);
                free(tree->models);
                tree->selected = tree->models = NULL;
            }
            free(search->trees);
            search->trees = NULL;
        }
        else {
            int *tempk = R_Calloc(k, int);
            for (int i = 0; i < k; ++i) tempk[i] = i;
            --tempk[k - 1]; // increment() advances to the first combination.

            int e = 0, h = k;
            Rboolean last = foundPI == k;
            const bool bound = possible && search->bounds.suffix_memberships && k > 1 &&
                nchoosek((int)foundPI, k) > 1;
            unsigned long long rank = 0, next_possible = 0;
            const int prefix_capacity = bound ? k : 1;
            unsigned long long prefix_end[prefix_capacity], checked_until[prefix_capacity];
            if (bound) Memzero(checked_until, k);
            while (keep_searching && ((tempk[0] != foundPI - k) || last)) {
                increment(k, &e, &h, foundPI + last, tempk, 0);
                last = false;

                if (bound && rank >= next_possible && rank >= checked_until[k - 2]) {
                    /* Count completions after this rank for each prefix. Use
                       increment() and charge every rank to max.comb as before. */
                    prefix_end[k - 1] = rank + 1;
                    for (int d = k - 1; d > 0; --d)
                        prefix_end[d - 1] = prefix_end[d] +
                            nchoosek((int)foundPI - tempk[d] - 1, k - d);
                    next_possible = qca_skip_failed_prefix(
                        &ctx, rank, tempk, prefix_end, checked_until);
                }
                if (possible && rank >= next_possible &&
                    qca_consistent_candidate_nonredundant(&ctx, tempk, k) &&
                    consistent_solution(&search->fit, k, tempk)) {
                    qca_reserve_consistent_models(search, (size_t)solfound + 1, &estimsol);
                    p_sol = search->models;
                    cksol = search->sizes;
                    ctx.p_sol = p_sol;
                    ctx.cksol = cksol;
                    for (int c = 0; c < k; ++c) p_sol[(size_t)solfound * maxk + c] = tempk[c];
                    cksol[solfound++] = k;
                }
                ++rank;
                if (maxcomb > 0) {
                    counter += 1.0;
                    if ((counter / 1000000000.0) >= maxcomb) keep_searching = false;
                }
            }
            R_Free(tempk);
        }

        // Every larger model has a prefix at this size. Stop only when all
        // such branches are redundant or cannot meet the fit thresholds.
        if (!larger_possible) break;
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

    R_Free(*solutions);
    *solutions = p_tempmat;
    return R_NilValue;
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
    // Any membership union has a representation using at most one maximizer
    // per case. A larger model has a proper subset with exactly the same fit.
    // Positive truth-table rows are not a valid bound for fuzzy memberships.
    int maxk = nrdata;
    if (foundPI < maxk) maxk = foundPI;
    if (soldepth > 0 && soldepth < maxk) maxk = soldepth;

    QCAConsistentSearch search = {
        .fit = {
            .data = p_data, .implicants = p_implicants, .depths = ck,
            .indices = indx, .fuzzy = p_fuzzy, .conditions = nconds, .cases = nrdata,
            .consistency = solcons, .coverage = solcov
        },
        .foundPI = foundPI, .maxk = maxk, .maxcomb = maxcomb,
        .solutions = solutions, .nr = nr, .nc = nc
    };
    R_ExecWithCleanup(search_consistent_models, &search, free_consistent_search, &search);
}
