#include "scp_solver.h"
#include "scp_relaxation.h"

#include <limits.h>
#include <math.h>
#include <R_ext/RS.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <time.h>
#include <string.h>

static void *qca_scp_calloc(size_t count, size_t size) {
    return R_chk_calloc((R_SIZE_T) count, (R_SIZE_T) size);
}

static void qca_scp_free(void *ptr) {
    if (ptr != NULL) {
        R_chk_free(ptr);
    }
}

typedef struct {
    qca_scp_state state;
    int *branch_cols;
    int *gains;
    int *candidate_solution;
    int64_t *dual;
} qca_scp_workspace;

typedef struct {
    const qca_scp_problem *problem;
    int *best_solution;
    int best_size;
    int proof_target_size;
    int found_target;
    int limit_reached;
    unsigned long long nodes;
    unsigned long long node_limit;
    unsigned long long hard_node_limit;
    double started;
    double time_limit_seconds;
    double hard_time_limit_seconds;
    int adaptive_limits;
    int checkpoint_best_size;
    unsigned long long checkpoint_root_branches;
    unsigned long long root_branches_total;
    unsigned long long root_branches_completed;
    qca_scp_workspace *workspaces;
    unsigned int *reduction_row_marks;
    int *reduction_row_support;
    int *reduction_gains;
    int *reduction_active_cols;
    double *reduction_dual_slack;
    int *branching_row_heuristic;
    int64_t *lagrangian_costs;
    int64_t *best_dual;
    int *gradient;
} qca_scp_search;

/* The incumbent cutoff tightens in the same tree, without target restarts. */
static int propagation_target(const qca_scp_search *search) {
    int target = search->proof_target_size;
    if (search->problem->target_propagation &&
        (target < 0 || search->best_size - 1 < target)) {
        target = search->best_size - 1;
    }
    return target;
}

static qca_scp_profile qca_profile = {0};

static double now_seconds(void) {
    return (double) clock() / (double) CLOCKS_PER_SEC;
}

static int popcount_ull(unsigned long long x) {
    return __builtin_popcountll(x);
}

static int bitset_count(const unsigned long long *bits, int nwords) {
    int total = 0;
    for (int i = 0; i < nwords; i++) {
        total += popcount_ull(bits[i]);
    }
    return total;
}

static int row_word(int row) {
    return row >> 6;
}

static unsigned long long row_mask(int row) {
    return 1ULL << (row & 63);
}

static int row_is_uncovered(const qca_scp_state *state, int row) {
    return (state->uncovered[row_word(row)] & row_mask(row)) != 0ULL;
}



static int row_conflict_packing_rec(
    const unsigned long long *adj,
    unsigned long long cand,
    int current,
    int *best
) {
    if (cand == 0ULL) {
        if (current > *best) {
            *best = current;
        }
        return *best;
    }

    if (current + popcount_ull(cand) <= *best) {
        return *best;
    }

    int v = 0;
    while (((cand >> v) & 1ULL) == 0ULL) {
        v++;
    }

    row_conflict_packing_rec(adj, cand & ~(1ULL << v), current, best);
    row_conflict_packing_rec(adj, cand & ~adj[v] & ~(1ULL << v), current + 1, best);

    return *best;
}

static int exact_row_packing_lb(
    const qca_scp_problem *problem,
    const qca_scp_state *state
) {
    const int col_words = (problem->nc + 63) / 64;
    unsigned long long *row_masks = NULL;
    unsigned long long *adj = NULL;
    int m = 0;
    int best = 0;

    if (problem->nr > 63) {
        return 0;
    }

    row_masks = (unsigned long long *) qca_scp_calloc(
        (size_t) problem->nr * col_words,
        sizeof(unsigned long long)
    );
    adj = (unsigned long long *) qca_scp_calloc(
        (size_t) problem->nr,
        sizeof(unsigned long long)
    );

    if (row_masks == NULL || adj == NULL) {
        qca_scp_free(row_masks);
        qca_scp_free(adj);
        return 0;
    }

    for (int row = 0; row < problem->nr; row++) {
        if (!row_is_uncovered(state, row)) {
            continue;
        }

        unsigned long long *mask = row_masks + (size_t) m * col_words;
        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (state->active[col]) {
                mask[col >> 6] |= (1ULL << (col & 63));
            }
        }
        m++;
    }

    for (int i = 0; i < m; i++) {
        const unsigned long long *mi = row_masks + (size_t) i * col_words;
        for (int j = i + 1; j < m; j++) {
            const unsigned long long *mj = row_masks + (size_t) j * col_words;
            for (int w = 0; w < col_words; w++) {
                if ((mi[w] & mj[w]) != 0ULL) {
                    adj[i] |= (1ULL << j);
                    adj[j] |= (1ULL << i);
                    break;
                }
            }
        }
    }

    row_conflict_packing_rec(adj, (m == 64 ? ~0ULL : ((1ULL << m) - 1ULL)), 0, &best);

    qca_scp_free(row_masks);
    qca_scp_free(adj);
    return best;
}

static void cover_with_column(
    const qca_scp_problem *problem,
    qca_scp_state *state,
    int col
) {
    state->chosen[col] = 1;
    state->active[col] = 0;
    state->chosen_count++;
    for (int w = 0; w < problem->nwords_rows; w++) {
        state->uncovered[w] &= ~problem->col_masks[col * problem->nwords_rows + w];
    }
}

static int new_rows_covered(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int col
) {
    int covered = 0;
    const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
    for (int w = 0; w < problem->nwords_rows; w++) {
        covered += popcount_ull(state->uncovered[w] & mask[w]);
    }
    return covered;
}

static int column_dominated(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int col_a,
    int col_b
) {
    const unsigned long long *a = problem->col_masks + col_a * problem->nwords_rows;
    const unsigned long long *b = problem->col_masks + col_b * problem->nwords_rows;

    for (int w = 0; w < problem->nwords_rows; w++) {
        unsigned long long only_a = state->uncovered[w] & a[w] & ~b[w];
        if (only_a != 0ULL) {
            return 0;
        }
    }
    return 1;
}

static int row_support_subset(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int row_sub,
    int row_super,
    unsigned int *marks,
    unsigned int stamp
) {
    int sub_support = 0;

    for (int i = problem->row_starts[row_super]; i < problem->row_starts[row_super + 1]; i++) {
        int col = problem->row_cols[i];
        if (state->active[col]) {
            marks[col] = stamp;
        }
    }

    for (int i = problem->row_starts[row_sub]; i < problem->row_starts[row_sub + 1]; i++) {
        int col = problem->row_cols[i];
        if (!state->active[col]) {
            continue;
        }
        sub_support++;
        if (marks[col] != stamp) {
            return 0;
        }
    }

    return sub_support > 0;
}

static int apply_reductions(
    const qca_scp_problem *problem,
    qca_scp_state *state,
    const qca_scp_search *search,
    int full_dominance
) {
    double started = now_seconds();
    int changed = 1;
    unsigned int *row_marks = search->reduction_row_marks;
    int *row_support = search->reduction_row_support;
    int *gains = search->reduction_gains;
    int *active_cols = search->reduction_active_cols;
    double *dual_slack = search->reduction_dual_slack;
    qca_profile.reduction_calls++;

    if (row_marks == NULL || row_support == NULL || gains == NULL ||
        active_cols == NULL || dual_slack == NULL) {
        qca_profile.reductions_seconds += now_seconds() - started;
        return 1;
    }

    memset(row_marks, 0, (size_t) problem->nc * sizeof(unsigned int));
    unsigned int row_mark_stamp = 0;

    while (changed) {
        changed = 0;
        int nactive = 0;

        for (int col = 0; col < problem->nc; col++) {
            gains[col] = 0;
            if (!state->active[col]) {
                continue;
            }
            gains[col] = new_rows_covered(problem, state, col);
            if (gains[col] == 0) {
                state->active[col] = 0;
                changed = 1;
                continue;
            }
            active_cols[nactive++] = col;
        }

        if (propagation_target(search) >= 0) {
            int slots = propagation_target(search) - state->chosen_count;
            int uncovered = bitset_count(state->uncovered, problem->nwords_rows);
            int max_gain = 0;
            for (int i = 0; i < nactive; ++i) {
                if (gains[active_cols[i]] > max_gain) max_gain = gains[active_cols[i]];
            }
            if (slots < 0 || (long long)slots * max_gain < uncovered) {
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }
            /* A selected column plus all remaining slots must cover every
               residual row. Apply this necessary condition before singleton
               and dominance propagation, not just when choosing branches. */
            long long minimum_gain = (long long)uncovered - (long long)(slots - 1) * max_gain;
            for (int i = 0; i < nactive; ++i) {
                int col = active_cols[i];
                if (gains[col] < minimum_gain) {
                    state->active[col] = 0;
                    changed = 1;
                }
            }
        }

        for (int row = 0; row < problem->nr; ++row) {
            int support = 0;
            if (row_is_uncovered(state, row)) {
                for (int i = problem->row_starts[row];
                     i < problem->row_starts[row + 1]; ++i) {
                    support += state->active[problem->row_cols[i]] != 0;
                }
            }
            row_support[row] = support;
        }

        int run_dominance = full_dominance || nactive <= 256;
        for (int ia = 0; run_dominance && ia < nactive; ia++) {
            int a = active_cols[ia];
            if (!state->active[a]) {
                continue;
            }

            /* Every column dominating a must cover every residual row of a.
               Scan only the candidates of a's rarest residual row rather
               than every active column. */
            int pivot_row = -1;
            int pivot_support = INT_MAX;
            const unsigned long long *a_mask =
                problem->col_masks + (size_t)a * problem->nwords_rows;
            for (int w = 0; w < problem->nwords_rows; ++w) {
                unsigned long long bits = state->uncovered[w] & a_mask[w];
                while (bits != 0ULL) {
                    int bit = __builtin_ctzll(bits);
                    int row = w * 64 + bit;
                    int support = row_support[row];
                    if (support < pivot_support) {
                        pivot_support = support;
                        pivot_row = row;
                    }
                    bits &= bits - 1ULL;
                }
            }

            if (pivot_row < 0) {
                continue;
            }

            for (int ib = problem->row_starts[pivot_row];
                 ib < problem->row_starts[pivot_row + 1]; ++ib) {
                int b = problem->row_cols[ib];
                if (a == b || !state->active[b]) {
                    continue;
                }
                if (gains[b] < gains[a]) {
                    continue;
                }
                if (column_dominated(problem, state, a, b)) {
                    state->active[a] = 0;
                    changed = 1;
                    break;
                }
            }
        }

        if (run_dominance) {
            /* Column dominance may have changed active supports. Refresh
               once; singleton detection can then reuse the counts. */
            for (int row = 0; row < problem->nr; ++row) {
                int support = 0;
                if (row_is_uncovered(state, row)) {
                    for (int i = problem->row_starts[row];
                         i < problem->row_starts[row + 1]; ++i) {
                        support += state->active[problem->row_cols[i]] != 0;
                    }
                }
                row_support[row] = support;
            }
        }

        for (int row = 0; row < problem->nr; row++) {
            if (!row_is_uncovered(state, row)) {
                continue;
            }

            int support = row_support[row];
            int forced_col = -1;
            if (support == 1) {
                support = 0;
                for (int i = problem->row_starts[row];
                     i < problem->row_starts[row + 1]; ++i) {
                    int col = problem->row_cols[i];
                    if (state->active[col]) {
                        ++support;
                        forced_col = col;
                        if (support > 1) break;
                    }
                }
            }

            if (support == 0) {
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }

            if (support == 1) {
                cover_with_column(problem, state, forced_col);
                changed = 1;
            }
        }

        if (search != NULL && propagation_target(search) >= 0) {
            int budget_left = propagation_target(search) - state->chosen_count;
            if (budget_left < 0) {
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }
            if (budget_left == 0 && bitset_count(state->uncovered, problem->nwords_rows) > 0) {
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }

            if (budget_left > 0) {
                double dual_lb = qca_scp_relaxation_dual_info(problem, state, dual_slack);
                double threshold = (double) budget_left - dual_lb + 1e-12;

                for (int col = 0; col < problem->nc; col++) {
                    if (!state->active[col]) {
                        continue;
                    }
                    if (dual_slack[col] > threshold) {
                        state->active[col] = 0;
                        changed = 1;
                    }
                }
            }
        }

        for (int a = 0; run_dominance && a < problem->nr; a++) {
            if (!row_is_uncovered(state, a)) {
                continue;
            }

            for (int b = 0; b < problem->nr; b++) {
                if (a == b || !row_is_uncovered(state, b)) {
                    continue;
                }

                ++row_mark_stamp;
                if (row_mark_stamp == 0) {
                    memset(row_marks, 0, (size_t) problem->nc * sizeof(unsigned int));
                    row_mark_stamp = 1;
                }

                if (row_support_subset(
                    problem, state, a, b, row_marks, row_mark_stamp
                )) {
                    state->uncovered[row_word(b)] &= ~row_mask(b);
                    changed = 1;
                }
            }
        }
    }

    qca_profile.reductions_seconds += now_seconds() - started;
    return 1;
}

/*
Chvatal-Gomory k-cover bound. A set S of uncovered rows such that every
active column covers at most k rows of S yields the valid inequality
sum(x) >= ceil(|S| / k) over the residual problem: summing the |S| row
constraints counts every active column at most k times. The k = 1 case is
the disjoint-row packing already computed elsewhere; k = 2 and 3 often cut
strictly deeper on dense charts, where large disjoint packings do not exist.
Greedy separation, admitting rows with small active support first.
*/
static int kcover_lb(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int k,
    int *hits,       /* nc scratch */
    int *rows_order, /* nr scratch */
    int *support     /* nr scratch */
) {
    const int nr = problem->nr;
    int n_unc = 0;

    if (k <= 1 || nr > 256) {
        return 0;
    }

    for (int row = 0; row < nr; row++) {
        if (!row_is_uncovered(state, row)) {
            continue;
        }
        int s = 0;
        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            if (state->active[problem->row_cols[i]]) {
                s++;
            }
        }
        support[row] = s;
        rows_order[n_unc++] = row;
    }

    if (n_unc == 0) {
        return 0;
    }

    /* ascending active support, ties by row id: deterministic */
    for (int i = 1; i < n_unc; i++) {
        int r = rows_order[i];
        int j = i - 1;
        while (j >= 0 && (support[rows_order[j]] > support[r] ||
               (support[rows_order[j]] == support[r] && rows_order[j] > r))) {
            rows_order[j + 1] = rows_order[j];
            j--;
        }
        rows_order[j + 1] = r;
    }

    memset(hits, 0, (size_t) problem->nc * sizeof(int));

    /*
    Connected components of the residual row-intersection graph (rows
    linked when an active column covers both). Components share no active
    columns, so the k-cover bound may be summed per component, which is
    strictly stronger than one global rounding: two disjoint odd holes
    give ceil(a/k) + ceil(b/k) instead of ceil((a+b)/k).
    Union-find over rows, path-halving.
    */
    int *comp = support; /* reuse: support values no longer needed */
    for (int i = 0; i < n_unc; i++) {
        comp[rows_order[i]] = rows_order[i];
    }
    for (int col = 0; col < problem->nc; col++) {
        if (!state->active[col]) {
            continue;
        }
        const unsigned long long *mask = problem->col_masks + (size_t) col * problem->nwords_rows;
        int first = -1;
        for (int w = 0; w < problem->nwords_rows; w++) {
            unsigned long long bits = mask[w] & state->uncovered[w];
            while (bits != 0ULL) {
                int b = 0;
                unsigned long long low = bits & (~bits + 1ULL);
                while (((low >> b) & 1ULL) == 0ULL) b++;
                int row = w * 64 + b;
                bits &= bits - 1ULL;

                if (first < 0) {
                    first = row;
                } else {
                    int ra = first, rb = row;
                    while (comp[ra] != ra) { comp[ra] = comp[comp[ra]]; ra = comp[ra]; }
                    while (comp[rb] != rb) { comp[rb] = comp[comp[rb]]; rb = comp[rb]; }
                    if (ra != rb) comp[rb] = ra;
                }
            }
        }
    }

    int in_s = 0;
    for (int idx = 0; idx < n_unc; idx++) {
        int row = rows_order[idx];
        int fits = 1;

        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (state->active[col] && hits[col] >= k) {
                fits = 0;
                break;
            }
        }
        if (!fits) {
            rows_order[idx] = -1; /* not in S */
            continue;
        }

        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (state->active[col]) {
                hits[col]++;
            }
        }
        in_s++;
    }

    if (in_s == 0) {
        return 0;
    }

    /* sum ceil(|S in component| / k) over components (nr <= 256) */
    int bound = 0;
    {
        int croots[257];
        int counted[257];
        int ncomp = 0;
        for (int idx = 0; idx < n_unc; idx++) {
            int row = rows_order[idx];
            if (row < 0) continue; /* not in S */
            int root = row;
            while (comp[root] != root) root = comp[root];
            int found = -1;
            for (int t = 0; t < ncomp; t++) {
                if (croots[t] == root) { found = t; break; }
            }
            if (found < 0) {
                croots[ncomp] = root;
                counted[ncomp] = 1;
                ncomp++;
            } else {
                counted[found]++;
            }
        }
        for (int t = 0; t < ncomp; t++) {
            bound += (counted[t] + k - 1) / k;
        }
    }

    return bound;
}

static int coverage_lower_bound(
    const qca_scp_problem *problem,
    const qca_scp_state *state
) {
    int uncovered = bitset_count(state->uncovered, problem->nwords_rows);
    int best_gain = 0;

    if (uncovered == 0) {
        return 0;
    }

    for (int col = 0; col < problem->nc; col++) {
        if (!state->active[col]) {
            continue;
        }
        int gain = new_rows_covered(problem, state, col);
        if (gain > best_gain) {
            best_gain = gain;
        }
    }

    if (best_gain <= 0) {
        return INT_MAX / 4;
    }
    return (uncovered + best_gain - 1) / best_gain;
}

static int lower_bound(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int max_useful_bound,
    int strong_bound
) {
    double started = now_seconds();
    int row_packing_lb = 0;
    int coverage_lb = coverage_lower_bound(problem, state);

    if (coverage_lb == 0 || coverage_lb >= INT_MAX / 4 ||
        (max_useful_bound >= 0 && coverage_lb > max_useful_bound) ||
        !strong_bound) {
        qca_profile.lower_bound_calls++;
        qca_profile.lower_bound_seconds += now_seconds() - started;
        return coverage_lb;
    }

    row_packing_lb = exact_row_packing_lb(problem, state);
    if (row_packing_lb <= 0) {
        unsigned char *blocked_cols = (unsigned char *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(unsigned char)
        );
        if (blocked_cols != NULL) {
            for (int row = 0; row < problem->nr; row++) {
                if (!row_is_uncovered(state, row)) {
                    continue;
                }

                int support = 0;
                int disjoint = 1;
                for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                    int col = problem->row_cols[i];
                    if (!state->active[col]) {
                        continue;
                    }
                    support++;
                    if (blocked_cols[col]) {
                        disjoint = 0;
                        break;
                    }
                }

                if (support == 0) {
                    qca_scp_free(blocked_cols);
                    qca_profile.lower_bound_calls++;
                    qca_profile.lower_bound_seconds += now_seconds() - started;
                    return INT_MAX / 4;
                }

                if (disjoint) {
                    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                        int col = problem->row_cols[i];
                        if (state->active[col]) {
                            blocked_cols[col] = 1;
                        }
                    }
                    row_packing_lb++;
                }
            }
            qca_scp_free(blocked_cols);
        }
    }

    if (row_packing_lb <= 0) {
        for (int row = 0; row < problem->nr; row++) {
            if (!row_is_uncovered(state, row)) {
                continue;
            }

            for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                int col = problem->row_cols[i];
                if (state->active[col]) {
                    row_packing_lb = 1;
                    break;
                }
            }
            if (row_packing_lb > 0) {
                break;
            }
        }
    }

    int relaxation_lb = qca_scp_relaxation_dual_lb(problem, state);
    if (relaxation_lb > row_packing_lb) {
        row_packing_lb = relaxation_lb;
    }

    /* Chvatal-Gomory k-cover cuts, the one bound family able to close
       integrality gaps that the dual (LP-strength) bound cannot reach */
    if (problem->nr <= 256) {
        int *hits = (int *) qca_scp_calloc((size_t) problem->nc, sizeof(int));
        int *rows_order = (int *) qca_scp_calloc((size_t) problem->nr, sizeof(int));
        int *support_s = (int *) qca_scp_calloc((size_t) problem->nr, sizeof(int));
        if (hits != NULL && rows_order != NULL && support_s != NULL) {
            for (int k = 2; k <= 3; k++) {
                int kc = kcover_lb(problem, state, k, hits, rows_order, support_s);
                if (kc > row_packing_lb) {
                    row_packing_lb = kc;
                }
            }
        }
        qca_scp_free(hits);
        qca_scp_free(rows_order);
        qca_scp_free(support_s);
    }

    qca_profile.lower_bound_calls++;
    qca_profile.lower_bound_seconds += now_seconds() - started;
    return (row_packing_lb > coverage_lb ? row_packing_lb : coverage_lb);
}



static int branch_col_score(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    const qca_scp_search *search,
    int col
) {
    int score = 0;
    const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
    /* Read the hybrid heuristic: starts as active row support, becomes incumbent coverage */
    const int *cached_support = search->branching_row_heuristic;
    for (int w = 0; w < problem->nwords_rows; w++) {
        unsigned long long bits = mask[w] & state->uncovered[w];
        while (bits) {
            int bit = __builtin_ctzll(bits);
            int row = w * 64 + bit;

            int support = cached_support[row];
            if (support > 0) {
                score += 1024 / support;
            }
            bits &= bits - 1ULL;
        }
    }

    if (problem->branch_priority != NULL) {
        double p = problem->branch_priority[col];
        if (p < 0.0) {
            score += (int) (-1000.0 * p);
        }
    }

    return score;
}

static int choose_branch_row(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    const qca_scp_search *search,
    const int *gains,
    int minimum_gain
) {
    int best_row = -1;
    int best_support = INT_MAX;
    int best_row_score = -1;

    for (int w = 0; w < problem->nwords_rows; w++) {
        unsigned long long bits = state->uncovered[w];
        while (bits) {
            int bit = __builtin_ctzll(bits);
            int row = w * 64 + bit;
            bits &= bits - 1ULL;

            int support = 0;
            for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                int col = problem->row_cols[i];
                if (state->active[col] && gains[col] >= minimum_gain) {
                    support++;
                }
            }

            if (support <= best_support) {
                int row_score = 0;
                for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                    int col = problem->row_cols[i];
                    if (state->active[col] && gains[col] >= minimum_gain) {
                        row_score += branch_col_score(problem, state, search, col);
                    }
                }

                if (support < best_support || (support == best_support && row_score > best_row_score)) {
                    best_support = support;
                    best_row_score = row_score;
                    best_row = row;
                    if (support <= 2) {
                        return best_row;
                    }
                }
            }
        }
    }

    return best_row;
}

static void copy_state(
    const qca_scp_problem *problem,
    const qca_scp_state *src,
    qca_scp_state *dst
) {
    memcpy(dst->uncovered, src->uncovered, (size_t) problem->nwords_rows * sizeof(unsigned long long));
    memcpy(dst->active, src->active, (size_t) problem->nc * sizeof(unsigned char));
    memcpy(dst->chosen, src->chosen, (size_t) problem->nc * sizeof(unsigned char));
    dst->chosen_count = src->chosen_count;
}

static void sort_branch_cols(int *branch_cols, int nbranch, const int *gains, const int *scores) {
    for (int i = 1; i < nbranch; i++) {
        int col = branch_cols[i];
        int j = i - 1;
        while (j >= 0) {
            int prev = branch_cols[j];
            if (gains[prev] > gains[col] ||
                (gains[prev] == gains[col] && scores[prev] > scores[col]) ||
                (scores[prev] == scores[col] && gains[prev] == gains[col] && prev < col)) {
                break;
            }
            branch_cols[j + 1] = prev;
            j--;
        }
        branch_cols[j + 1] = col;
    }
}

static int greedy_cover_in_place(
    const qca_scp_problem *problem,
    qca_scp_state *state,
    int *candidate_solution
) {
    double started = now_seconds();
    while (bitset_count(state->uncovered, problem->nwords_rows) > 0) {
        int best_col = -1;
        int best_gain = 0;

        for (int col = 0; col < problem->nc; col++) {
            if (!state->active[col]) {
                continue;
            }
            int gain = 0;
            const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
            for (int w = 0; w < problem->nwords_rows; w++) {
                gain += popcount_ull(state->uncovered[w] & mask[w]);
            }
            if (gain > best_gain ||
                (gain == best_gain && best_col >= 0 &&
                 problem->branch_priority != NULL &&
                 problem->branch_priority[col] < problem->branch_priority[best_col])) {
                best_gain = gain;
                best_col = col;
            }
        }

        if (best_col < 0 || best_gain == 0) {
            qca_profile.greedy_calls++;
            qca_profile.greedy_seconds += now_seconds() - started;
            return INT_MAX / 4;
        }

        cover_with_column(problem, state, best_col);
    }

    for (int col = 0; col < problem->nc; col++) {
        candidate_solution[col] = state->chosen[col];
    }

    qca_profile.greedy_calls++;
    qca_profile.greedy_seconds += now_seconds() - started;
    return state->chosen_count;
}

static int cleanup_solution(
    const qca_scp_problem *problem,
    int *solution,
    int *row_coverage
) {
    int changed = 1;
    int size = 0;

    memset(row_coverage, 0, (size_t) problem->nr * sizeof(int));
    for (int col = 0; col < problem->nc; col++) {
        if (!solution[col]) {
            continue;
        }
        size++;
        const unsigned long long *mask =
            problem->col_masks + (size_t) col * problem->nwords_rows;
        for (int row = 0; row < problem->nr; row++) {
            if ((mask[row_word(row)] & row_mask(row)) != 0ULL) {
                row_coverage[row]++;
            }
        }
    }

    while (changed) {
        changed = 0;
        for (int col = 0; col < problem->nc; col++) {
            if (!solution[col]) {
                continue;
            }

            int redundant = 1;
            const unsigned long long *mask =
                problem->col_masks + (size_t) col * problem->nwords_rows;
            for (int row = 0; row < problem->nr; row++) {
                if ((mask[row_word(row)] & row_mask(row)) != 0ULL &&
                    row_coverage[row] <= 1) {
                    redundant = 0;
                    break;
                }
            }
            if (!redundant) {
                continue;
            }

            solution[col] = 0;
            size--;
            changed = 1;
            for (int row = 0; row < problem->nr; row++) {
                if ((mask[row_word(row)] & row_mask(row)) != 0ULL) {
                    row_coverage[row]--;
                }
            }
        }
    }

    return size;
}

static int keep_greedy_candidate(
    qca_scp_search *search,
    int *candidate_solution,
    int candidate_size
) {
    if (candidate_size >= INT_MAX / 4) {
        return 0;
    }

    candidate_size = cleanup_solution(
        search->problem,
        candidate_solution,
        search->branching_row_heuristic
    );
    if (candidate_size >= search->best_size) {
        return 0;
    }

    search->best_size = candidate_size;
    memcpy(
        search->best_solution,
        candidate_solution,
        (size_t) search->problem->nc * sizeof(int)
    );
    if (search->proof_target_size >= 0 &&
        candidate_size <= search->proof_target_size) {
        search->found_target = 1;
        return 1;
    }
    return 0;
}

static int search_limit_is_reached(qca_scp_search *search) {
    double elapsed = now_seconds() - search->started;
    int hard_node_reached = search->hard_node_limit > 0 &&
        search->nodes >= search->hard_node_limit;
    int hard_time_reached = search->hard_time_limit_seconds > 0.0 &&
        elapsed >= search->hard_time_limit_seconds;
    if (hard_node_reached || hard_time_reached) {
        return 1;
    }

    int probe_node_reached = search->node_limit > 0 &&
        search->nodes >= search->node_limit;
    int probe_time_reached = search->time_limit_seconds > 0.0 &&
        elapsed >= search->time_limit_seconds;
    if (!probe_node_reached && !probe_time_reached) {
        return 0;
    }
    if (!search->adaptive_limits) {
        return 1;
    }

    int progress = search->best_size < search->checkpoint_best_size ||
        search->root_branches_completed > search->checkpoint_root_branches;
    if (!progress) {
        return 1;
    }

    search->checkpoint_best_size = search->best_size;
    search->checkpoint_root_branches = search->root_branches_completed;
    qca_profile.adaptive_extensions++;

    if (search->node_limit > 0) {
        unsigned long long extended = search->node_limit > ULLONG_MAX / 3ULL
            ? ULLONG_MAX
            : search->node_limit * 3ULL;
        if (search->hard_node_limit > 0 && extended > search->hard_node_limit) {
            extended = search->hard_node_limit;
        }
        search->node_limit = extended;
    }
    if (search->time_limit_seconds > 0.0) {
        double extended = search->time_limit_seconds * 3.0;
        if (
            search->hard_time_limit_seconds > 0.0 &&
            extended > search->hard_time_limit_seconds
        ) {
            extended = search->hard_time_limit_seconds;
        }
        search->time_limit_seconds = extended;
    }
    return 0;
}

static void search_exact(
    qca_scp_search *search,
    qca_scp_state *state,
    int depth
) {
    const qca_scp_problem *problem = search->problem;
    qca_scp_workspace *workspace = search->workspaces + depth;

    if (search->found_target || search->limit_reached) {
        return;
    }

    if (search_limit_is_reached(search)) {
        search->limit_reached = 1;
        return;
    }

    search->nodes++;
    qca_profile.nodes++;
    if ((search->nodes & 1023ULL) == 0ULL) {
        R_CheckUserInterrupt();
    }

    if (propagation_target(search) >= 0) {
        int fast_lb = coverage_lower_bound(problem, state);
        if (state->chosen_count + fast_lb > propagation_target(search)) {
            qca_profile.leaves++;
            return;
        }
    }

    if (!apply_reductions(problem, state, search, depth == 0)) {
        qca_profile.leaves++;
        return;
    }

    if (state->chosen_count >= search->best_size) {
        qca_profile.leaves++;
        return;
    }

    if (problem->lagrangian_iterations > 0) {
        double bound_started = now_seconds();
        int fixed_columns = 0;
        int target = search->best_size - state->chosen_count - 1;
        if (search->proof_target_size >= 0 &&
            search->proof_target_size - state->chosen_count < target) {
            target = search->proof_target_size - state->chosen_count;
        }
        int residual_lb = qca_scp_lagrangian_lb(
            problem, state, workspace->dual, target,
            depth % 3 == 0 ? problem->lagrangian_iterations : 1,
            search->lagrangian_costs, search->gradient, search->best_dual,
            &fixed_columns
        );
        ++qca_profile.lagrangian_calls;
        qca_profile.lagrangian_fixed_columns += fixed_columns;
        qca_profile.lagrangian_seconds += now_seconds() - bound_started;
        if (residual_lb > target) {
            ++qca_profile.lagrangian_prunes;
            ++qca_profile.leaves;
            return;
        }
    }

    if (search->best_size > state->chosen_count + 1) {
        int *candidate_solution = workspace->candidate_solution;
        if (candidate_solution != NULL) {
            qca_scp_state *greedy_state =
                &search->workspaces[depth + 1].state;
            copy_state(problem, state, greedy_state);
            int greedy_size = greedy_cover_in_place(
                problem,
                greedy_state,
                candidate_solution
            );
            if (keep_greedy_candidate(search, candidate_solution, greedy_size)) {
                qca_profile.leaves++;
                return;
            }
        }
    }

    if (bitset_count(state->uncovered, problem->nwords_rows) == 0) {
        qca_profile.leaves++;
        search->best_size = state->chosen_count;
        for (int c = 0; c < problem->nc; c++) {
            search->best_solution[c] = state->chosen[c];
        }
        if (search->proof_target_size >= 0 && state->chosen_count <= search->proof_target_size) {
            search->found_target = 1;
        }
        return;
    }

    int max_useful_bound = search->best_size - state->chosen_count - 1;
    if (search->proof_target_size >= 0) {
        int proof_limit = search->proof_target_size - state->chosen_count;
        if (proof_limit < max_useful_bound) {
            max_useful_bound = proof_limit;
        }
    }
    int lb = lower_bound(
        problem,
        state,
        max_useful_bound,
        depth == 0 || depth % 3 == 0
    );
    if (search->proof_target_size >= 0 &&
        state->chosen_count + lb > search->proof_target_size) {
        qca_profile.leaves++;
        return;
    }
    if (state->chosen_count + lb >= search->best_size) {
        qca_profile.leaves++;
        return;
    }

    int *branch_cols = workspace->branch_cols;
    int *gains = workspace->gains;
    int *scores = workspace->candidate_solution;
    if (branch_cols == NULL || gains == NULL || scores == NULL) {
        return;
    }

    int max_gain = 0;
    for (int col = 0; col < problem->nc; col++) {
        gains[col] = state->active[col]
            ? new_rows_covered(problem, state, col)
            : 0;
        if (gains[col] > max_gain) {
            max_gain = gains[col];
        }
    }

    int columns_left = search->proof_target_size >= 0
        ? search->proof_target_size - state->chosen_count
        : search->best_size - state->chosen_count - 1;
    int uncovered_count = bitset_count(
        state->uncovered,
        problem->nwords_rows
    );
    long long required_gain = (long long) uncovered_count -
        (long long) (columns_left - 1) * max_gain;
    int minimum_gain = required_gain > 1 ? (int) required_gain : 1;

    int row = choose_branch_row(problem, state, search, gains, minimum_gain);
    if (row < 0) {
        return;
    }

    int nbranch = 0;
    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
        int col = problem->row_cols[i];
        if (state->active[col] && gains[col] >= minimum_gain) {
            branch_cols[nbranch++] = col;
            scores[col] = branch_col_score(problem, state, search, col);
        }
    }

    {
        double started = now_seconds();
        sort_branch_cols(branch_cols, nbranch, gains, scores);
        qca_profile.branching_seconds += now_seconds() - started;
    }

    if (depth == 0) {
        search->root_branches_total = (unsigned long long)nbranch;
        qca_profile.root_branches_total += (unsigned long long)nbranch;
    }

    for (int i = 0; i < nbranch; i++) {
        int col = branch_cols[i];
        qca_scp_state *child = &(search->workspaces[depth + 1].state);

        /* Partition the search space: sibling i contains solutions that
           select branch_cols[i] and none of the earlier alternatives for
           this row. Without these exclusions, the same column set is
           revisited in many selection orders. */
        copy_state(problem, state, child);
        cover_with_column(problem, child, col);
        if (problem->lagrangian_iterations > 0) {
            memcpy(search->workspaces[depth + 1].dual, workspace->dual,
                (size_t)problem->nr * sizeof(int64_t));
        }
        search_exact(search, child, depth + 1);
        if (search->limit_reached) {
            for (int j = 0; j < nbranch; j++) {
                state->active[branch_cols[j]] = 1;
            }
            return;
        }
        if (depth == 0) {
            search->root_branches_completed++;
            qca_profile.root_branches_completed++;
        }
        if (search->found_target) {
            for (int j = 0; j < nbranch; j++) {
                state->active[branch_cols[j]] = 1;
            }
            return;
        }
        state->active[col] = 0;
    }

    for (int i = 0; i < nbranch; i++) {
        state->active[branch_cols[i]] = 1;
    }
}

typedef struct {
    qca_scp_search *search;
    qca_scp_state *state;
    int max_depth;
} qca_scp_run;

static SEXP run_exact_search(void *data) {
    qca_scp_run *run = (qca_scp_run *)data;
    search_exact(run->search, run->state, 0);
    return R_NilValue;
}

static void free_exact_search(void *data) {
    qca_scp_run *run = (qca_scp_run *)data;
    qca_scp_free(run->state->uncovered);
    qca_scp_free(run->state->active);
    qca_scp_free(run->state->chosen);
    for (int depth = 0; depth < run->max_depth; depth++) {
        qca_scp_free(run->search->workspaces[depth].state.uncovered);
        qca_scp_free(run->search->workspaces[depth].state.active);
        qca_scp_free(run->search->workspaces[depth].state.chosen);
        qca_scp_free(run->search->workspaces[depth].branch_cols);
        qca_scp_free(run->search->workspaces[depth].gains);
        qca_scp_free(run->search->workspaces[depth].candidate_solution);
        qca_scp_free(run->search->workspaces[depth].dual);
    }
    qca_scp_free(run->search->workspaces);
    qca_scp_free(run->search->reduction_row_marks);
    qca_scp_free(run->search->reduction_row_support);
    qca_scp_free(run->search->reduction_gains);
    qca_scp_free(run->search->reduction_active_cols);
    qca_scp_free(run->search->reduction_dual_slack);
    qca_scp_free(run->search->branching_row_heuristic);
    qca_scp_free(run->search->lagrangian_costs);
    qca_scp_free(run->search->best_dual);
    qca_scp_free(run->search->gradient);
}

qca_scp_result qca_scp_solve_exact(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size
) {
    return qca_scp_solve_exact_with_incumbent(problem, solution, solution_size, NULL, 0, -1);
}

qca_scp_result qca_scp_solve_exact_with_incumbent(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size
) {
    return qca_scp_solve_exact_with_incumbent_bounded(
        problem,
        solution,
        solution_size,
        initial_solution,
        initial_size,
        proof_target_size,
        0,
        0.0
    );
}

qca_scp_result qca_scp_solve_exact_with_incumbent_bounded(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size,
    unsigned long long node_limit,
    double time_limit_seconds
) {
    return qca_scp_solve_exact_with_incumbent_adaptive(
        problem,
        solution,
        solution_size,
        initial_solution,
        initial_size,
        proof_target_size,
        node_limit,
        time_limit_seconds,
        node_limit,
        time_limit_seconds
    );
}

qca_scp_result qca_scp_solve_exact_with_incumbent_adaptive(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size,
    unsigned long long probe_node_limit,
    double probe_time_limit_seconds,
    unsigned long long hard_node_limit,
    double hard_time_limit_seconds
) {
    double started = now_seconds();
    qca_scp_state state = {0};
    qca_scp_search search = {0};
    qca_scp_workspace *workspaces = NULL;

    if (solution_size == NULL) {
        return QCA_SCP_ERROR;
    }
    *solution_size = 0;

    if (
        problem == NULL ||
        solution == NULL ||
        problem->nr <= 0 ||
        problem->nc <= 0 ||
        problem->nwords_rows <= 0 ||
        problem->row_starts == NULL ||
        problem->row_cols == NULL ||
        problem->col_masks == NULL
    ) {
        return QCA_SCP_ERROR;
    }

    /*
    Search depth is bounded by the incumbent / proof target: each branching
    level adds at least one chosen column, and nodes with chosen_count
    beyond the bound are pruned before branching. The former nc+1 sizing
    is kept only as a fallback (and upper cap) when no bound is known.
    */
    int max_depth = problem->nc + 1;
    if (proof_target_size >= 0) {
        max_depth = proof_target_size + 2;
    } else if (initial_solution != NULL && initial_size > 0) {
        max_depth = initial_size + 2;
    }
    if (max_depth > problem->nc + 1) {
        max_depth = problem->nc + 1;
    }

    state.uncovered = (unsigned long long *) qca_scp_calloc(
        (size_t) problem->nwords_rows,
        sizeof(unsigned long long)
    );
    state.active = (unsigned char *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(unsigned char)
    );
    state.chosen = (unsigned char *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(unsigned char)
    );
    workspaces = (qca_scp_workspace *) qca_scp_calloc(
        (size_t) max_depth,
        sizeof(qca_scp_workspace)
    );
    search.reduction_row_marks = (unsigned int *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(unsigned int)
    );
    search.reduction_row_support = (int *) qca_scp_calloc(
        (size_t) problem->nr,
        sizeof(int)
    );
    search.reduction_gains = (int *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(int)
    );
    search.reduction_active_cols = (int *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(int)
    );
    search.reduction_dual_slack = (double *) qca_scp_calloc(
        (size_t) problem->nc,
        sizeof(double)
    );
    search.branching_row_heuristic = (int *) qca_scp_calloc(
        (size_t) problem->nr,
        sizeof(int)
    );

    if (problem->lagrangian_iterations > 0) {
        search.lagrangian_costs = qca_scp_calloc((size_t)problem->nc, sizeof(int64_t));
        search.best_dual = qca_scp_calloc((size_t)problem->nr, sizeof(int64_t));
        search.gradient = qca_scp_calloc((size_t)problem->nr, sizeof(int));
    }

    if (state.uncovered == NULL || state.active == NULL || state.chosen == NULL ||
        workspaces == NULL || search.reduction_row_marks == NULL ||
        search.reduction_row_support == NULL || search.reduction_gains == NULL ||
        search.reduction_active_cols == NULL ||
        search.reduction_dual_slack == NULL ||
        search.branching_row_heuristic == NULL ||
        (problem->lagrangian_iterations > 0 && (!search.lagrangian_costs ||
            !search.best_dual || !search.gradient))) {
        qca_scp_free(state.uncovered);
        qca_scp_free(state.active);
        qca_scp_free(state.chosen);
        qca_scp_free(workspaces);
        qca_scp_free(search.reduction_row_marks);
        qca_scp_free(search.reduction_gains);
        qca_scp_free(search.reduction_active_cols);
        qca_scp_free(search.reduction_dual_slack);
        qca_scp_free(search.reduction_row_support);
        qca_scp_free(search.branching_row_heuristic);
        qca_scp_free(search.lagrangian_costs);
        qca_scp_free(search.best_dual);
        qca_scp_free(search.gradient);
        return QCA_SCP_ERROR;
    }

    for (int depth = 0; depth < max_depth; depth++) {
        workspaces[depth].state.uncovered = (unsigned long long *) qca_scp_calloc(
            (size_t) problem->nwords_rows,
            sizeof(unsigned long long)
        );
        workspaces[depth].state.active = (unsigned char *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(unsigned char)
        );
        workspaces[depth].state.chosen = (unsigned char *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(unsigned char)
        );
        workspaces[depth].branch_cols = (int *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(int)
        );
        workspaces[depth].gains = (int *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(int)
        );
        if (problem->lagrangian_iterations > 0) {
            workspaces[depth].dual = qca_scp_calloc((size_t)problem->nr, sizeof(int64_t));
        }
        workspaces[depth].candidate_solution = (int *) qca_scp_calloc(
            (size_t) problem->nc,
            sizeof(int)
        );

        if (workspaces[depth].state.uncovered == NULL ||
            workspaces[depth].state.active == NULL ||
            workspaces[depth].state.chosen == NULL ||
            workspaces[depth].branch_cols == NULL ||
            workspaces[depth].gains == NULL ||
            workspaces[depth].candidate_solution == NULL ||
            (problem->lagrangian_iterations > 0 && !workspaces[depth].dual)) {
            for (int i = 0; i <= depth; i++) {
                qca_scp_free(workspaces[i].state.uncovered);
                qca_scp_free(workspaces[i].state.active);
                qca_scp_free(workspaces[i].state.chosen);
                qca_scp_free(workspaces[i].branch_cols);
                qca_scp_free(workspaces[i].gains);
                qca_scp_free(workspaces[i].candidate_solution);
                qca_scp_free(workspaces[i].dual);
            }
            qca_scp_free(state.uncovered);
            qca_scp_free(state.active);
            qca_scp_free(state.chosen);
            qca_scp_free(workspaces);
            qca_scp_free(search.reduction_row_marks);
            qca_scp_free(search.reduction_row_support);
            qca_scp_free(search.reduction_gains);
            qca_scp_free(search.reduction_active_cols);
            qca_scp_free(search.reduction_dual_slack);
            qca_scp_free(search.branching_row_heuristic);
            qca_scp_free(search.lagrangian_costs);
            qca_scp_free(search.best_dual);
            qca_scp_free(search.gradient);
            return QCA_SCP_ERROR;
        }
    }

    for (int row = 0; row < problem->nr; row++) {
        state.uncovered[row_word(row)] |= row_mask(row);
    }
    memset(state.active, 1, (size_t) problem->nc * sizeof(unsigned char));
    memset(state.chosen, 0, (size_t) problem->nc * sizeof(unsigned char));
    state.chosen_count = 0;

    if (problem->lagrangian_iterations > 0) {
        for (int row = 0; row < problem->nr; ++row) {
            double y = problem->initial_row_dual ? problem->initial_row_dual[row] : 0.0;
            workspaces[0].dual[row] = !isfinite(y) || y <= 0.0 ? 0 : y >= 1.0
                ? QCA_SCP_DUAL_SCALE : (int64_t)(y * QCA_SCP_DUAL_SCALE);
        }
    }

    search.problem = problem;
    search.best_solution = solution;
    search.best_size = problem->nc + 1;
    search.proof_target_size = proof_target_size;
    search.found_target = 0;
    search.limit_reached = 0;
    search.nodes = 0;
    search.node_limit = probe_node_limit;
    search.hard_node_limit = hard_node_limit;
    search.started = started;
    search.time_limit_seconds = probe_time_limit_seconds;
    search.hard_time_limit_seconds = hard_time_limit_seconds;
    search.adaptive_limits =
        (hard_node_limit == 0 || hard_node_limit > probe_node_limit) ||
        (
            hard_time_limit_seconds <= 0.0 ||
            hard_time_limit_seconds > probe_time_limit_seconds
        );
    search.workspaces = workspaces;
    memset(solution, 0, (size_t) problem->nc * sizeof(int));

    if (initial_solution != NULL &&
        initial_size > 0 &&
        initial_size <= problem->nc) {
        memcpy(solution, initial_solution, (size_t) problem->nc * sizeof(int));
        search.best_size = initial_size;
        if (proof_target_size >= 0 && initial_size <= proof_target_size) {
            search.found_target = 1;
        }
    }
    search.checkpoint_best_size = search.best_size;

    qca_scp_run run = {&search, &state, max_depth};
    /* An unlimited proof can be interrupted from R. Release its recursive
       workspaces on that exit as well as on normal completion. */
    R_ExecWithCleanup(run_exact_search, &run, free_exact_search, &run);


    if (search.limit_reached) {
        if (search.best_size <= problem->nc) {
            *solution_size = search.best_size;
        }
        qca_profile.total_seconds += now_seconds() - started;
        return QCA_SCP_LIMIT;
    }

    if (proof_target_size >= 0 && !search.found_target) {
        qca_profile.total_seconds += now_seconds() - started;
        return QCA_SCP_NO_SOLUTION;
    }

    if (search.best_size > problem->nc) {
        qca_profile.total_seconds += now_seconds() - started;
        return QCA_SCP_NO_SOLUTION;
    }

    *solution_size = search.best_size;
    qca_profile.total_seconds += now_seconds() - started;
    return QCA_SCP_SOLUTION;
}

void qca_scp_profile_reset(void) {
    memset(&qca_profile, 0, sizeof(qca_profile));
}

qca_scp_profile qca_scp_profile_get(void) {
    return qca_profile;
}
