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

#include "scp_solver.h"
#include "scp_relaxation.h"
#include <limits.h>
#include <R_ext/RS.h>
#include <time.h>
#include <string.h>
typedef struct {
    qca_scp_state state;
    int *branch_cols;
    int *gains;
    int *candidate_solution;
} qca_scp_workspace;
typedef struct {
    const qca_scp_problem *problem;
    int *best_solution;
    int best_size;
    int proof_target_size;
    int found_target;
    qca_scp_workspace *workspaces;
} qca_scp_search;
static qca_scp_profile qca_profile = {0};
static double now_seconds(void) {
    return (double) clock() / (double) CLOCKS_PER_SEC;
}
static int popcount_ull(unsigned long long x) {
    int n = 0;
    while (x != 0ULL) {
        x &= (x - 1ULL);
        n++;
    }
    return n;
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
    row_masks = (unsigned long long *) R_Calloc((size_t) problem->nr * col_words, unsigned long long);
    adj = (unsigned long long *) R_Calloc((size_t) problem->nr, unsigned long long);
    if (row_masks == NULL || adj == NULL) {
        if (row_masks) R_Free(row_masks);
        if (adj) R_Free(adj);
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
    R_Free(row_masks);
    R_Free(adj);
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
    unsigned char *marks
) {
    int sub_support = 0;
    memset(marks, 0, (size_t) problem->nc * sizeof(unsigned char));
    for (int i = problem->row_starts[row_super]; i < problem->row_starts[row_super + 1]; i++) {
        int col = problem->row_cols[i];
        if (state->active[col]) {
            marks[col] = 1;
        }
    }
    for (int i = problem->row_starts[row_sub]; i < problem->row_starts[row_sub + 1]; i++) {
        int col = problem->row_cols[i];
        if (!state->active[col]) {
            continue;
        }
        sub_support++;
        if (!marks[col]) {
            return 0;
        }
    }
    return sub_support > 0;
}
static int apply_reductions(
    const qca_scp_problem *problem,
    qca_scp_state *state,
    const qca_scp_search *search
) {
    double started = now_seconds();
    int changed = 1;
    unsigned char *row_marks = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
    int *gains = (int *) R_Calloc((size_t) problem->nc, int);
    int *active_cols = (int *) R_Calloc((size_t) problem->nc, int);
    double *dual_slack = NULL;
    qca_profile.reduction_calls++;
    if (row_marks == NULL || gains == NULL || active_cols == NULL) {
        if (row_marks) R_Free(row_marks);
        if (gains) R_Free(gains);
        if (active_cols) R_Free(active_cols);
        qca_profile.reductions_seconds += now_seconds() - started;
        return 1;
    }
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
        for (int ia = 0; ia < nactive; ia++) {
            int a = active_cols[ia];
            if (!state->active[a]) {
                continue;
            }
            for (int ib = 0; ib < nactive; ib++) {
                int b = active_cols[ib];
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
        for (int row = 0; row < problem->nr; row++) {
            if (!row_is_uncovered(state, row)) {
                continue;
            }
            int support = 0;
            int forced_col = -1;
            for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
                int col = problem->row_cols[i];
                if (state->active[col]) {
                    support++;
                    forced_col = col;
                    if (support > 1) {
                        break;
                    }
                }
            }
            if (support == 0) {
                return 0;
            }
            if (support == 1) {
                cover_with_column(problem, state, forced_col);
                changed = 1;
            }
        }
        if (search != NULL && search->proof_target_size >= 0) {
            int budget_left = search->proof_target_size - state->chosen_count;
            if (budget_left < 0) {
                R_Free(active_cols);
                R_Free(gains);
                R_Free(row_marks);
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }
            if (budget_left == 0 && bitset_count(state->uncovered, problem->nwords_rows) > 0) {
                R_Free(active_cols);
                R_Free(gains);
                R_Free(row_marks);
                qca_profile.reductions_seconds += now_seconds() - started;
                return 0;
            }
            dual_slack = (double *) R_Calloc((size_t) problem->nc, double);
            if (dual_slack != NULL && budget_left > 0) {
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
                R_Free(dual_slack);
                dual_slack = NULL;
            }
        }
        for (int a = 0; a < problem->nr; a++) {
            if (!row_is_uncovered(state, a)) {
                continue;
            }
            for (int b = 0; b < problem->nr; b++) {
                if (a == b || !row_is_uncovered(state, b)) {
                    continue;
                }
                if (row_support_subset(problem, state, a, b, row_marks)) {
                    state->uncovered[row_word(b)] &= ~row_mask(b);
                    changed = 1;
                }
            }
        }
    }
    R_Free(active_cols);
    R_Free(gains);
    R_Free(row_marks);
    qca_profile.reductions_seconds += now_seconds() - started;
    return 1;
}
static int lower_bound(
    const qca_scp_problem *problem,
    const qca_scp_state *state
) {
    double started = now_seconds();
    int uncovered = bitset_count(state->uncovered, problem->nwords_rows);
    int best_gain = 0;
    int row_packing_lb = 0;
    if (uncovered == 0) {
        qca_profile.lower_bound_calls++;
        qca_profile.lower_bound_seconds += now_seconds() - started;
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
        qca_profile.lower_bound_calls++;
        qca_profile.lower_bound_seconds += now_seconds() - started;
        return INT_MAX / 4;
    }
    row_packing_lb = exact_row_packing_lb(problem, state);
    if (row_packing_lb <= 0) {
        unsigned char *blocked_cols = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
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
                    R_Free(blocked_cols);
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
            R_Free(blocked_cols);
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
    int coverage_lb = (uncovered + best_gain - 1) / best_gain;
    int relaxation_lb = qca_scp_relaxation_dual_lb(problem, state);
    if (relaxation_lb > row_packing_lb) {
        row_packing_lb = relaxation_lb;
    }
    qca_profile.lower_bound_calls++;
    qca_profile.lower_bound_seconds += now_seconds() - started;
    return (row_packing_lb > coverage_lb ? row_packing_lb : coverage_lb);
}
static int active_row_support(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int row
) {
    int support = 0;
    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
        if (state->active[problem->row_cols[i]]) {
            support++;
        }
    }
    return support;
}
static int branch_col_score(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int col
) {
    int score = 0;
    const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
    for (int row = 0; row < problem->nr; row++) {
        if (!row_is_uncovered(state, row)) {
            continue;
        }
        const int word = row_word(row);
        const unsigned long long bit = row_mask(row);
        if ((mask[word] & bit) == 0ULL) {
            continue;
        }
        int support = active_row_support(problem, state, row);
        if (support > 0) {
            score += 1024 / support;
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
    const qca_scp_state *state
) {
    int best_row = -1;
    int best_support = INT_MAX;
    int best_row_score = -1;
    for (int row = 0; row < problem->nr; row++) {
        if (!row_is_uncovered(state, row)) {
            continue;
        }
        int support = 0;
        int row_score = 0;
        for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
            int col = problem->row_cols[i];
            if (state->active[col]) {
                support++;
                row_score += branch_col_score(problem, state, col);
            }
        }
        if (support < best_support || (support == best_support && row_score > best_row_score)) {
            best_support = support;
            best_row_score = row_score;
            best_row = row;
            if (support <= 2) {
                break;
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
            if (scores[prev] > scores[col] ||
                (scores[prev] == scores[col] && gains[prev] > gains[col]) ||
                (scores[prev] == scores[col] && gains[prev] == gains[col] && prev < col)) {
                break;
            }
            branch_cols[j + 1] = prev;
            j--;
        }
        branch_cols[j + 1] = col;
    }
}
static int greedy_cover(
    const qca_scp_problem *problem,
    const qca_scp_state *state,
    int *candidate_solution
) {
    double started = now_seconds();
    unsigned long long *uncovered = (unsigned long long *) R_Calloc((size_t) problem->nwords_rows, unsigned long long);
    unsigned char *active = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
    unsigned char *chosen = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
    if (uncovered == NULL || active == NULL || chosen == NULL) {
        if (uncovered) R_Free(uncovered);
        if (active) R_Free(active);
        if (chosen) R_Free(chosen);
        qca_profile.greedy_calls++;
        qca_profile.greedy_seconds += now_seconds() - started;
        return INT_MAX / 4;
    }
    memcpy(uncovered, state->uncovered, (size_t) problem->nwords_rows * sizeof(unsigned long long));
    memcpy(active, state->active, (size_t) problem->nc * sizeof(unsigned char));
    memcpy(chosen, state->chosen, (size_t) problem->nc * sizeof(unsigned char));
    int chosen_count = state->chosen_count;
    while (bitset_count(uncovered, problem->nwords_rows) > 0) {
        int best_col = -1;
        int best_gain = 0;
        for (int col = 0; col < problem->nc; col++) {
            if (!active[col]) {
                continue;
            }
            int gain = 0;
            const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
            for (int w = 0; w < problem->nwords_rows; w++) {
                gain += popcount_ull(uncovered[w] & mask[w]);
            }
            if (gain > best_gain) {
                best_gain = gain;
                best_col = col;
            }
        }
        if (best_col < 0 || best_gain == 0) {
            R_Free(uncovered);
            R_Free(active);
            R_Free(chosen);
            qca_profile.greedy_calls++;
            qca_profile.greedy_seconds += now_seconds() - started;
            return INT_MAX / 4;
        }
        chosen[best_col] = 1;
        active[best_col] = 0;
        chosen_count++;
        for (int w = 0; w < problem->nwords_rows; w++) {
            uncovered[w] &= ~problem->col_masks[best_col * problem->nwords_rows + w];
        }
    }
    for (int col = 0; col < problem->nc; col++) {
        candidate_solution[col] = chosen[col];
    }
    R_Free(uncovered);
    R_Free(active);
    R_Free(chosen);
    qca_profile.greedy_calls++;
    qca_profile.greedy_seconds += now_seconds() - started;
    return chosen_count;
}
static int solution_still_covers_without(
    const qca_scp_problem *problem,
    const int *solution,
    int drop_col
) {
    for (int row = 0; row < problem->nr; row++) {
        int covered = 0;
        const unsigned long long drop_mask = row_mask(row);
        const int drop_word = row_word(row);
        for (int col = 0; col < problem->nc; col++) {
            if (!solution[col] || col == drop_col) {
                continue;
            }
            const unsigned long long *mask = problem->col_masks + col * problem->nwords_rows;
            if ((mask[drop_word] & drop_mask) != 0ULL) {
                covered = 1;
                break;
            }
        }
        if (!covered) {
            return 0;
        }
    }
    return 1;
}
static int cleanup_solution(
    const qca_scp_problem *problem,
    int *solution
) {
    int changed = 1;
    while (changed) {
        changed = 0;
        for (int col = 0; col < problem->nc; col++) {
            if (!solution[col]) {
                continue;
            }
            if (solution_still_covers_without(problem, solution, col)) {
                solution[col] = 0;
                changed = 1;
            }
        }
    }
    int size = 0;
    for (int col = 0; col < problem->nc; col++) {
        size += solution[col] != 0;
    }
    return size;
}
static void search_exact(
    qca_scp_search *search,
    qca_scp_state *state,
    int depth
) {
    const qca_scp_problem *problem = search->problem;
    qca_scp_workspace *workspace = search->workspaces + depth;
    if (search->found_target) {
        return;
    }
    qca_profile.nodes++;
    if (!apply_reductions(problem, state, search)) {
        qca_profile.leaves++;
        return;
    }
    if (state->chosen_count >= search->best_size) {
        qca_profile.leaves++;
        return;
    }
    if (search->best_size > state->chosen_count + 1) {
        int *candidate_solution = workspace->candidate_solution;
        if (candidate_solution != NULL) {
            int greedy_size = greedy_cover(problem, state, candidate_solution);
            if (greedy_size < INT_MAX / 4) {
                greedy_size = cleanup_solution(problem, candidate_solution);
            }
            if (greedy_size < search->best_size) {
                search->best_size = greedy_size;
                memcpy(search->best_solution, candidate_solution, (size_t) problem->nc * sizeof(int));
                if (search->proof_target_size >= 0 && greedy_size <= search->proof_target_size) {
                    search->found_target = 1;
                    qca_profile.leaves++;
                    return;
                }
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
    int lb = lower_bound(problem, state);
    if (search->proof_target_size >= 0 &&
        state->chosen_count + lb > search->proof_target_size) {
        qca_profile.leaves++;
        return;
    }
    if (state->chosen_count + lb >= search->best_size) {
        qca_profile.leaves++;
        return;
    }
    int row = choose_branch_row(problem, state);
    if (row < 0) {
        return;
    }
    int *branch_cols = workspace->branch_cols;
    int *gains = workspace->gains;
    int *scores = workspace->candidate_solution;
    if (branch_cols == NULL || gains == NULL || scores == NULL) {
        return;
    }
    int nbranch = 0;
    for (int i = problem->row_starts[row]; i < problem->row_starts[row + 1]; i++) {
        int col = problem->row_cols[i];
        if (state->active[col]) {
            branch_cols[nbranch++] = col;
            gains[col] = new_rows_covered(problem, state, col);
            scores[col] = branch_col_score(problem, state, col);
        }
    }
    {
        double started = now_seconds();
        sort_branch_cols(branch_cols, nbranch, gains, scores);
        qca_profile.branching_seconds += now_seconds() - started;
    }
    for (int i = 0; i < nbranch; i++) {
        int col = branch_cols[i];
        qca_scp_state *child = &(search->workspaces[depth + 1].state);
        copy_state(problem, state, child);
        cover_with_column(problem, child, col);
        search_exact(search, child, depth + 1);
        if (search->found_target) {
            return;
        }
    }
}
int qca_scp_solve_exact(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size
) {
    return qca_scp_solve_exact_with_incumbent(problem, solution, solution_size, NULL, 0, -1);
}
int qca_scp_solve_exact_with_incumbent(
    const qca_scp_problem *problem,
    int *solution,
    int *solution_size,
    const int *initial_solution,
    int initial_size,
    int proof_target_size
) {
    double started = now_seconds();
    qca_scp_state state;
    qca_scp_search search;
    qca_scp_workspace *workspaces = NULL;
    const int max_depth = problem->nc + 1;
    state.uncovered = (unsigned long long *) R_Calloc((size_t) problem->nwords_rows, unsigned long long);
    state.active = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
    state.chosen = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
    workspaces = (qca_scp_workspace *) R_Calloc((size_t) max_depth, qca_scp_workspace);
    if (state.uncovered == NULL || state.active == NULL || state.chosen == NULL || workspaces == NULL) {
        if (state.uncovered) R_Free(state.uncovered);
        if (state.active) R_Free(state.active);
        if (state.chosen) R_Free(state.chosen);
        if (workspaces) R_Free(workspaces);
        return 0;
    }
    for (int depth = 0; depth < max_depth; depth++) {
        workspaces[depth].state.uncovered = (unsigned long long *) R_Calloc((size_t) problem->nwords_rows, unsigned long long);
        workspaces[depth].state.active = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
        workspaces[depth].state.chosen = (unsigned char *) R_Calloc((size_t) problem->nc, unsigned char);
        workspaces[depth].branch_cols = (int *) R_Calloc((size_t) problem->nc, int);
        workspaces[depth].gains = (int *) R_Calloc((size_t) problem->nc, int);
        workspaces[depth].candidate_solution = (int *) R_Calloc((size_t) problem->nc, int);
        if (workspaces[depth].state.uncovered == NULL ||
            workspaces[depth].state.active == NULL ||
            workspaces[depth].state.chosen == NULL ||
            workspaces[depth].branch_cols == NULL ||
            workspaces[depth].gains == NULL ||
            workspaces[depth].candidate_solution == NULL) {
            for (int i = 0; i <= depth; i++) {
                if (workspaces[i].state.uncovered) R_Free(workspaces[i].state.uncovered);
                if (workspaces[i].state.active) R_Free(workspaces[i].state.active);
                if (workspaces[i].state.chosen) R_Free(workspaces[i].state.chosen);
                if (workspaces[i].branch_cols) R_Free(workspaces[i].branch_cols);
                if (workspaces[i].gains) R_Free(workspaces[i].gains);
                if (workspaces[i].candidate_solution) R_Free(workspaces[i].candidate_solution);
            }
            R_Free(state.uncovered);
            R_Free(state.active);
            R_Free(state.chosen);
            R_Free(workspaces);
            return 0;
        }
    }
    for (int row = 0; row < problem->nr; row++) {
        state.uncovered[row_word(row)] |= row_mask(row);
    }
    memset(state.active, 1, (size_t) problem->nc * sizeof(unsigned char));
    memset(state.chosen, 0, (size_t) problem->nc * sizeof(unsigned char));
    state.chosen_count = 0;
    search.problem = problem;
    search.best_solution = solution;
    search.best_size = problem->nc + 1;
    search.proof_target_size = proof_target_size;
    search.found_target = 0;
    search.workspaces = workspaces;
    memset(solution, 0, (size_t) problem->nc * sizeof(int));
    if (initial_solution != NULL &&
        initial_size > 0 &&
        initial_size <= problem->nc) {
        memcpy(solution, initial_solution, (size_t) problem->nc * sizeof(int));
        search.best_size = initial_size;
    }
    search_exact(&search, &state, 0);
    R_Free(state.uncovered);
    R_Free(state.active);
    R_Free(state.chosen);
    for (int depth = 0; depth < max_depth; depth++) {
        R_Free(workspaces[depth].state.uncovered);
        R_Free(workspaces[depth].state.active);
        R_Free(workspaces[depth].state.chosen);
        R_Free(workspaces[depth].branch_cols);
        R_Free(workspaces[depth].gains);
        R_Free(workspaces[depth].candidate_solution);
    }
    R_Free(workspaces);
    if (search.best_size > problem->nc) {
        qca_profile.total_seconds += now_seconds() - started;
        return 0;
    }
    *solution_size = search.best_size;
    qca_profile.total_seconds += now_seconds() - started;
    return 1;
}
void qca_scp_profile_reset(void) {
    memset(&qca_profile, 0, sizeof(qca_profile));
}
qca_scp_profile qca_scp_profile_get(void) {
    return qca_profile;
}
