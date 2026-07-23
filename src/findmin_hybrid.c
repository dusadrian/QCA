#include "qca_rinternals.h"
#include <R_ext/RS.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include "findmin_hybrid.h"
#include "findmin_lpsolve.h"
#include "pichart_presolve.h"
#include "scp_solver/scp_solver.h"
#include "findmin_lagrangian.h"

/*
All limits are node counts: node budgets are deterministic, so the same
chart always takes the same path and returns the same cover, regardless
of machine load. Worst-case time is controlled instead by scaling the
hard budgets with the core width (per-node cost grows with the number of
active columns), which is equally deterministic.
*/
#define HYBRID_SCP_PROBE_NODE_LIMIT 5000ULL
#define HYBRID_SCP_WIDE_PROBE_NODE_LIMIT 10000ULL
#define HYBRID_SCP_GAP_ONE_NODE_LIMIT 1500000ULL
#define HYBRID_SCP_GAP_TWO_NODE_LIMIT 600000ULL
#define HYBRID_SCP_LARGER_GAP_NODE_LIMIT 250000ULL
#define HYBRID_SCP_NODE_SCALE_COLUMNS 2000
#define HYBRID_SCP_VERY_WIDE_RATIO 10
#define DENSE_MASK_ROWS 20

typedef struct {
    int original_columns;
    int presolved_columns;
    int core_columns;
    int scp_attempted;
    int scp_limited;
    int lpsolve_fallback;
    int incumbent_requested;
    int incumbent_accepted;
} qca_hybrid_profile;

static qca_hybrid_profile hybrid_profile = {0};

typedef struct {
    int column;
    int coverage_count;
    uint64_t hash;
} qca_mask_column;

static int compare_mask_columns(const void *left, const void *right) {
    const qca_mask_column *a = (const qca_mask_column *)left;
    const qca_mask_column *b = (const qca_mask_column *)right;
    if (a->coverage_count != b->coverage_count) {
        return a->coverage_count < b->coverage_count ? 1 : -1;
    }
    return a->column < b->column ? -1 : a->column > b->column;
}

static Rboolean masks_equal(
    const uint64_t *left,
    const uint64_t *right,
    int nwords
) {
    for (int w = 0; w < nwords; ++w) {
        if (left[w] != right[w]) return FALSE;
    }
    return TRUE;
}

static Rboolean mask_is_subset(
    const uint64_t *candidate,
    const uint64_t *possible_superset,
    int nwords
) {
    for (int w = 0; w < nwords; ++w) {
        if (candidate[w] & ~possible_superset[w]) return FALSE;
    }
    return TRUE;
}

static int reduce_active_generic(
    const int *chart,
    int nrows,
    int ncols,
    unsigned char *active,
    int active_count
) {
    int result = active_count;
    int nwords_rows = (nrows + 63) / 64;
    int nwords_columns = (ncols + 63) / 64;
    qca_mask_column *columns = NULL;
    uint64_t *masks = NULL;
    uint64_t *row_index = NULL;
    int *row_retained = NULL;
    int *hash_table = NULL;
    size_t *touched_slots = NULL;
    size_t hash_size = 1;

    if (
        (size_t)ncols > SIZE_MAX / (size_t)nwords_rows / sizeof(uint64_t) ||
        (size_t)nrows > SIZE_MAX / (size_t)nwords_columns / sizeof(uint64_t)
    ) {
        return active_count;
    }

    columns = (qca_mask_column *)malloc(
        (size_t)active_count * sizeof(qca_mask_column)
    );
    masks = (uint64_t *)calloc(
        (size_t)ncols * (size_t)nwords_rows, sizeof(uint64_t)
    );
    row_index = (uint64_t *)calloc(
        (size_t)nrows * (size_t)nwords_columns, sizeof(uint64_t)
    );
    row_retained = (int *)calloc((size_t)nrows, sizeof(int));
    if (!columns || !masks || !row_index || !row_retained) goto cleanup;

    for (int c = 0, cc = 0; c < ncols; ++c) {
        if (!active[c]) continue;
        uint64_t hash = 1469598103934665603ULL;
        int count = 0;
        uint64_t *mask = &masks[(size_t)c * nwords_rows];
        for (int r = 0; r < nrows; ++r) {
            if (!chart[(size_t)c * nrows + r]) continue;
            mask[r >> 6] |= 1ULL << (r & 63);
            ++count;
        }
        for (int w = 0; w < nwords_rows; ++w) {
            hash ^= mask[w];
            hash *= 1099511628211ULL;
        }
        columns[cc++] = (qca_mask_column){
            .column = c,
            .coverage_count = count,
            .hash = hash
        };
    }
    qsort(
        columns, (size_t)active_count,
        sizeof(qca_mask_column), compare_mask_columns
    );

    while (hash_size < (size_t)active_count * 2U) {
        if (hash_size > SIZE_MAX / 2U) goto cleanup;
        hash_size <<= 1;
    }
    hash_table = (int *)malloc(hash_size * sizeof(int));
    touched_slots = (size_t *)malloc((size_t)active_count * sizeof(size_t));
    if (!hash_table || !touched_slots) goto cleanup;
    for (size_t slot = 0; slot < hash_size; ++slot) hash_table[slot] = -1;

    for (int group_start = 0; group_start < active_count;) {
        int group_end = group_start + 1;
        int coverage_count = columns[group_start].coverage_count;
        while (
            group_end < active_count &&
            columns[group_end].coverage_count == coverage_count
        ) {
            ++group_end;
        }

        /* Equal-cardinality coverage can dominate only when it is identical. */
        size_t touched_count = 0;
        for (int i = group_start; i < group_end; ++i) {
            int c = columns[i].column;
            size_t slot = (size_t)columns[i].hash & (hash_size - 1U);
            while (hash_table[slot] >= 0) {
                int previous_i = hash_table[slot];
                int previous = columns[previous_i].column;
                if (
                    columns[previous_i].hash == columns[i].hash &&
                    masks_equal(
                        &masks[(size_t)c * nwords_rows],
                        &masks[(size_t)previous * nwords_rows],
                        nwords_rows
                    )
                ) {
                    active[c] = 0;
                    --result;
                    break;
                }
                slot = (slot + 1U) & (hash_size - 1U);
            }
            if (active[c]) {
                hash_table[slot] = i;
                touched_slots[touched_count++] = slot;
            }
        }
        for (size_t i = 0; i < touched_count; ++i) {
            hash_table[touched_slots[i]] = -1;
        }

        /* Only strictly larger coverage groups are in row_index. Start from
           the candidate's rarest covered row, then verify exact containment. */
        for (int i = group_start; i < group_end; ++i) {
            int c = columns[i].column;
            if (!active[c]) continue;
            const uint64_t *candidate = &masks[(size_t)c * nwords_rows];
            int rarest_row = -1;
            int rarest_count = INT_MAX;
            for (int r = 0; r < nrows; ++r) {
                if (
                    (candidate[r >> 6] & (1ULL << (r & 63))) &&
                    row_retained[r] < rarest_count
                ) {
                    rarest_row = r;
                    rarest_count = row_retained[r];
                }
            }

            Rboolean dominated = FALSE;
            if (rarest_row < 0) {
                dominated = group_start > 0;
            }
            else if (rarest_count > 0) {
                const uint64_t *possible = &row_index[
                    (size_t)rarest_row * nwords_columns
                ];
                for (int w = 0; w < nwords_columns && !dominated; ++w) {
                    uint64_t bits = possible[w];
                    while (bits && !dominated) {
                        int bit = __builtin_ctzll(bits);
                        int previous = (w << 6) + bit;
                        if (
                            previous < ncols && active[previous] &&
                            mask_is_subset(
                                candidate,
                                &masks[(size_t)previous * nwords_rows],
                                nwords_rows
                            )
                        ) {
                            dominated = TRUE;
                        }
                        bits &= bits - 1ULL;
                    }
                }
            }
            if (dominated) {
                active[c] = 0;
                --result;
            }
        }

        /* Insert the surviving group only after every same-size column was
           queried, so the index always contains strict supersets only. */
        for (int i = group_start; i < group_end; ++i) {
            int c = columns[i].column;
            if (!active[c]) continue;
            const uint64_t *mask = &masks[(size_t)c * nwords_rows];
            for (int r = 0; r < nrows; ++r) {
                if (!(mask[r >> 6] & (1ULL << (r & 63)))) continue;
                row_index[(size_t)r * nwords_columns + (c >> 6)] |=
                    1ULL << (c & 63);
                ++row_retained[r];
            }
        }
        group_start = group_end;
    }

cleanup:
    free(columns);
    free(masks);
    free(row_index);
    free(row_retained);
    free(hash_table);
    free(touched_slots);
    return result;
}

int qca_reduce_active_columns(
    const int *chart,
    int nrows,
    int ncols,
    unsigned char *active
) {
    int active_count = 0;
    if (chart == NULL || active == NULL || nrows <= 0 || ncols <= 0) {
        return -1;
    }

    for (int c = 0; c < ncols; ++c) {
        active_count += active[c] != 0;
    }
    if (active_count <= 1) {
        return active_count;
    }
    if (nrows > DENSE_MASK_ROWS) {
        return reduce_active_generic(
            chart, nrows, ncols, active, active_count
        );
    }

    size_t mask_count = (size_t)1U << nrows;
    int *representative = (int *)malloc(mask_count * sizeof(int));
    unsigned char *has_superset = (unsigned char *)calloc(mask_count, sizeof(unsigned char));
    if (representative == NULL || has_superset == NULL) {
        free(representative);
        free(has_superset);
        return -1;
    }
    for (size_t mask = 0; mask < mask_count; ++mask) {
        representative[mask] = -1;
    }

    for (int c = 0; c < ncols; ++c) {
        if (!active[c]) continue;
        size_t coverage = 0;
        for (int r = 0; r < nrows; ++r) {
            if (chart[(size_t)c * nrows + r]) {
                coverage |= (size_t)1U << r;
            }
        }
        if (representative[coverage] < 0) {
            representative[coverage] = c;
            has_superset[coverage] = 1;
        }
    }

    /* Superset zeta transform. Afterwards has_superset[m] says whether an
       active coverage mask containing m exists. */
    for (int bit = 0; bit < nrows; ++bit) {
        size_t bit_mask = (size_t)1U << bit;
        for (size_t mask = 0; mask < mask_count; ++mask) {
            if (!(mask & bit_mask) && has_superset[mask | bit_mask]) {
                has_superset[mask] = 1;
            }
        }
    }

    memset(active, 0, (size_t)ncols * sizeof(unsigned char));
    active_count = 0;
    for (size_t mask = 0; mask < mask_count; ++mask) {
        int column = representative[mask];
        if (column < 0) continue;

        int dominated = 0;
        for (int bit = 0; bit < nrows && !dominated; ++bit) {
            size_t bit_mask = (size_t)1U << bit;
            if (!(mask & bit_mask) && has_superset[mask | bit_mask]) {
                dominated = 1;
            }
        }
        if (!dominated) {
            active[column] = 1;
            ++active_count;
        }
    }

    free(representative);
    free(has_superset);
    return active_count;
}

static int solve_scp_from_int_matrix(
    const int *p_chart,
    const int nr,
    const int nc,
    const int *initial_indices,
    int initial_solmin,
    int *solution
) {
    int ok = 0;
    int solution_size = 0;
    int incumbent_size = -1;
    double lagr_lb = -1e308;
    int core_nc = 0;
    int nnz = 0;
    int nwords_rows = (nr + 63) / 64;

    int *incumbent_indices = NULL;
    int *core_map = NULL;
    int *reverse_map = NULL;
    int *core_incumbent = NULL;
    int *core_solution = NULL;
    int *core_lp_indices = NULL;
    int *row_counts = NULL;
    int *row_starts = NULL;
    int *row_cols = NULL;
    int *core_chart = NULL;
    unsigned char *improving_core = NULL;
    unsigned char *include_column = NULL;
    unsigned long long *col_masks = NULL;
    double *lagr_scores = NULL;
    double *core_priority = NULL;

    if (p_chart == NULL || solution == NULL || nr <= 0 || nc <= 0) {
        return 0;
    }
    memset(solution, 0, (size_t)nc * sizeof(int));
    memset(&hybrid_profile, 0, sizeof(hybrid_profile));
    hybrid_profile.original_columns = nc;
    hybrid_profile.presolved_columns = nc;
    hybrid_profile.incumbent_requested =
        initial_indices != NULL && initial_solmin > 0;
    hybrid_profile.incumbent_accepted =
        hybrid_profile.incumbent_requested;

    incumbent_indices = (int *)calloc((size_t)nc, sizeof(int));
    lagr_scores = (double *)calloc((size_t)nc, sizeof(double));
    improving_core = (unsigned char *)calloc((size_t)nc, sizeof(unsigned char));
    include_column = (unsigned char *)calloc((size_t)nc, sizeof(unsigned char));
    reverse_map = (int *)malloc((size_t)nc * sizeof(int));
    if (!incumbent_indices || !lagr_scores || !improving_core || !include_column || !reverse_map) {
        goto cleanup;
    }
    for (int c = 0; c < nc; ++c) reverse_map[c] = -1;

    solvePIchart_lagrangian_prepare_with_incumbent(
        (int *)p_chart,
        nc,
        nr,
        NULL,
        initial_indices,
        initial_solmin,
        incumbent_indices,
        &incumbent_size,
        &lagr_lb,
        lagr_scores,
        improving_core,
        NULL
    );

    if (incumbent_size <= 0 || incumbent_size > nc) {
        goto cleanup;
    }

    for (int i = 0; i < incumbent_size; ++i) {
        int col = incumbent_indices[i];
        if (col < 0 || col >= nc) goto cleanup;
        include_column[col] = 1;
    }

    if (
        lagr_lb > -1e307 &&
        (double)incumbent_size <= ceil(lagr_lb - 1e-12) + 1e-12
    ) {
        for (int i = 0; i < incumbent_size; ++i) {
            solution[incumbent_indices[i]] = 1;
        }
        ok = 1;
        goto cleanup;
    }

    for (int c = 0; c < nc; ++c) {
        if (improving_core[c]) include_column[c] = 1;
        if (include_column[c]) ++core_nc;
    }
    hybrid_profile.core_columns = core_nc;
    if (core_nc <= 0) goto cleanup;

    core_map = (int *)malloc((size_t)core_nc * sizeof(int));
    core_incumbent = (int *)calloc((size_t)core_nc, sizeof(int));
    core_solution = (int *)calloc((size_t)core_nc, sizeof(int));
    core_priority = (double *)calloc((size_t)core_nc, sizeof(double));
    core_chart = (int *)calloc((size_t)nr * (size_t)core_nc, sizeof(int));
    if (!core_map || !core_incumbent || !core_solution || !core_priority || !core_chart) {
        goto cleanup;
    }

    for (int c = 0, cc = 0; c < nc; ++c) {
        if (!include_column[c]) continue;
        core_map[cc] = c;
        reverse_map[c] = cc;
        core_priority[cc] = lagr_scores[c];
        for (int r = 0; r < nr; ++r) {
            int covered = p_chart[c * nr + r] != 0;
            core_chart[cc * nr + r] = covered;
            nnz += covered;
        }
        ++cc;
    }
    for (int i = 0; i < incumbent_size; ++i) {
        int cc = reverse_map[incumbent_indices[i]];
        if (cc < 0) goto cleanup;
        core_incumbent[cc] = 1;
    }

    row_counts = (int *)calloc((size_t)nr, sizeof(int));
    row_starts = (int *)calloc((size_t)nr + 1, sizeof(int));
    row_cols = (int *)calloc((size_t)(nnz > 0 ? nnz : 1), sizeof(int));
    col_masks = (unsigned long long *)calloc(
        (size_t)core_nc * (size_t)nwords_rows,
        sizeof(unsigned long long)
    );
    if (!row_counts || !row_starts || !row_cols || !col_masks) goto cleanup;

    for (int c = 0; c < core_nc; ++c) {
        for (int r = 0; r < nr; ++r) {
            if (!core_chart[c * nr + r]) continue;
            row_counts[r]++;
            col_masks[(size_t)c * nwords_rows + (r >> 6)] |= 1ULL << (r & 63);
        }
    }
    row_starts[0] = 0;
    for (int r = 0; r < nr; ++r) {
        row_starts[r + 1] = row_starts[r] + row_counts[r];
    }
    memset(row_counts, 0, (size_t)nr * sizeof(int));
    for (int c = 0; c < core_nc; ++c) {
        for (int r = 0; r < nr; ++r) {
            if (!core_chart[c * nr + r]) continue;
            row_cols[row_starts[r] + row_counts[r]++] = c;
        }
    }

    /*
    Probe the internal finisher on every uncertified core: small cores
    usually resolve within the probe budget (cheaper than an lp_solve
    setup), wide cores extend geometrically only while making proof
    progress. All budgets are node counts (deterministic); the hard caps
    shrink with the core width so the worst case stays time-bounded.
    */
    {
        qca_scp_problem problem = {
            .nr = nr,
            .nc = core_nc,
            .nwords_rows = nwords_rows,
            .row_starts = row_starts,
            .row_cols = row_cols,
            .col_masks = col_masks,
            .branch_priority = core_priority
        };
        int very_wide = nc >= nr * HYBRID_SCP_VERY_WIDE_RATIO;
        int lower_bound_size = lagr_lb > -1e307
            ? (int)ceil(lagr_lb - 1e-12)
            : 0;
        int optimality_gap = incumbent_size - lower_bound_size;
        unsigned long long hard_node_limit = optimality_gap <= 1
            ? HYBRID_SCP_GAP_ONE_NODE_LIMIT
            : optimality_gap == 2
                ? HYBRID_SCP_GAP_TWO_NODE_LIMIT
                : HYBRID_SCP_LARGER_GAP_NODE_LIMIT;
        unsigned long long node_scale = core_nc > HYBRID_SCP_NODE_SCALE_COLUMNS
            ? (unsigned long long)core_nc / HYBRID_SCP_NODE_SCALE_COLUMNS
            : 1ULL;
        hard_node_limit /= node_scale;
        if (hard_node_limit < HYBRID_SCP_PROBE_NODE_LIMIT) {
            hard_node_limit = HYBRID_SCP_PROBE_NODE_LIMIT;
        }
        hybrid_profile.scp_attempted = 1;
        qca_scp_result result = qca_scp_solve_exact_with_incumbent_adaptive(
            &problem,
            core_solution,
            &solution_size,
            core_incumbent,
            incumbent_size,
            -1,
            very_wide
                ? HYBRID_SCP_WIDE_PROBE_NODE_LIMIT
                : HYBRID_SCP_PROBE_NODE_LIMIT,
            0.0, /* node budgets only: deterministic */
            hard_node_limit,
            0.0
        );
        if (result == QCA_SCP_SOLUTION) {
            for (int c = 0; c < core_nc; ++c) {
                if (core_solution[c]) solution[core_map[c]] = 1;
            }
            ok = 1;
            goto cleanup;
        }
        hybrid_profile.scp_limited = result == QCA_SCP_LIMIT;
    }

    hybrid_profile.lpsolve_fallback = 1;
    core_lp_indices = (int *)calloc((size_t)core_nc, sizeof(int));
    if (!core_lp_indices) goto cleanup;
    if (!solvePIchart_lpsolve(core_chart, nr, core_nc, core_lp_indices, &solution_size)) {
        goto cleanup;
    }
    for (int i = 0; i < solution_size; ++i) {
        int cc = core_lp_indices[i];
        if (cc < 0 || cc >= core_nc) goto cleanup;
        solution[core_map[cc]] = 1;
    }
    ok = 1;

cleanup:
    free(incumbent_indices);
    free(core_map);
    free(reverse_map);
    free(core_incumbent);
    free(core_solution);
    free(core_lp_indices);
    free(row_counts);
    free(row_starts);
    free(row_cols);
    free(core_chart);
    free(improving_core);
    free(include_column);
    free(col_masks);
    free(lagr_scores);
    free(core_priority);
    return ok;
}

static Rboolean solvePIchart_hybrid_active_impl(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    const int *initial_indices,
    int initial_solmin,
    int *indices,
    int *solmin
) {
    int compact_ncols = 0;
    for (int c = 0; c < ncols; ++c) {
        compact_ncols += active == NULL || active[c] != 0;
    }
    if (compact_ncols <= 0) return FALSE;

    int *solution = (int *)calloc((size_t)compact_ncols, sizeof(int));
    int *compact_map = NULL;
    int *reverse_map = NULL;
    int *compact_initial = NULL;
    int *compact_chart = NULL;
    const int *solver_chart = chart;
    if (solution == NULL) {
        return FALSE;
    }

    if (compact_ncols != ncols) {
        compact_map = (int *)malloc((size_t)compact_ncols * sizeof(int));
        reverse_map = (int *)malloc((size_t)ncols * sizeof(int));
        compact_chart = (int *)malloc(
            (size_t)nrows * (size_t)compact_ncols * sizeof(int)
        );
        if (compact_map == NULL || reverse_map == NULL || compact_chart == NULL) {
            free(solution);
            free(compact_map);
            free(reverse_map);
            free(compact_chart);
            return FALSE;
        }
        for (int c = 0; c < ncols; ++c) reverse_map[c] = -1;
        for (int c = 0, cc = 0; c < ncols; ++c) {
            if (active != NULL && !active[c]) continue;
            compact_map[cc] = c;
            reverse_map[c] = cc;
            memcpy(
                &compact_chart[(size_t)cc * nrows],
                &chart[(size_t)c * nrows],
                (size_t)nrows * sizeof(int)
            );
            ++cc;
        }
        solver_chart = compact_chart;
    }

    if (initial_indices && initial_solmin > 0) {
        compact_initial = (int *)malloc((size_t)initial_solmin * sizeof(int));
        if (!compact_initial) {
            free(solution);
            free(compact_map);
            free(reverse_map);
            free(compact_chart);
            return FALSE;
        }
        for (int i = 0; i < initial_solmin; ++i) {
            int col = initial_indices[i];
            int compact_col = compact_map == NULL
                ? col
                : (col >= 0 && col < ncols ? reverse_map[col] : -1);
            if (compact_col < 0 || compact_col >= compact_ncols) {
                free(solution);
                free(compact_map);
                free(reverse_map);
                free(compact_chart);
                free(compact_initial);
                return FALSE;
            }
            compact_initial[i] = compact_col;
        }
    }

    if (!solve_scp_from_int_matrix(
        solver_chart, nrows, compact_ncols,
        compact_initial, initial_solmin, solution
    )) {
        free(solution);
        free(compact_map);
        free(reverse_map);
        free(compact_chart);
        free(compact_initial);
        return FALSE;
    }
    hybrid_profile.original_columns = ncols;
    hybrid_profile.presolved_columns = compact_ncols;

    *solmin = 0;
    for (int c = 0; c < compact_ncols; c++) {
        if (solution[c] > 0) {
            indices[*solmin] = compact_map == NULL ? c : compact_map[c];
            (*solmin)++;
        }
    }

    free(solution);
    free(compact_map);
    free(reverse_map);
    free(compact_chart);
    free(compact_initial);
    return TRUE;
}

Rboolean solvePIchart_hybrid_active(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    int *indices,
    int *solmin
) {
    return solvePIchart_hybrid_active_impl(
        chart, nrows, ncols, active, NULL, 0, indices, solmin
    );
}

Rboolean solvePIchart_hybrid_active_with_incumbent(
    const int *chart,
    int nrows,
    int ncols,
    const unsigned char *active,
    const int *initial_indices,
    int initial_solmin,
    int *indices,
    int *solmin
) {
    return solvePIchart_hybrid_active_impl(
        chart, nrows, ncols, active,
        initial_indices, initial_solmin, indices, solmin
    );
}

Rboolean solvePIchart_hybrid(
    const int *chart,
    int nrows,
    int ncols,
    int *indices,
    int *solmin
) {
    unsigned char *active = (unsigned char *)malloc((size_t)ncols);
    if (active == NULL) return FALSE;
    memset(active, 1, (size_t)ncols);

    int active_count = qca_reduce_active_columns(chart, nrows, ncols, active);
    if (active_count <= 0) {
        free(active);
        return FALSE;
    }
    Rboolean ok = solvePIchart_hybrid_active_impl(
        chart, nrows, ncols, active, NULL, 0, indices, solmin
    );
    free(active);
    return ok;
}

SEXP C_getScpProfile(void) {
    qca_scp_profile profile = qca_scp_profile_get();
    SEXP out = PROTECT(allocVector(VECSXP, 21));
    SEXP names = PROTECT(allocVector(STRSXP, 21));

    SET_STRING_ELT(names, 0, mkChar("total_seconds"));
    SET_STRING_ELT(names, 1, mkChar("reductions_seconds"));
    SET_STRING_ELT(names, 2, mkChar("lower_bound_seconds"));
    SET_STRING_ELT(names, 3, mkChar("greedy_seconds"));
    SET_STRING_ELT(names, 4, mkChar("branching_seconds"));
    SET_STRING_ELT(names, 5, mkChar("nodes"));
    SET_STRING_ELT(names, 6, mkChar("leaves"));
    SET_STRING_ELT(names, 7, mkChar("reduction_calls"));
    SET_STRING_ELT(names, 8, mkChar("lower_bound_calls"));
    SET_STRING_ELT(names, 9, mkChar("greedy_calls"));
    SET_STRING_ELT(names, 10, mkChar("original_columns"));
    SET_STRING_ELT(names, 11, mkChar("presolved_columns"));
    SET_STRING_ELT(names, 12, mkChar("core_columns"));
    SET_STRING_ELT(names, 13, mkChar("scp_attempted"));
    SET_STRING_ELT(names, 14, mkChar("scp_limited"));
    SET_STRING_ELT(names, 15, mkChar("lpsolve_fallback"));
    SET_STRING_ELT(names, 16, mkChar("root_branches_total"));
    SET_STRING_ELT(names, 17, mkChar("root_branches_completed"));
    SET_STRING_ELT(names, 18, mkChar("adaptive_extensions"));
    SET_STRING_ELT(names, 19, mkChar("incumbent_requested"));
    SET_STRING_ELT(names, 20, mkChar("incumbent_accepted"));

    SET_VECTOR_ELT(out, 0, ScalarReal(profile.total_seconds));
    SET_VECTOR_ELT(out, 1, ScalarReal(profile.reductions_seconds));
    SET_VECTOR_ELT(out, 2, ScalarReal(profile.lower_bound_seconds));
    SET_VECTOR_ELT(out, 3, ScalarReal(profile.greedy_seconds));
    SET_VECTOR_ELT(out, 4, ScalarReal(profile.branching_seconds));
    SET_VECTOR_ELT(out, 5, ScalarReal((double) profile.nodes));
    SET_VECTOR_ELT(out, 6, ScalarReal((double) profile.leaves));
    SET_VECTOR_ELT(out, 7, ScalarReal((double) profile.reduction_calls));
    SET_VECTOR_ELT(out, 8, ScalarReal((double) profile.lower_bound_calls));
    SET_VECTOR_ELT(out, 9, ScalarReal((double) profile.greedy_calls));
    SET_VECTOR_ELT(out, 10, ScalarInteger(hybrid_profile.original_columns));
    SET_VECTOR_ELT(out, 11, ScalarInteger(hybrid_profile.presolved_columns));
    SET_VECTOR_ELT(out, 12, ScalarInteger(hybrid_profile.core_columns));
    SET_VECTOR_ELT(out, 13, ScalarLogical(hybrid_profile.scp_attempted));
    SET_VECTOR_ELT(out, 14, ScalarLogical(hybrid_profile.scp_limited));
    SET_VECTOR_ELT(out, 15, ScalarLogical(hybrid_profile.lpsolve_fallback));
    SET_VECTOR_ELT(out, 16, ScalarReal((double) profile.root_branches_total));
    SET_VECTOR_ELT(out, 17, ScalarReal((double) profile.root_branches_completed));
    SET_VECTOR_ELT(out, 18, ScalarReal((double) profile.adaptive_extensions));
    SET_VECTOR_ELT(out, 19, ScalarLogical(hybrid_profile.incumbent_requested));
    SET_VECTOR_ELT(out, 20, ScalarLogical(hybrid_profile.incumbent_accepted));
    setAttrib(out, R_NamesSymbol, names);

    UNPROTECT(2);
    return out;
}

SEXP C_resetScpProfile(void) {
    qca_scp_profile_reset();
    memset(&hybrid_profile, 0, sizeof(hybrid_profile));
    return R_NilValue;
}

SEXP C_findminHybridInternal(SEXP chart) {
    if (!isMatrix(chart) || TYPEOF(chart) != LGLSXP) {
        error("C_findminHybridInternal expects a logical matrix.");
    }

    const int nr = nrows(chart);
    const int nc = ncols(chart);
    const int *p_chart = LOGICAL(chart);
    int *solution = (int *) R_Calloc((size_t) nc, int);
    if (solution == NULL) {
        error("Failed to allocate SCP solver solution vector.");
    }

    int *indices = (int *)R_Calloc((size_t)nc, int);
    int solmin = 0;
    if (indices == NULL || !solvePIchart_hybrid(p_chart, nr, nc, indices, &solmin)) {
        if (indices != NULL) R_Free(indices);
        R_Free(solution);
        error("Internal SCP solver failed to find a feasible exact cover.");
    }
    for (int i = 0; i < solmin; ++i) {
        solution[indices[i]] = 1;
    }
    R_Free(indices);

    SEXP out = PROTECT(allocVector(REALSXP, nc));
    double *p_out = REAL(out);
    for (int c = 0; c < nc; c++) {
        p_out[c] = (double) solution[c];
    }

    R_Free(solution);
    UNPROTECT(1);
    return out;
}
