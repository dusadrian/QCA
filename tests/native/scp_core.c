/* Independent R test harness for the actual native sources. No root hybrid
   preparation or external optimizer can discharge these test instances. */
#include <R.h>
#include <Rinternals.h>
#include "scp_solver/scp_solver.c"
#include "scp_solver/scp_relaxation.c"

SEXP test_scp_core_seeded(SEXP chart, SEXP method, SEXP target, SEXP weights, SEXP seed) {
    int nr = nrows(chart), nc = ncols(chart), words = (nr + 63) / 64;
    const int *a = LOGICAL(chart);
    int *starts = (int *)R_alloc(nr + 1, sizeof(int));
    int *cols = (int *)R_alloc((size_t)nr * nc, sizeof(int));
    unsigned long long *masks = (unsigned long long *)R_alloc((size_t)nc * words, sizeof(unsigned long long));
    int *initial = (int *)R_alloc(nc, sizeof(int));
    int *solution = (int *)R_alloc(nc, sizeof(int));
    memset(masks, 0, (size_t)nc * words * sizeof(unsigned long long));
    int count = 0;
    for (int r = 0; r < nr; ++r) {
        starts[r] = count;
        for (int c = 0; c < nc; ++c) if (a[c * nr + r]) {
            cols[count++] = c;
            masks[(size_t)c * words + (r >> 6)] |= 1ULL << (r & 63);
        }
    }
    starts[nr] = count;
    int initial_size = 0;
    for (int c = 0; c < nc; ++c) {
        initial[c] = seed == R_NilValue ? 1 : INTEGER(seed)[c];
        initial_size += initial[c];
    }
    qca_scp_problem problem = {
        .nr = nr, .nc = nc, .nwords_rows = words,
        .row_starts = starts, .row_cols = cols, .col_masks = masks,
        .target_propagation = asInteger(method) >= 1,
        .initial_row_dual = REAL(weights),
        .lagrangian_iterations = asInteger(method) >= 2 ? 8 : 0,
        .selective_bounds = asInteger(method) == 3
    };
    int size = 0;
    qca_scp_profile_reset();
    qca_scp_result result = qca_scp_solve_exact_with_incumbent(
        &problem, solution, &size, initial, initial_size, asInteger(target)
    );
    SEXP out = PROTECT(allocVector(VECSXP, 3));
    SEXP selected = PROTECT(allocVector(INTSXP, nc));
    for (int c = 0; c < nc; ++c) INTEGER(selected)[c] = solution[c];
    SET_VECTOR_ELT(out, 0, ScalarInteger(result));
    SET_VECTOR_ELT(out, 1, selected);
    SET_VECTOR_ELT(out, 2, ScalarReal((double)qca_scp_profile_get().lagrangian_calls));
    UNPROTECT(2);
    return out;
}

SEXP test_scp_core(SEXP chart, SEXP method, SEXP target, SEXP weights) {
    return test_scp_core_seeded(chart, method, target, weights, R_NilValue);
}

SEXP test_scp_statistics(void) {
    qca_scp_profile p = qca_scp_profile_get();
    SEXP out = PROTECT(allocVector(REALSXP, 4));
    REAL(out)[0] = (double)p.nodes;
    REAL(out)[1] = (double)p.lagrangian_iterations;
    REAL(out)[2] = (double)p.incumbent_improvements;
    REAL(out)[3] = (double)p.last_improvement_node;
    UNPROTECT(1);
    return out;
}
