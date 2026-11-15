/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "cover_encoding.h"
#include "native_sat_support.h"

typedef struct {
    NsatSolver *solver;
    const NsatVerifier *checker;
    NsatBudget budget;
    NsatClauseId emitted;
} Encoder;

int nsat_cover_selector(int inputs, int term, int input, int value)
{
    return 2 * (inputs * term + input) + value;
}

static bool emit(Encoder *e, const int *lits, size_t count)
{
    const int *expected;
    size_t width;
    if (e->emitted == UINT64_MAX ||
        !nsat_verify_clause(e->checker, ++e->emitted, &expected, &width) || width != count)
        return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_PROOF);
    for (size_t i = 0; i < count; ++i) {
        if (!nsat_work(&e->budget, 1)) return false;
        if (lits[i] != expected[i])
            return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_PROOF);
    }
    NsatStatus status = nsat_add_original(e->solver, lits, count);
    if (status != NSAT_OK)
        return nsat_stop(&e->budget, status, nsat_stats(e->solver).stop);
    return true;
}

static bool write_target(Encoder *e, const NsatCoverInput *q, int *scratch)
{
    const int inputs = (int)q->inputs, target = (int)q->target;
    const int literal_vars = 2 * inputs * target;
    const int used_base = literal_vars + (int)q->positives * target;

    for (size_t row = 0; row < q->positives; ++row) {
        size_t source = q->row_order ? q->row_order[row] : row;
        if (source >= q->positives)
            return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        const int *point = q->on + source * q->inputs;
        const int first = (int)row * target + 1;
        for (int t = 0; t < target; ++t) {
            if (!nsat_work(&e->budget, 1)) return false;
            scratch[t] = literal_vars + first + t;
        }
        if (!emit(e, scratch, q->target)) return false;
        for (int t = 0; t < target; ++t) {
            int assigned = literal_vars + first + t;
            int used = used_base + first + t;
            for (int i = 0; i < inputs; ++i) {
                if (point[i] != 1 && point[i] != 2)
                    return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_INPUT);
                int opposite = nsat_cover_selector(inputs, t, i, 3 - point[i]);
                if (!emit(e, (int[]){-assigned, -opposite}, 2)) return false;
            }
            if (!emit(e, (int[]){-assigned, used}, 2)) return false;
            if (!row) {
                if (!emit(e, (int[]){-used, assigned}, 2) ||
                    (t && !emit(e, (int[]){-assigned}, 1))) return false;
            } else {
                int previous = used - target;
                if (!emit(e, (int[]){-previous, used}, 2) ||
                    !emit(e, (int[]){-used, previous, assigned}, 3) ||
                    (t && !emit(e, (int[]){-assigned, previous - 1}, 2))) return false;
            }
        }
    }
    for (int t = 0; t < target; ++t) {
        for (int i = 0; i < inputs; ++i) {
            int zero = nsat_cover_selector(inputs, t, i, 1);
            if (!emit(e, (int[]){-zero, -(zero + 1)}, 2)) return false;
        }
        for (size_t row = 0; row < q->negatives; ++row) {
            const int *point = q->off + row * q->inputs;
            for (int i = 0; i < inputs; ++i) {
                if (!nsat_work(&e->budget, 1)) return false;
                if (point[i] != 1 && point[i] != 2)
                    return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_INPUT);
                scratch[i] = nsat_cover_selector(inputs, t, i, 3 - point[i]);
            }
            if (!emit(e, scratch, q->inputs)) return false;
        }
    }
    return true;
}

NsatStatus nsat_cover_load_solver(NsatSolver **out, const NsatCoverInput *q,
                                  const NsatVerifier *checker, NsatLimits limits,
                                  NsatStopReason *stop)
{
    *out = NULL;
    *stop = NSAT_STOP_INPUT;
    size_t max = (size_t)INT_MAX / 2, product;
    if (!checker || !q->on || !q->off || !q->inputs || !q->target ||
        !q->positives || !q->negatives || q->inputs > max ||
        q->positives > max - q->inputs ||
        q->target > max / (q->inputs + q->positives) ||
        (q->row_order ? q->order_count != q->positives : q->order_count != 0) ||
        !nsat_mul(q->positives, q->inputs, &product) || product > SIZE_MAX / sizeof(int) ||
        !nsat_mul(q->negatives, q->inputs, &product) || product > SIZE_MAX / sizeof(int))
        return NSAT_ERROR;
    int variables = (int)(2 * (q->inputs + q->positives) * q->target);
    if (variables != nsat_verify_variables(checker)) {
        *stop = NSAT_STOP_PROOF;
        return NSAT_ERROR;
    }
    NsatStatus status = nsat_init(out, variables, limits);
    if (status != NSAT_OK) {
        *stop = *out ? nsat_stats(*out).stop : NSAT_STOP_NONE;
        return status;
    }
    Encoder e = {.solver = *out, .checker = checker, .budget = {.limits = limits}};
    size_t width = q->inputs > q->target ? q->inputs : q->target;
    int *scratch = nsat_alloc(&e.budget, width, sizeof(*scratch));
    if (scratch && write_target(&e, q, scratch) &&
        e.emitted != nsat_verify_original_count(checker))
        nsat_stop(&e.budget, NSAT_ERROR, NSAT_STOP_PROOF);
    nsat_free(&e.budget, scratch, width * sizeof(*scratch));
    *stop = e.budget.stats.stop;
    return e.budget.status;
}
