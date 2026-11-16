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

/* A slot contains at most one equality per condition. An assignment to a
 * positive row forbids every differing equality. Each negative row must
 * differ on at least one selected equality in every slot. First-use order
 * only removes permutations of slots; an unused slot can duplicate a used
 * legal term, so satisfiability means a cover with AT MOST target terms. */
static bool write_multivalue(Encoder *e, const NsatMultivalueInput *q,
                             int values, int *scratch)
{
    const int k = (int)q->target;
    const int selectors = values * k;
    const int prefixes = selectors + (int)q->positives * k;
    for (size_t r = 0; r < q->positives; ++r) {
        int first = (int)r * k + 1;
        for (int t = 0; t < k; ++t) {
            if (!nsat_work(&e->budget, 1)) return false;
            scratch[t] = selectors + first + t;
        }
        if (!emit(e, scratch, q->target)) return false;
        for (int t = 0; t < k; ++t) {
            int a = selectors + first + t, u = prefixes + first + t;
            int base = t * values;
            for (size_t c = 0; c < q->inputs; ++c) {
                int point = q->on[r * q->inputs + c];
                if (point < 1 || point > q->levels[c])
                    return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_INPUT);
                for (int value = 1; value <= q->levels[c]; ++value) {
                    if (!nsat_work(&e->budget, 1)) return false;
                    if (value != point && !emit(e, (int[]){-a, -(base + value)}, 2))
                        return false;
                }
                base += q->levels[c];
            }
            if (!emit(e, (int[]){-a, u}, 2)) return false;
            if (!r) {
                if (!emit(e, (int[]){-u, a}, 2) ||
                    (t && !emit(e, (int[]){-a}, 1))) return false;
            } else {
                int previous = u - k;
                if (!emit(e, (int[]){-previous, u}, 2) ||
                    !emit(e, (int[]){-u, previous, a}, 3) ||
                    (t && !emit(e, (int[]){-a, previous - 1}, 2))) return false;
            }
        }
    }
    for (int t = 0; t < k; ++t) {
        int base = t * values;
        for (size_t c = 0; c < q->inputs; ++c) {
            for (int a = 1; a <= q->levels[c]; ++a) {
                if (!nsat_work(&e->budget, 1)) return false;
                for (int b = a + 1; b <= q->levels[c]; ++b)
                    if (!emit(e, (int[]){-(base + a), -(base + b)}, 2)) return false;
            }
            base += q->levels[c];
        }
        for (size_t r = 0; r < q->negatives; ++r) {
            size_t count = 0;
            base = t * values;
            for (size_t c = 0; c < q->inputs; ++c) {
                int point = q->off[r * q->inputs + c];
                if (point < 1 || point > q->levels[c])
                    return nsat_stop(&e->budget, NSAT_ERROR, NSAT_STOP_INPUT);
                for (int value = 1; value <= q->levels[c]; ++value) {
                    if (!nsat_work(&e->budget, 1)) return false;
                    if (value != point) scratch[count++] = base + value;
                }
                base += q->levels[c];
            }
            if (!emit(e, scratch, count)) return false;
        }
    }
    return true;
}

NsatStatus nsat_multivalue_load_solver(NsatSolver **out, const NsatMultivalueInput *q,
                                      const NsatVerifier *checker, NsatLimits limits,
                                      NsatStopReason *stop)
{
    if (!out || !stop) return NSAT_ERROR;
    *out = NULL;
    *stop = NSAT_STOP_INPUT;
    size_t cells;
    if (!q || !checker || !q->on || !q->off || !q->levels || !q->inputs ||
        !q->positives || !q->negatives || !q->target || q->inputs > INT_MAX ||
        q->positives > (size_t)INT_MAX / 2 || q->target > INT_MAX ||
        !nsat_mul(q->inputs, q->positives, &cells) || cells > SIZE_MAX / sizeof(int) ||
        !nsat_mul(q->inputs, q->negatives, &cells) || cells > SIZE_MAX / sizeof(int))
        return NSAT_ERROR;
    Encoder e = {.checker = checker, .budget = {.limits = limits}};
    size_t values = 0, width = 0;
    int *scratch = NULL;
    for (size_t c = 0; c < q->inputs; ++c) {
        if (!nsat_work(&e.budget, 1)) goto done;
        if (q->levels[c] < 1 || (size_t)q->levels[c] > (size_t)INT_MAX - values)
            goto invalid;
        values += (size_t)q->levels[c];
    }
    if (values > (size_t)INT_MAX - 2 * q->positives ||
        q->target > (size_t)INT_MAX / (values + 2 * q->positives)) goto invalid;
    int variables = (int)((values + 2 * q->positives) * q->target);
    if (variables != nsat_verify_variables(checker)) {
        nsat_stop(&e.budget, NSAT_ERROR, NSAT_STOP_PROOF);
        goto done;
    }
    NsatStatus status = nsat_init(out, variables, limits);
    if (status != NSAT_OK) {
        nsat_stop(&e.budget, status, *out ? nsat_stats(*out).stop : NSAT_STOP_ALLOCATION);
        goto done;
    }
    e.solver = *out;
    width = values > q->target ? values : q->target;
    scratch = nsat_alloc(&e.budget, width, sizeof(*scratch));
    if (scratch && write_multivalue(&e, q, (int)values, scratch) &&
        e.emitted != nsat_verify_original_count(checker))
        nsat_stop(&e.budget, NSAT_ERROR, NSAT_STOP_PROOF);
    goto done;
invalid:
    nsat_stop(&e.budget, NSAT_ERROR, NSAT_STOP_INPUT);
done:
    nsat_free(&e.budget, scratch, width * sizeof(*scratch));
    *stop = e.budget.stats.stop;
    return e.budget.status;
}

NsatStatus nsat_multivalue_solve(const NsatMultivalueInput *q, NsatLimits limits,
                                NsatSolver **solver, NsatVerifier **checker,
                                NsatStopReason *stop)
{
    if (!solver || !checker || !stop) return NSAT_ERROR;
    *solver = NULL;
    *checker = NULL;
    *stop = NSAT_STOP_INPUT;
    if (!q) return NSAT_ERROR;
    NsatStatus status = nsat_verify_init_multivalue(checker, q->on, q->off, q->levels,
        q->inputs, q->positives, q->negatives, q->target, limits);
    if (status != NSAT_OK) {
        if (*checker) *stop = nsat_verify_stats(*checker).stop;
        else if (limits.shared && limits.shared->stop) *stop = limits.shared->stop;
        return status;
    }
    status = nsat_multivalue_load_solver(solver, q, *checker, limits, stop);
    if (status == NSAT_OK) {
        /* The final 2*p*k variables are assignment and prefix blocks. */
        int assignments = (int)(q->positives * q->target);
        NsatSearchOptions options = {
            .priority_first = nsat_verify_variables(*checker) - 2 * assignments + 1,
            .priority_count = assignments, .row_width = (int)q->target,
            .sticky_priority = true, .cache_root_units = true
        };
        status = nsat_solve(*solver, *checker, options);
        *stop = nsat_stats(*solver).stop;
        if (!*stop) *stop = nsat_verify_stats(*checker).stop;
    }
    return status;
}
