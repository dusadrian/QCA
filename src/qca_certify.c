#include "qca_certify.h"
#include "native_sat/cover_encoding.h"
#include "native_sat/ipc_weighted_cover.h"
#include "native_sat/native_sat_support.h"

static bool available(const QcaCertification *table)
{
    return table->conditions > 0 && table->positives > 0 && table->negatives > 0 &&
        (size_t)table->conditions <= SIZE_MAX / sizeof(int) / (size_t)table->positives &&
        (size_t)table->conditions <= SIZE_MAX / sizeof(int) / (size_t)table->negatives &&
        table->resources.work < table->resources.work_limit;
}

/* Fix every condition on which the positive rows agree. If this term
 * contains a negative row, every one-term cover contains it as well. */
static NsatStatus one_term(const QcaCertification *table, NsatBudget *budget)
{
    const size_t width = (size_t)table->conditions;
    int *fixed = nsat_alloc(budget, width, sizeof(*fixed));
    if (!fixed) return budget->status;
    int count = 0;
    NsatStatus status = NSAT_UNKNOWN;
    if (!nsat_work(budget, width)) goto done;
    for (int c = 0; c < table->conditions; ++c) {
        const int *column = table->on + (size_t)c * table->positives;
        int r = 1;
        for (; r < table->positives; ++r) {
            if (!nsat_work(budget, 1)) goto done;
            if (column[r] != column[0]) break;
        }
        if (r == table->positives) fixed[count++] = c;
    }
    for (int r = 0; r < table->negatives; ++r) {
        int i = 0;
        for (; i < count; ++i) {
            if (!nsat_work(budget, 1)) goto done;
            int c = fixed[i];
            if (table->off[(size_t)c * table->negatives + r] !=
                table->on[(size_t)c * table->positives]) break;
        }
        if (i == count) {
            status = NSAT_UNSAT;
            goto done;
        }
    }
    status = NSAT_SAT;
done:
    nsat_free(budget, fixed, width * sizeof(*fixed));
    return status;
}

NsatStatus qca_check_smaller_model(QcaCertification *table, int incumbent)
{
    if (incumbent >= 2 && table->verified_minimum == incumbent) return NSAT_UNSAT;
    if (!available(table) || incumbent < 2) return NSAT_UNKNOWN;
    NsatLimits limits = {.shared = &table->resources};
    NsatBudget budget = {.limits = limits};
    NsatStatus single = one_term(table, &budget);
    if (incumbent == 2 || single != NSAT_UNSAT) return single;

    size_t width = (size_t)table->conditions, on_cells = (size_t)table->positives * width;
    size_t cells;
    if (!nsat_mul(width, (size_t)table->positives + (size_t)table->negatives, &cells))
        return NSAT_UNKNOWN;
    int *points = nsat_alloc(&budget, cells, sizeof(*points));
    if (!points) return budget.status;
    NsatStatus status = NSAT_UNKNOWN;
    for (int side = 0; side < 2; ++side) {
        const int *source = side ? table->off : table->on;
        int rows = side ? table->negatives : table->positives;
        int *destination = points + (side ? on_cells : 0);
        for (int r = 0; r < rows; ++r) {
            for (size_t c = 0; c < width; ++c) {
                if (!nsat_work(&budget, 1)) goto done;
                int value = source[c * rows + r];
                if (value != 0 && value != 1) goto done;
                destination[(size_t)r * width + c] = value + 1;
            }
        }
    }
    NsatCoverInput query = {.on = points, .off = points + on_cells,
        .inputs = width, .positives = (size_t)table->positives,
        .negatives = (size_t)table->negatives, .target = (size_t)incumbent};
    const uint64_t left = table->resources.work_limit - table->resources.work;
    if (left >= 2) {
        /* Let weighted verification use up to half the remaining work.
         * Reserve the rest for SAT and chart checks in this minimization. */
        const uint64_t weighted_cap = UINT64_C(8000000);
        NsatLimits pricing = limits;
        pricing.work = left / 2 < weighted_cap ? left / 2 : weighted_cap;
        IpcWeightedCoverResult weighted = ipc_weighted_cover_seed_then_reprice(
            &query, 1, pricing, 4, 256);
        if (weighted.status != NSAT_ERROR) {
            const bool checked = weighted.check_status == NSAT_UNSAT;
            if (checked && weighted.lower == incumbent) {
                status = NSAT_UNSAT;
            } else if (weighted.cover && weighted.terms < incumbent) {
                /* The chart still has to find its own model. Reuse this proof
                 * if a later plateau reaches the independently checked bound. */
                if (checked && weighted.terms == weighted.lower)
                    table->verified_minimum = weighted.lower;
                status = NSAT_SAT;
            }
        }
        free(weighted.cover);
        if (status != NSAT_UNKNOWN) goto done;
    }
    query.target = (size_t)(incumbent - 1);
    NsatSolver *solver = NULL;
    NsatVerifier *checker = NULL;
    NsatStopReason stop;
    NsatCoverDomainStats stats;
    status = nsat_cover_solve(&query, NSAT_COVER_ORDERING, limits, &solver, &checker, &stop, &stats);
    nsat_destroy(solver);
    nsat_verify_destroy(checker);
done:
    nsat_free(&budget, points, cells * sizeof(*points));
    return status;
}

static bool add_clause(NsatVerifier *checker, NsatSolver *solver, const int *lits, size_t count)
{
    return nsat_verify_add_original(checker, lits, count) == NSAT_OK &&
        nsat_add_original(solver, lits, count) == NSAT_OK;
}

/*
 * One variable per condition records whether the new term fixes that
 * condition to the anchor row's value. Each OFF row supplies a clause
 * requiring a differing fixed condition. A known PI supplies a clause
 * requiring at least one of its fixed conditions to be left free.
 */
static bool anchor_complete(const QcaCertification *table, int row, const int *chart,
                            const int *implicants, int columns, int *anchor, int *lits,
                            NsatLimits limits)
{
    NsatVerifier *checker = NULL;
    NsatSolver *solver = NULL;
    bool complete = false;
    int width = table->conditions;
    NsatBudget budget = {.limits = limits};
    if (nsat_verify_init_cnf(&checker, width, limits) != NSAT_OK ||
        nsat_init(&solver, width, limits) != NSAT_OK ||
        !nsat_work(&budget, (uint64_t)width)) goto done;
    for (int c = 0; c < width; ++c)
        anchor[c] = table->on[(size_t)c * table->positives + row];
    for (int r = 0; r < table->negatives; ++r) {
        if (!nsat_work(&budget, (uint64_t)width)) goto done;
        size_t count = 0;
        for (int c = 0; c < width; ++c)
            if (table->off[(size_t)c * table->negatives + r] != anchor[c]) lits[count++] = c + 1;
        if (!add_clause(checker, solver, lits, count)) goto done;
    }
    for (int pi = 0; pi < columns; ++pi) {
        if (!nsat_work(&budget, 1)) goto done;
        if (!chart[(size_t)pi * table->positives + row]) continue;
        size_t count = 0;
        int c = 0;
        for (; c < width; ++c) {
            int literal = implicants[(size_t)pi * width + c];
            if (!nsat_work(&budget, 1)) goto done;
            if (!literal) continue;
            if (literal - 1 != anchor[c]) break;
            lits[count++] = -(c + 1);
        }
        if (c == width && !add_clause(checker, solver, lits, count)) goto done;
    }
    complete = nsat_verify_seal_originals(checker) == NSAT_OK &&
        nsat_solve(solver, checker, (NsatSearchOptions){0}) == NSAT_UNSAT;
done:
    nsat_destroy(solver);
    nsat_verify_destroy(checker);
    return complete;
}

bool qca_check_complete_chart(QcaCertification *table, const int *chart,
                              const int *implicants, int columns)
{
    if (!available(table) || columns <= 0) return false;
    const size_t width = (size_t)table->conditions;
    if ((size_t)columns > SIZE_MAX / sizeof(int) / width ||
        (size_t)columns > SIZE_MAX / sizeof(int) / (size_t)table->positives)
        return false;
    NsatLimits limits = {.shared = &table->resources};
    NsatBudget budget = {.limits = limits};
    int *scratch = nsat_alloc(&budget, width, 2 * sizeof(*scratch));
    if (!scratch) return false;
    int r = 0;
    for (; r < table->positives; ++r)
        if (!anchor_complete(table, r, chart, implicants, columns, scratch,
                             scratch + width, limits)) break;
    nsat_free(&budget, scratch, 2 * width * sizeof(*scratch));
    return r == table->positives;
}

/* Slot zero is distinguished. The assignment variables ensure each ON row
 * is covered and slot zero covers an ON row. The remaining slots are free
 * to use any legal terms, including terms already present in the chart. */
static bool add_minimum_cover(const QcaCertification *table, int minimum, int literal_vars,
                              NsatVerifier *checker, NsatSolver *solver, int *lits,
                              NsatBudget *budget)
{
    int width = table->conditions, positives = table->positives;
    for (int t = 0; t < minimum; ++t) {
        for (int c = 0; c < width; ++c) {
            int zero = nsat_cover_selector(width, t, c, 1);
            if (!add_clause(checker, solver, (int[]){-zero, -(zero + 1)}, 2))
                return false;
        }
        for (int r = 0; r < table->negatives; ++r) {
            if (!nsat_work(budget, (uint64_t)width)) return false;
            for (int c = 0; c < width; ++c) {
                int value = table->off[(size_t)c * table->negatives + r];
                if (value != 0 && value != 1) return false;
                lits[c] = nsat_cover_selector(width, t, c, 2 - value);
            }
            if (!add_clause(checker, solver, lits, width)) return false;
        }
    }
    for (int r = 0; r < positives; ++r) {
        int first = literal_vars + r * minimum + 1;
        for (int t = 0; t < minimum; ++t) {
            for (int c = 0; c < width; ++c) {
                int value = table->on[(size_t)c * positives + r];
                if (value != 0 && value != 1) return false;
                int opposite = nsat_cover_selector(width, t, c, 2 - value);
                if (!add_clause(checker, solver, (int[]){-(first + t), -opposite}, 2))
                    return false;
            }
            lits[t] = first + t;
        }
        if (!add_clause(checker, solver, lits, minimum)) return false;
    }
    if (!nsat_work(budget, (uint64_t)positives)) return false;
    for (int r = 0; r < positives; ++r) lits[r] = literal_vars + r * minimum + 1;
    return add_clause(checker, solver, lits, positives);
}

bool qca_check_minimum_models(QcaCertification *table, const int *implicants,
                              int columns, int minimum)
{
    if (!available(table) || columns <= 0 || minimum < 1 || minimum > table->positives)
        return false;
    const int width = table->conditions, positives = table->positives;
    if ((size_t)columns > SIZE_MAX / sizeof(int) / (size_t)width ||
        width > (INT_MAX - positives) / 2 ||
        minimum > INT_MAX / (2 * width + positives))
        return false;
    const int literal_vars = 2 * width * minimum;
    const int variables = literal_vars + positives * minimum;
    const size_t capacity = width > positives ? (size_t)width : (size_t)positives;
    NsatLimits limits = {.shared = &table->resources};
    NsatBudget budget = {.limits = limits};
    int *lits = nsat_alloc(&budget, capacity, sizeof(*lits));
    NsatVerifier *checker = NULL;
    NsatSolver *solver = NULL;
    bool complete = false;
    if (!lits || nsat_verify_init_cnf(&checker, variables, limits) != NSAT_OK ||
        nsat_init(&solver, variables, limits) != NSAT_OK ||
        !add_minimum_cover(table, minimum, literal_vars, checker, solver, lits, &budget))
        goto done;

    /* A missing PI contains no known PI as a subterm. Conversely, a legal
     * slot-zero term satisfying these clauses can be expanded to a missing
     * PI by deleting literals. This preserves the cover at the same size. */
    for (int pi = 0; pi < columns; ++pi) {
        if (!nsat_work(&budget, (uint64_t)width)) goto done;
        size_t count = 0;
        for (int c = 0; c < width; ++c) {
            int value = implicants[(size_t)pi * width + c];
            if (!value) continue;
            if (value != 1 && value != 2) goto done;
            lits[count++] = -nsat_cover_selector(width, 0, c, value);
        }
        if (!add_clause(checker, solver, lits, count)) goto done;
    }
    complete = nsat_verify_seal_originals(checker) == NSAT_OK &&
        nsat_solve(solver, checker, (NsatSearchOptions){0}) == NSAT_UNSAT;
done:
    nsat_destroy(solver);
    nsat_verify_destroy(checker);
    nsat_free(&budget, lits, capacity * sizeof(*lits));
    return complete;
}
