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
    if (!nsat_mul(width, (size_t)table->positives + (size_t)table->negatives + 1, &cells))
        return NSAT_UNKNOWN;
    int *points = nsat_alloc(&budget, cells, sizeof(*points));
    if (!points) return budget.status;
    NsatStatus status = NSAT_UNKNOWN;
    int *levels = points + cells - width;
    bool multivalue = false;
    if (!nsat_work(&budget, width)) goto done;
    for (size_t c = 0; c < width; ++c) levels[c] = 1;
    for (int side = 0; side < 2; ++side) {
        const int *source = side ? table->off : table->on;
        int rows = side ? table->negatives : table->positives;
        int *destination = points + (side ? on_cells : 0);
        for (int r = 0; r < rows; ++r) {
            for (size_t c = 0; c < width; ++c) {
                if (!nsat_work(&budget, 1)) goto done;
                int value = source[c * rows + r];
                if (value < 0 || value == INT_MAX) goto done;
                if (value > 1) multivalue = true;
                if (levels[c] < value + 1) levels[c] = value + 1;
                destination[(size_t)r * width + c] = value + 1;
            }
        }
    }
    if (multivalue) {
        /* Categories absent from the care rows need not be represented:
         * fixing such a value cannot cover an ON row. Values are equality
         * literals, never unions of categories within a condition. */
        NsatMultivalueInput query = {.on = points, .off = points + on_cells,
            .levels = levels, .inputs = width, .positives = (size_t)table->positives,
            .negatives = (size_t)table->negatives, .target = (size_t)(incumbent - 1)};
        NsatStopReason stop;
        status = nsat_multivalue_certify(&query, limits, &table->pair_work, &stop);
        goto done;
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
    status = nsat_cover_solve_default(&query, limits, &solver, &checker, &stop, &stats);
    nsat_destroy(solver);
    nsat_verify_destroy(checker);
done:
    nsat_free(&budget, points, cells * sizeof(*points));
    return status;
}

/* These encoders already supply the same originals to both consumers. Seal
 * the index first so its growth never overlaps preallocated solver arrays. */
static bool load_originals(NsatVerifier *checker, NsatSolver *solver) {
    return nsat_verify_seal_originals(checker) == NSAT_OK &&
        nsat_load_originals(solver, checker) == NSAT_OK;
}

/*
 * One variable per condition records whether the new term fixes that
 * condition to the anchor row's value. Each OFF row supplies a clause
 * requiring a differing fixed condition. A known PI supplies a clause
 * requiring at least one of its fixed conditions to be left free.
 */
static bool anchor_complete(const QcaCertification *table, int row, const int *chart,
                            const int *implicants, int columns, size_t originals,
                            int *anchor, int *lits, NsatLimits limits) {
    NsatVerifier *checker = NULL;
    NsatSolver *solver = NULL;
    bool complete = false;
    int width = table->conditions;
    NsatBudget budget = {.limits = limits};
    if (nsat_verify_init_cnf(&checker, width, limits) != NSAT_OK ||
        nsat_init(&solver, width, limits) != NSAT_OK ||
        nsat_verify_reserve_originals(checker, originals) != NSAT_OK ||
        !nsat_work(&budget, (uint64_t)width)) goto done;
    for (int c = 0; c < width; ++c)
        anchor[c] = table->on[(size_t)c * table->positives + row];
    for (int r = 0; r < table->negatives; ++r) {
        if (!nsat_work(&budget, (uint64_t)width)) goto done;
        size_t count = 0;
        for (int c = 0; c < width; ++c)
            if (table->off[(size_t)c * table->negatives + r] != anchor[c]) lits[count++] = c + 1;
        if (nsat_verify_add_original(checker, lits, count) != NSAT_OK) goto done;
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
        if (c == width && nsat_verify_add_original(checker, lits, count) != NSAT_OK) goto done;
    }
    complete = load_originals(checker, solver) &&
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
    /* Every OFF row emits one clause. Reserve the capacity ordinary growth
     * reaches for those clauses, so later PI clauses grow at the same points.
     * If even all columns fit, the upper bound needs no further growth. */
    size_t originals = 8;
    while (originals < (size_t)table->negatives) originals *= 2;
    size_t maximum = (size_t)table->negatives + (size_t)columns;
    if (maximum < originals) originals = maximum;
    int r = 0;
    for (; r < table->positives; ++r)
        if (!anchor_complete(table, r, chart, implicants, columns, originals, scratch,
                             scratch + width, limits)) break;
    nsat_free(&budget, scratch, 2 * width * sizeof(*scratch));
    return r == table->positives;
}

/* Slot zero is distinguished. The assignment variables ensure each ON row
 * is covered and slot zero covers an ON row. The remaining slots are free
 * to use any legal terms, including terms already present in the chart. */
static bool add_minimum_cover(const QcaCertification *table, int minimum, int literal_vars,
                              NsatVerifier *checker, int *lits, NsatBudget *budget) {
    int width = table->conditions, positives = table->positives;
    for (int t = 0; t < minimum; ++t) {
        for (int c = 0; c < width; ++c) {
            int zero = nsat_cover_selector(width, t, c, 1);
            if (nsat_verify_add_original(checker, (int[]){-zero, -(zero + 1)}, 2) != NSAT_OK)
                return false;
        }
        for (int r = 0; r < table->negatives; ++r) {
            if (!nsat_work(budget, (uint64_t)width)) return false;
            for (int c = 0; c < width; ++c) {
                int value = table->off[(size_t)c * table->negatives + r];
                if (value != 0 && value != 1) return false;
                lits[c] = nsat_cover_selector(width, t, c, 2 - value);
            }
            if (nsat_verify_add_original(checker, lits, width) != NSAT_OK) return false;
        }
    }
    for (int r = 0; r < positives; ++r) {
        int first = literal_vars + r * minimum + 1;
        for (int t = 0; t < minimum; ++t) {
            for (int c = 0; c < width; ++c) {
                int value = table->on[(size_t)c * positives + r];
                if (value != 0 && value != 1) return false;
                int opposite = nsat_cover_selector(width, t, c, 2 - value);
                if (nsat_verify_add_original(checker, (int[]){-(first + t), -opposite}, 2) != NSAT_OK)
                    return false;
            }
            lits[t] = first + t;
        }
        if (nsat_verify_add_original(checker, lits, minimum) != NSAT_OK) return false;
    }
    if (!nsat_work(budget, (uint64_t)positives)) return false;
    for (int r = 0; r < positives; ++r) lits[r] = literal_vars + r * minimum + 1;
    return nsat_verify_add_original(checker, lits, positives) == NSAT_OK;
}

/* Order first assigned ON rows among slots 1..minimum-1. Slot zero keeps
 * its separate missing-PI role and need not precede any other slot.
 * A cover can assign each row to one covering term, retain an assignment
 * to slot zero, and permute the remaining used slots by their first row.
 * Unused slots follow them. Thus these original encoding constraints
 * preserve satisfiability; they are not RUP consequences of the old CNF. */
static bool order_minimum_slots(const QcaCertification *table, int minimum,
                                int literal_vars, NsatVerifier *checker) {
    if (minimum < 3) return true;
    const int prefixes = minimum - 2;
    const int first_prefix = literal_vars + table->positives * minimum + 1;
    for (int r = 0; r < table->positives; ++r) {
        for (int t = 1; t < minimum; ++t) {
            int assigned = literal_vars + r * minimum + t + 1;
            if (t < minimum - 1) {
                int used = first_prefix + r * prefixes + t - 1;
                if (nsat_verify_add_original(checker, (int[]){-assigned, used}, 2) != NSAT_OK)
                    return false;
                if (!r) {
                    if (nsat_verify_add_original(checker, (int[]){-used, assigned}, 2) != NSAT_OK)
                        return false;
                } else {
                    int previous = used - prefixes;
                    if (nsat_verify_add_original(checker, (int[]){-previous, used}, 2) != NSAT_OK ||
                        nsat_verify_add_original(checker, (int[]){-used, previous, assigned}, 3) != NSAT_OK)
                        return false;
                }
            }
            if (t > 1) {
                if (!r) {
                    if (nsat_verify_add_original(checker, (int[]){-assigned}, 1) != NSAT_OK) return false;
                } else {
                    int earlier = first_prefix + (r - 1) * prefixes + t - 2;
                    if (nsat_verify_add_original(checker, (int[]){-assigned, earlier}, 2) != NSAT_OK)
                        return false;
                }
            }
        }
    }
    return true;
}

static NsatStatus minimum_distance_one(const QcaCertification *table, int minimum,
                                       NsatSolver *solver, NsatVerifier *checker,
                                       NsatLimits limits)
{
    size_t cells, width = (size_t)table->conditions;
    if (!nsat_mul(width, (size_t)table->positives + (size_t)table->negatives, &cells))
        return NSAT_ERROR;
    NsatBudget budget = {.limits = limits};
    int *points = nsat_alloc(&budget, cells, sizeof(*points));
    if (!points) return budget.status;
    size_t on_cells = width * (size_t)table->positives;
    NsatStatus status = NSAT_OK;
    for (int side = 0; side < 2; ++side) {
        const int *source = side ? table->off : table->on;
        int rows = side ? table->negatives : table->positives;
        int *destination = points + (side ? on_cells : 0);
        for (int r = 0; r < rows; ++r) {
            for (size_t c = 0; c < width; ++c) {
                if (!nsat_work(&budget, 1)) {
                    status = budget.status;
                    goto done;
                }
                destination[(size_t)r * width + c] = source[c * rows + r] + 1;
            }
        }
    }
    NsatCoverInput query = {.on = points, .off = points + on_cells, .inputs = width,
        .positives = (size_t)table->positives, .negatives = (size_t)table->negatives,
        .target = (size_t)minimum};
    NsatCoverDistanceStats stats;
    status = nsat_cover_distance_one(&query, solver, checker, limits, &stats);
done:
    nsat_free(&budget, points, cells * sizeof(*points));
    return status;
}

static const uint64_t slot_zero_work_limit = UINT64_C(4000000);

static bool minimum_slot_zero_complete(QcaCertification *table, const int *chart,
                                       const int *implicants, int columns)
{
    if (table->slot_zero_work >= slot_zero_work_limit) return false;

    /* Share the memory ceiling, but keep optional work and stops out of
     * the main certification allowance. All probe allocations are freed
     * before the fallback can run. */
    QcaCertification probe = *table;
    probe.resources.work = table->slot_zero_work;
    probe.resources.work_limit = slot_zero_work_limit;
    probe.resources.stop = NSAT_STOP_NONE;
    bool complete = qca_check_complete_chart(&probe, chart, implicants, columns);
    table->slot_zero_work = probe.resources.work;
    table->resources.peak_memory = probe.resources.peak_memory;
    return complete;
}

#define QCA_MINIMUM_PAIR_WORK_LIMIT UINT64_C(4000000)
#define QCA_MINIMUM_ROW_RETRY_WORK_LIMIT UINT64_C(16000000)
#define QCA_MINIMUM_ROW_RETRY_MIN_REMAINING UINT64_C(64000000)
enum { MAX_MINIMUM_PAIR_ON = 96, MAX_MINIMUM_PAIR_OFF = 128,
       MAX_MINIMUM_PAIR_SCANS = 100000, MAX_MINIMUM_PAIR_CLAUSES = 4096 };

static bool pack_minimum_pair_rows(const int *points, int width, int rows,
                                   uint64_t *packed, NsatBudget *budget) {
    for (int row = 0; row < rows; ++row) {
        if (!nsat_work(budget, (uint64_t)width)) return false;
        uint64_t bits = 0;
        for (int c = 0; c < width; ++c) {
            bits |= (uint64_t)points[(size_t)c * (size_t)rows + (size_t)row] << c;
        }
        packed[row] = bits;
    }
    return true;
}

static NsatStatus minimum_pair_clauses(QcaCertification *table, int minimum,
                                       NsatSolver *solver, NsatVerifier *checker) {
    int width = table->conditions, positives = table->positives, negatives = table->negatives;
    size_t added = 0;

    NsatBudget budget = {.limits = {.shared = &table->resources}};
    uint64_t on[MAX_MINIMUM_PAIR_ON], off[MAX_MINIMUM_PAIR_OFF];
    if (!pack_minimum_pair_rows(table->on, width, positives, on, &budget) ||
        !pack_minimum_pair_rows(table->off, width, negatives, off, &budget)) return budget.status;

    /* An OFF row in the smallest cube containing both ON rows rules out
     * assigning both to one slot. On each condition, an ON implication
     * rules out the selector that would exclude the witnessed OFF row. */
    int literal_vars = 2 * width * minimum;
    NsatClauseId slot_clauses = (NsatClauseId)width + (NsatClauseId)negatives;
    NsatClauseId first_on = 1 + (NsatClauseId)minimum * slot_clauses;
    NsatClauseId row_clauses = 1 + (NsatClauseId)minimum * (NsatClauseId)width;
    size_t scanned = 0;
    for (int first = 0; first < positives - 1; ++first) {
        for (int second = first + 1; second < positives; ++second) {
            int witness = 0;
            for (; witness < negatives; ++witness) {
                if (scanned == MAX_MINIMUM_PAIR_SCANS) return NSAT_OK;
                ++scanned;
                if (!nsat_work(&budget, 1)) return budget.status;
                if (((on[first] ^ off[witness]) & (on[second] ^ off[witness])) == 0) break;
            }
            if (witness == negatives) continue;

            /* Ordinary slot t needs t-1 earlier first uses. Slot zero is
             * distinguished and can have its first use at any ON row. */
            int slots = first + 2;
            if (slots > minimum) slots = minimum;
            uint64_t differs = on[first] ^ off[witness];
            for (int slot = 0; slot < slots; ++slot) {
                if (!nsat_work(&budget, (uint64_t)width)) return budget.status;
                NsatClauseId hints[65];
                for (int c = 0; c < width; ++c) {
                    int row = differs & (UINT64_C(1) << c) ? second : first;
                    hints[c] = first_on + (NsatClauseId)row * row_clauses +
                        (NsatClauseId)slot * (NsatClauseId)width + (NsatClauseId)c;
                }
                hints[width] = 1 + (NsatClauseId)slot * slot_clauses +
                    (NsatClauseId)width + (NsatClauseId)witness;
                int clause[] = {
                    -(literal_vars + first * minimum + slot + 1),
                    -(literal_vars + second * minimum + slot + 1)
                };
                NsatStatus status = nsat_add_root_rup(solver, checker, clause, 2,
                    hints, (size_t)width + 1, NULL);
                if (status != NSAT_OK) return status;
                if (++added == MAX_MINIMUM_PAIR_CLAUSES) return NSAT_OK;
            }
        }
    }
    return NSAT_OK;
}

static NsatStatus search_minimum_models(QcaCertification *table, const int *implicants,
                                        int columns, int minimum, int variables,
                                        bool ordering, bool distance, bool pairs, bool row_guided) {
    int width = table->conditions, positives = table->positives;
    const int literal_vars = 2 * width * minimum;
    const size_t capacity = width > positives ? (size_t)width : (size_t)positives;
    NsatLimits limits = {.shared = &table->resources};
    NsatBudget budget = {.limits = limits};
    int *lits = nsat_alloc(&budget, capacity, sizeof(*lits));
    NsatVerifier *checker = NULL;
    NsatSolver *solver = NULL;
    NsatStatus status = NSAT_UNKNOWN;
    if (!lits || nsat_verify_init_cnf(&checker, variables, limits) != NSAT_OK ||
        nsat_init(&solver, variables, limits) != NSAT_OK)
        goto done;
    /* Slot legality, row coverage, a nonempty slot zero, and the known PIs. */
    uint64_t slots = (uint64_t)minimum, rows = (uint64_t)positives;
    uint64_t originals = slots * ((uint64_t)width + table->negatives) +
        rows * (slots * width + 1) + 1 + columns;
    if (ordering && minimum > 2) originals += (4 * rows - 1) * (slots - 2);
    if (originals <= SIZE_MAX / sizeof(void *) &&
        nsat_verify_reserve_originals(checker, (size_t)originals) != NSAT_OK)
        goto done;
    if (!add_minimum_cover(table, minimum, literal_vars, checker, lits, &budget))
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
        if (nsat_verify_add_original(checker, lits, count) != NSAT_OK) goto done;
    }
    nsat_free(&budget, lits, capacity * sizeof(*lits));
    lits = NULL;
    if (ordering && !order_minimum_slots(table, minimum, literal_vars, checker))
        goto done;
    if (!load_originals(checker, solver) ||
        (distance && minimum_distance_one(table, minimum, solver, checker, limits) != NSAT_OK))
        goto done;
    if (pairs && minimum_pair_clauses(table, minimum, solver, checker) != NSAT_OK)
        goto done;
    NsatSearchOptions options = {.cache_root_units = true};
    if (row_guided) {
        options.priority_first = literal_vars + 1;
        options.priority_count = positives * minimum;
        options.row_width = minimum;
    }
    status = nsat_solve(solver, checker, options);
done:
    nsat_destroy(solver);
    nsat_verify_destroy(checker);
    nsat_free(&budget, lits, capacity * sizeof(*lits));
    return status;
}

/* Each probe owns its work counter. The caller keeps its live memory and
 * peak accounting, and the ordinary search keeps all its remaining work. */
static NsatStatus probe_minimum_pairs(QcaCertification *table, const int *implicants,
                                    int columns, int minimum, int variables,
                                    bool row_guided) {
    uint64_t *spent = row_guided ? &table->row_retry_work : &table->pair_work;
    uint64_t allowance = row_guided ? QCA_MINIMUM_ROW_RETRY_WORK_LIMIT : QCA_MINIMUM_PAIR_WORK_LIMIT;
    QcaCertification probe = *table;
    probe.resources.work = *spent;
    probe.resources.work_limit = allowance;
    probe.resources.stop = NSAT_STOP_NONE;
    NsatStatus status = search_minimum_models(&probe, implicants, columns, minimum,
        variables, true, true, true, row_guided);
    *spent = probe.resources.stop == NSAT_STOP_WORK ? allowance : probe.resources.work;
    table->resources.peak_memory = probe.resources.peak_memory;
    return status;
}

static bool check_minimum_models(QcaCertification *table, const int *implicants,
                                 int columns, int minimum, bool ordering, bool distance,
                                 const int *chart) {
    if (!available(table) || columns <= 0 || minimum < 1 || minimum > table->positives)
        return false;
    const int width = table->conditions, positives = table->positives;
    if ((size_t)columns > SIZE_MAX / sizeof(int) / (size_t)width ||
        width > (INT_MAX - positives) / 2 ||
        minimum > INT_MAX / (2 * width + positives))
        return false;
    const int literal_vars = 2 * width * minimum;
    int variables = literal_vars + positives * minimum;
    if (ordering && minimum > 2) {
        size_t extra = (size_t)positives * (size_t)(minimum - 2);
        if (extra > (size_t)(INT_MAX - variables)) return false;
        variables += (int)extra;
    }
    if (chart && minimum > 1 && minimum_slot_zero_complete(table, chart, implicants, columns))
        return true;

    /* Only a checked answer replaces the ordinary search. The trial shares
     * the memory ceiling and releases both solver and checker before retry. */
    if (ordering && distance && minimum >= 2 && width <= 64 &&
        positives <= MAX_MINIMUM_PAIR_ON && table->negatives <= MAX_MINIMUM_PAIR_OFF &&
        table->pair_work < QCA_MINIMUM_PAIR_WORK_LIMIT) {
        NsatStatus status = probe_minimum_pairs(table, implicants, columns, minimum, variables, false);
        if (status == NSAT_UNSAT || status == NSAT_SAT) return status == NSAT_UNSAT;
        if (status == NSAT_UNKNOWN && table->row_retry_work < QCA_MINIMUM_ROW_RETRY_WORK_LIMIT &&
            table->resources.work_limit - table->resources.work >= QCA_MINIMUM_ROW_RETRY_MIN_REMAINING) {
            status = probe_minimum_pairs(table, implicants, columns, minimum, variables, true);
            if (status == NSAT_UNSAT || status == NSAT_SAT) return status == NSAT_UNSAT;
        }
    }
    return search_minimum_models(table, implicants, columns, minimum, variables,
        ordering, distance, false, false) == NSAT_UNSAT;
}

bool qca_check_minimum_models(QcaCertification *table, const int *chart,
                              const int *implicants, int columns, int minimum)
{
    return check_minimum_models(table, implicants, columns, minimum, true, true, chart);
}
