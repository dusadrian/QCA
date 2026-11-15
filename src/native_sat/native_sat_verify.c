#include "native_sat_support.h"

typedef struct {
    size_t count;
    bool target_bound;
    int lits[];
} Clause;

typedef struct {
    NsatClauseId id;
    Clause *clause;
} ClauseEntry;

/* id=0 is unused. A nonzero ID with a NULL clause is a tombstone. */
struct NsatVerifier {
    NsatBudget budget;
    int variables;
    NsatClauseId last_id, originals;
    enum { LOADING, SEALED, SAT_CONCLUDED, CONCLUDED } phase;
    ClauseEntry *table;
    size_t capacity, occupied, live;
    signed char *values;
    int *assigned;
    size_t assigned_count;
    bool started;
    size_t cover_slots, cover_target;
};

static size_t hash(NsatClauseId id)
{
    id ^= id >> 33;
    id *= UINT64_C(0xff51afd7ed558ccd);
    return (size_t)(id ^ (id >> 33));
}

static ClauseEntry *lookup(const NsatVerifier *v, NsatClauseId id, NsatBudget *budget)
{
    if (!id || !v->capacity) return NULL;
    size_t i = hash(id) & (v->capacity - 1);
    while (!budget || nsat_work(budget, 1)) {
        ClauseEntry *e = &v->table[i];
        if (!e->id) return NULL;
        if (e->id == id) return e->clause ? e : NULL;
        i = (i + 1) & (v->capacity - 1);
    }
    return NULL;
}

static bool insert(NsatBudget *budget, ClauseEntry *table, size_t capacity,
                   ClauseEntry entry, size_t *occupied)
{
    size_t i = hash(entry.id) & (capacity - 1);
    while (nsat_work(budget, 1)) {
        if (!table[i].clause) {
            if (!table[i].id) ++*occupied;
            table[i] = entry;
            return true;
        }
        i = (i + 1) & (capacity - 1);
    }
    return false;
}

static bool rebuild(NsatVerifier *v, size_t capacity)
{
    NsatBudget *b = &v->budget;
    ClauseEntry *table = nsat_alloc(b, capacity, sizeof(*table));
    if (!table) return false;
    size_t occupied = 0;
    if (!nsat_copy(b, table, NULL, capacity * sizeof(*table))) goto fail;
    for (size_t i = 0; i < v->capacity; ++i) {
        if (!nsat_work(b, 1)) goto fail;
        if (v->table[i].clause && !insert(b, table, capacity, v->table[i], &occupied))
            goto fail;
    }
    nsat_free(b, v->table, v->capacity * sizeof(*v->table));
    v->table = table;
    v->capacity = capacity;
    v->occupied = occupied;
    ++b->stats.index_rebuilds;
    return true;
fail:
    nsat_free(b, table, capacity * sizeof(*table));
    return false;
}

static bool save(NsatVerifier *v, NsatClauseId id, const int *lits, size_t count)
{
    NsatBudget *b = &v->budget;
    if (!v->capacity || v->occupied >= v->capacity / 2) {
        size_t capacity = v->capacity ? v->capacity : 16;
        if (v->live >= capacity / 2) {
            if (!nsat_mul(capacity, 2, &capacity))
                return nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
        }
        if (!rebuild(v, capacity)) return false;
    }
    if (count > (SIZE_MAX - sizeof(Clause)) / sizeof(int))
        return nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
    size_t bytes = sizeof(Clause) + count * sizeof(int);
    Clause *clause = nsat_alloc(b, bytes, 1);
    if (!clause) return false;
    clause->count = count;
    clause->target_bound = false;
    if (!nsat_copy(b, clause->lits, lits, count * sizeof(int)) ||
        !nsat_poll(b) ||
        !insert(b, v->table, v->capacity, (ClauseEntry){id, clause}, &v->occupied)) {
        nsat_free(b, clause, bytes);
        return false;
    }
    ++v->live;
    v->last_id = id;
    return true;
}

NsatStatus nsat_verify_init_cnf(NsatVerifier **out, int variables, NsatLimits limits)
{
    if (!out) return NSAT_ERROR;
    *out = NULL;
    if (variables < 0) return NSAT_ERROR;
    NsatBudget budget = {.limits = limits};
    NsatVerifier *v = nsat_alloc(&budget, 1, sizeof(*v));
    if (!v) return budget.status;
    *v = (NsatVerifier){.budget = budget, .variables = variables};
    *out = v;
    size_t count = (size_t)variables + 1;
    v->values = nsat_alloc(&v->budget, count, sizeof(*v->values));
    if (!v->values || !nsat_copy(&v->budget, v->values, NULL, count))
        return v->budget.status;
    v->assigned = nsat_alloc(&v->budget, count, sizeof(*v->assigned));
    return v->budget.status;
}

static bool require_phase(NsatVerifier *v, int phase)
{
    if (v->budget.status != NSAT_OK) return false;
    if ((int)v->phase != phase)
        return nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_STATE);
    return nsat_poll(&v->budget);
}

static void assign(NsatVerifier *v, int var, int value)
{
    v->values[var] = (signed char)value;
    v->assigned[v->assigned_count++] = var;
}

static bool assume_false(NsatVerifier *v, const int *lits, size_t count)
{
    if (count > (size_t)v->variables || (count && !lits))
        return nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_INPUT);
    for (size_t i = 0; i < count; ++i) {
        if (!nsat_work(&v->budget, 1)) return false;
        int var = nsat_var(lits[i], v->variables);
        if (!var || v->values[var])
            return nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        assign(v, var, lits[i] < 0 ? 1 : -1);
    }
    return true;
}

static bool clear_values(NsatVerifier *v)
{
    while (v->assigned_count) {
        if (!nsat_work(&v->budget, 1)) return false;
        v->values[v->assigned[--v->assigned_count]] = 0;
    }
    return true;
}

NsatStatus nsat_verify_add_original(NsatVerifier *v, const int *lits, size_t count)
{
    if (!v) return NSAT_ERROR;
    if (!require_phase(v, LOADING)) return v->budget.status;
    if (v->last_id == UINT64_MAX)
        nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_INPUT);
    else if (assume_false(v, lits, count) && clear_values(v) &&
             save(v, v->last_id + 1, lits, count))
        v->originals = v->last_id;
    return v->budget.status;
}

NsatStatus nsat_verify_seal_originals(NsatVerifier *v)
{
    if (!v) return NSAT_ERROR;
    if (require_phase(v, LOADING)) v->phase = SEALED;
    return v->budget.status;
}

NsatStatus nsat_verify_begin(NsatVerifier *v)
{
    if (!v) return NSAT_ERROR;
    if (!require_phase(v, SEALED)) return v->budget.status;
    if (v->started || v->last_id != v->originals)
        nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_STATE);
    else
        v->started = true;
    return v->budget.status;
}

NsatStatus nsat_verify_derive(NsatVerifier *v, NsatClauseId id,
                              const int *lits, size_t count,
                              const NsatClauseId *hints, size_t hint_count)
{
    if (!v) return NSAT_ERROR;
    NsatBudget *b = &v->budget;
    if (!require_phase(v, SEALED)) return b->status;
    if (id <= v->last_id || !hints || !hint_count || hint_count > v->live) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_PROOF);
        return b->status;
    }
    if (!assume_false(v, lits, count)) return b->status;

    for (size_t h = 0; h < hint_count; ++h) {
        ClauseEntry *entry = lookup(v, hints[h], b);
        if (!entry) goto invalid;
        const Clause *clause = entry->clause;
        int unit = 0;
        for (size_t i = 0; i < clause->count; ++i) {
            if (!nsat_work(b, 1)) return b->status;
            int lit = clause->lits[i];
            int value = lit > 0 ? v->values[lit] : -v->values[-lit];
            if (value > 0 || (!value && unit)) goto invalid;
            if (!value) unit = lit;
        }
        if (h == hint_count - 1) {
            if (unit) goto invalid;
        } else {
            if (!unit) goto invalid;
            assign(v, unit > 0 ? unit : -unit, unit > 0 ? 1 : -1);
        }
    }
    if (clear_values(v) && save(v, id, lits, count)) ++b->stats.lemmas_checked;
    return b->status;
invalid:
    nsat_stop(b, NSAT_ERROR, NSAT_STOP_PROOF);
    return b->status;
}

NsatStatus nsat_verify_delete(NsatVerifier *v, NsatClauseId id)
{
    if (!v) return NSAT_ERROR;
    NsatBudget *b = &v->budget;
    if (!require_phase(v, SEALED)) return b->status;
    ClauseEntry *entry = id > v->originals ? lookup(v, id, b) : NULL;
    if (!entry || entry->clause->target_bound) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_PROOF);
        return b->status;
    }
    if (!nsat_poll(b)) return b->status;
    nsat_free(b, entry->clause, sizeof(Clause) + entry->clause->count * sizeof(int));
    entry->clause = NULL;
    --v->live;
    ++b->stats.clauses_deleted;
    if (v->capacity > 16 && v->live <= v->capacity / 8)
        rebuild(v, v->capacity / 2);
    return b->status;
}

NsatStatus nsat_verify_conclude_unsat(NsatVerifier *v, NsatClauseId id)
{
    if (!v) return NSAT_ERROR;
    if (!require_phase(v, SEALED)) return v->budget.status;
    ClauseEntry *entry = lookup(v, id, &v->budget);
    if (!entry || entry->clause->count)
        nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_PROOF);
    if (!nsat_poll(&v->budget)) return v->budget.status;
    v->phase = CONCLUDED;
    return NSAT_UNSAT;
}

NsatStatus nsat_verify_conclude_sat(NsatVerifier *v, const int *model, size_t count)
{
    if (!v) return NSAT_ERROR;
    NsatBudget *b = &v->budget;
    if (!require_phase(v, SEALED)) return b->status;
    if (count != (size_t)v->variables || (count && !model)) goto invalid;
    for (size_t i = 0; i < count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int lit = model[i];
        int var = nsat_var(lit, v->variables);
        if (!var || v->values[var]) goto invalid;
        assign(v, var, lit > 0 ? 1 : -1);
    }
    for (size_t i = 0; i < v->capacity; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        ClauseEntry entry = v->table[i];
        if (!entry.clause || (entry.id > v->originals && !entry.clause->target_bound))
            continue;
        bool satisfied = false;
        for (size_t j = 0; j < entry.clause->count; ++j) {
            if (!nsat_work(b, 1)) return b->status;
            int lit = entry.clause->lits[j];
            if ((lit > 0 ? v->values[lit] : -v->values[-lit]) > 0) {
                satisfied = true;
                break;
            }
        }
        if (!satisfied) goto invalid;
    }
    if (!clear_values(v) || !nsat_poll(b)) return b->status;
    v->phase = SAT_CONCLUDED;
    return NSAT_SAT;
invalid:
    nsat_stop(b, NSAT_ERROR, NSAT_STOP_PROOF);
    return b->status;
}

NsatStatus nsat_verify_tighten_cover(NsatVerifier *v, size_t target,
                                      NsatClauseId *id, int *literal)
{
    if (!v) return NSAT_ERROR;
    NsatBudget *b = &v->budget;
    if (!require_phase(v, SAT_CONCLUDED)) return b->status;
    if (!v->started || !v->cover_slots || !target || target >= v->cover_target ||
        v->last_id == UINT64_MAX || !id || !literal) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
        return b->status;
    }
    /* The last K variable IDs are the final row's prefix-use variables.
     * Forbid zero-based slot target. First use forces all later slots off. */
    int unit = -(v->variables - (int)v->cover_slots + 1 + (int)target);
    NsatClauseId next = v->last_id + 1;
    if (!save(v, next, &unit, 1)) return b->status;
    ClauseEntry *entry = lookup(v, next, b);
    if (!entry) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
        return b->status;
    }
    entry->clause->target_bound = true;
    v->cover_target = target;
    v->phase = SEALED;
    *id = next;
    *literal = unit;
    return NSAT_OK;
}

size_t nsat_verify_cover_target(const NsatVerifier *v) { return v->cover_target; }
size_t nsat_verify_cover_slots(const NsatVerifier *v) { return v->cover_slots; }

NsatStats nsat_verify_stats(const NsatVerifier *v) { return v->budget.stats; }
int nsat_verify_variables(const NsatVerifier *v) { return v->variables; }
NsatClauseId nsat_verify_original_count(const NsatVerifier *v) { return v->originals; }

bool nsat_verify_clause(const NsatVerifier *v, NsatClauseId id, const int **lits,
                        size_t *count)
{
    ClauseEntry *entry = v->budget.status == NSAT_OK ? lookup(v, id, NULL) : NULL;
    if (!entry || !lits || !count) return false;
    *lits = entry->clause->lits;
    *count = entry->clause->count;
    return true;
}

void nsat_verify_destroy(NsatVerifier *v)
{
    if (!v) return;
    for (size_t i = 0; i < v->capacity; ++i) free(v->table[i].clause);
    free(v->table);
    free(v->values);
    free(v->assigned);
    nsat_release_budget(&v->budget);
    free(v);
}

static bool valid_care(NsatVerifier *v, const int *on, const int *off,
                       size_t n, size_t on_count, size_t off_count)
{
    NsatBudget *b = &v->budget;
    size_t capacity = 16;
    while (off_count > capacity / 2) {
        if (!nsat_mul(capacity, 2, &capacity))
            return nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
    }
    const int **table = nsat_alloc(b, capacity, sizeof(*table));
    if (!table) return false;
    bool ok = nsat_copy(b, table, NULL, capacity * sizeof(*table));
    for (size_t side = 0; ok && side < 2; ++side) {
        const int *points = side ? on : off;
        size_t rows = side ? on_count : off_count;
        for (size_t row = 0; ok && row < rows; ++row) {
            const int *point = points + row * n;
            uint64_t key = UINT64_C(14695981039346656037);
            for (size_t i = 0; i < n; ++i) {
                if (!nsat_work(b, 1)) goto done;
                if (point[i] != 1 && point[i] != 2) goto invalid;
                key = (key ^ (unsigned)point[i]) * UINT64_C(1099511628211);
            }
            size_t pos = hash(key) & (capacity - 1);
            while (table[pos]) {
                if (!nsat_work(b, 1)) goto done;
                size_t i = 0;
                while (i < n && table[pos][i] == point[i]) {
                    if (!nsat_work(b, 1)) goto done;
                    ++i;
                }
                if (i == n) {
                    if (side) goto invalid;
                    break;
                }
                pos = (pos + 1) & (capacity - 1);
            }
            if (!side) table[pos] = point;
        }
    }
    goto done;
invalid:
    nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
done:
    nsat_free(b, table, capacity * sizeof(*table));
    return b->status == NSAT_OK;
}

/* Encoding v1. Number selector pairs by slot/condition, then assignment
 * and prefix variables by position/slot. Emit all ON blocks before OFF. */
static void cover_originals(NsatVerifier *v, const int *on, const int *off,
                            size_t n, size_t p, size_t q, size_t k,
                            const size_t *order, int *scratch)
{
    int selectors = (int)(2 * n * k);
    int prefix_offset = (int)(p * k);
    for (size_t r = 0; r < p; ++r) {
        const int *point = on + n * (order ? order[r] : r);
        int first = selectors + (int)(r * k) + 1;
        for (size_t t = 0; t < k; ++t) {
            if (!nsat_work(&v->budget, 1)) return;
            scratch[t] = first + (int)t;
        }
        if (nsat_verify_add_original(v, scratch, k) != NSAT_OK) return;
        for (size_t t = 0; t < k; ++t) {
            int a = first + (int)t;
            int u = a + prefix_offset;
            int base = (int)(2 * t * n);
            for (size_t i = 0; i < n; ++i) {
                int opposite = base + (int)(2 * i) + (3 - point[i]);
                if (nsat_verify_add_original(v, (int[]){-a, -opposite}, 2) != NSAT_OK)
                    return;
            }
            if (nsat_verify_add_original(v, (int[]){-a, u}, 2) != NSAT_OK)
                return;
            if (!r) {
                if (nsat_verify_add_original(v, (int[]){-u, a}, 2) != NSAT_OK)
                    return;
                if (t && nsat_verify_add_original(v, (int[]){-a}, 1) != NSAT_OK)
                    return;
            } else {
                int previous = u - (int)k;
                if (nsat_verify_add_original(v, (int[]){-previous, u}, 2) != NSAT_OK ||
                    nsat_verify_add_original(v, (int[]){-u, previous, a}, 3) != NSAT_OK)
                    return;
                if (t && nsat_verify_add_original(v, (int[]){-a, previous - 1}, 2) != NSAT_OK)
                    return;
            }
        }
    }
    for (size_t t = 0; t < k; ++t) {
        int base = (int)(2 * t * n);
        for (size_t i = 0; i < n; ++i) {
            int zero = base + (int)(2 * i) + 1;
            if (nsat_verify_add_original(v, (int[]){-zero, -(zero + 1)}, 2) != NSAT_OK)
                return;
        }
        for (size_t r = 0; r < q; ++r) {
            for (size_t i = 0; i < n; ++i) {
                if (!nsat_work(&v->budget, 1)) return;
                scratch[i] = base + (int)(2 * i) + (3 - off[r * n + i]);
            }
            if (nsat_verify_add_original(v, scratch, n) != NSAT_OK) return;
        }
    }
}

NsatStatus nsat_verify_init_cover(NsatVerifier **out, const int *on, const int *off,
                                  size_t n, size_t on_count, size_t off_count,
                                  size_t target, const size_t *row_order,
                                  size_t order_count, NsatLimits limits)
{
    if (!out) return NSAT_ERROR;
    *out = NULL;
    size_t cells;
    size_t half = (size_t)INT_MAX / 2;
    if (!on || !off || !n || !on_count || !off_count || !target ||
        (row_order ? order_count != on_count : order_count != 0) ||
        n > half || on_count > half - n || target > half / (n + on_count) ||
        !nsat_mul(on_count, n, &cells) || cells > SIZE_MAX / sizeof(int) ||
        !nsat_mul(off_count, n, &cells) || cells > SIZE_MAX / sizeof(int))
        return NSAT_ERROR;

    int variables = (int)(2 * (n + on_count) * target);
    NsatStatus status = nsat_verify_init_cnf(out, variables, limits);
    if (status != NSAT_OK) return status;
    NsatVerifier *v = *out;
    NsatBudget *b = &v->budget;
    if (row_order) {
        /* Use the otherwise empty scratch assignment array as a seen array. */
        for (size_t i = 0; i < on_count; ++i) {
            if (!nsat_work(b, 1)) return b->status;
            size_t row = row_order[i];
            if (row >= on_count || v->values[row]) {
                nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
                return b->status;
            }
            v->values[row] = 1;
        }
        if (!nsat_copy(b, v->values, NULL, on_count)) return b->status;
    }
    if (!valid_care(v, on, off, n, on_count, off_count)) return b->status;
    size_t width = n > target ? n : target;
    int *scratch = nsat_alloc(b, width, sizeof(*scratch));
    if (!scratch) return b->status;
    cover_originals(v, on, off, n, on_count, off_count, target, row_order, scratch);
    nsat_free(b, scratch, width * sizeof(*scratch));
    v->cover_slots = v->cover_target = target;
    return nsat_verify_seal_originals(v);
}
