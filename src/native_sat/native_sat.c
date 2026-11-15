#include "native_sat.h"
#include "native_sat_support.h"

typedef struct {
    size_t next;
    int mate; /* Present for binary watches. */
} Watch;

typedef struct {
    size_t offset, count;
    NsatClauseId id;
    Watch watch[2];
    unsigned glue;
    bool used;
} SolverClause;

enum {
    ANALYSED = 1, SOURCE = 2, REMOVABLE = 4, BLOCKED = 8, NEEDED = 16
};

typedef struct {
    int level, prev, next;
    size_t reason; /* One-based clause reference, zero at root and for decisions. */
    uint64_t stamp;
    signed char phase;
    unsigned char marks;
} Variable;

typedef struct {
    int variable;
    size_t next_literal;
} ReasonFrame;

struct NsatSolver {
    NsatBudget budget;
    int variables;
    bool started, model_valid;
    NsatVerifier *checker;
    NsatSearchOptions options;
    size_t pinned_count;
    SolverClause *clauses;
    size_t count, clause_capacity;
    int *words;
    size_t word_count, word_capacity, max_width;
    signed char *values; /* Adjacent entries for opposite literals. */
    NsatClauseId *units, *hints;
    NsatClauseId last_id;
    int *trail;
    size_t trail_size, propagated;
    size_t *heads;

    Variable *vars;
    size_t *level_start;
    NsatClauseId *grade;
    int level;
    int *learned, *touched;
    size_t learned_count, touched_count;
    ReasonFrame *frames;
    int queue_front[2], queue_scan[2];
    int priority_first, priority_count;
    uint64_t stamp;
    NsatSearchStats search;
};

static int variable(int lit)
{
    return lit > 0 ? lit : -lit;
}

static size_t lit_index(int lit)
{
    return 2 * (size_t)variable(lit) + (lit < 0);
}

static int lit_value(const NsatSolver *s, int lit)
{
    return s->values[lit_index(lit)];
}

static void set_value(NsatSolver *s, int lit, int value)
{
    size_t i = lit_index(lit);
    s->values[i] = (signed char)value;
    s->values[i ^ 1] = (signed char)-value;
}

static void *reserve(NsatBudget *b, void *old, size_t *capacity, size_t used,
                     size_t additional, size_t width)
{
    if (additional <= *capacity - used) return old;
    if (additional > SIZE_MAX - used) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
        return NULL;
    }
    size_t required = used + additional;
    size_t next = *capacity > SIZE_MAX / 2 ? required : *capacity * 2;
    if (next < 16) next = 16;
    if (next < required) next = required;
    void *buffer = nsat_alloc(b, next, width);
    if (!buffer) return NULL;
    if (!nsat_copy(b, buffer, old, used * width)) {
        nsat_free(b, buffer, next * width);
        return NULL;
    }
    nsat_free(b, old, *capacity * width);
    *capacity = next;
    return buffer;
}

NsatStatus nsat_init(NsatSolver **out, int variables, NsatLimits limits)
{
    if (!out) return NSAT_ERROR;
    *out = NULL;
    size_t count;
    if (variables < 0 || !nsat_mul((size_t)variables + 1, 2, &count)) return NSAT_ERROR;
    NsatBudget b = {.limits = limits};
    NsatSolver *s = nsat_alloc(&b, 1, sizeof(*s));
    if (!s) return b.status;
    *s = (NsatSolver){.budget = b, .variables = variables};
    *out = s;
    s->values = nsat_alloc(&s->budget, count, sizeof(*s->values));
    if (s->values) nsat_copy(&s->budget, s->values, NULL, count);
    return s->budget.status;
}

static NsatStatus store_clause(NsatSolver *s, NsatClauseId id,
                               const int *lits, size_t count)
{
    NsatBudget *b = &s->budget;
    /* A watch link encodes the clause number, side and binary tag. */
    if (s->count >= SIZE_MAX / 4) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
        return b->status;
    }
    SolverClause *clauses = reserve(b, s->clauses, &s->clause_capacity, s->count, 1,
                                  sizeof(*s->clauses));
    if (!clauses) return b->status;
    s->clauses = clauses;
    int *words = reserve(b, s->words, &s->word_capacity, s->word_count, count,
                         sizeof(*s->words));
    if (count && !words) return b->status;
    s->words = words;
    if (count && !nsat_copy(b, words + s->word_count, lits, count * sizeof(int)))
        return b->status;
    s->clauses[s->count++] = (SolverClause){.offset = s->word_count, .count = count, .id = id};
    s->word_count += count;
    if (count > s->max_width) s->max_width = count;
    return b->status;
}

NsatStatus nsat_add_original(NsatSolver *s, const int *lits, size_t count)
{
    if (!s) return NSAT_ERROR;
    NsatBudget *b = &s->budget;
    if (!nsat_poll(b)) return b->status;
    if (s->started || s->checker) nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
    if (count > (size_t)s->variables || (count && !lits))
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
    if (b->status != NSAT_OK) return b->status;
    for (size_t i = 0; i < count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int var = nsat_var(lits[i], s->variables);
        if (!var || lit_value(s, var)) {
            nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
            return b->status;
        }
        set_value(s, var, 1);
    }
    for (size_t i = 0; i < count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        set_value(s, variable(lits[i]), 0);
    }
    return store_clause(s, (NsatClauseId)s->count + 1, lits, count);
}

static NsatStatus checker_status(NsatSolver *s, NsatVerifier *v, NsatStatus status)
{
    if (status == NSAT_ERROR || status == NSAT_UNKNOWN)
        nsat_stop(&s->budget, status, nsat_verify_stats(v).stop);
    return status;
}

/* The reason contains the forced literal or is entirely false. Cite the
 * already checked root units, then the intact reason under its own ID. */
static NsatStatus root_event(NsatSolver *s, NsatVerifier *v, size_t index, int forced)
{
    NsatBudget *b = &s->budget;
    const SolverClause *reason = &s->clauses[index];
    NsatClauseId id = reason->id;
    if (forced && lit_value(s, forced) > 0) return NSAT_OK;
    if (!forced && !reason->count)
        return checker_status(s, v, nsat_verify_conclude_unsat(v, id));

    bool reuse = forced && reason->count == 1 && !lit_value(s, forced);
    if (!reuse) {
        size_t used = 0;
        int implied = forced && !lit_value(s, forced) ? forced : 0;
        for (size_t i = 0; i < reason->count; ++i) {
            if (!nsat_work(b, 1)) return b->status;
            int lit = s->words[reason->offset + i];
            if (lit == implied) continue;
            int var = lit > 0 ? lit : -lit;
            if (lit_value(s, lit) >= 0 || !s->units[var]) {
                nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
                return b->status;
            }
            s->hints[used++] = s->units[var];
        }
        s->hints[used++] = id;
        if (s->last_id == UINT64_MAX) {
            nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
            return b->status;
        }
        id = ++s->last_id;
        NsatStatus status = nsat_verify_derive(v, id, &implied, implied != 0,
                                               s->hints, used);
        if (status != NSAT_OK) return checker_status(s, v, status);
        if (!implied)
            return checker_status(s, v, nsat_verify_conclude_unsat(v, id));
    } else {
        const int *lits;
        size_t count;
        if (!nsat_verify_clause(v, id, &lits, &count) || count != 1 || lits[0] != forced) {
            nsat_stop(b, NSAT_ERROR, NSAT_STOP_PROOF);
            return b->status;
        }
    }
    int var = forced > 0 ? forced : -forced;
    set_value(s, forced, 1);
    s->units[var] = id;
    s->trail[s->trail_size++] = forced;
    return NSAT_OK;
}


static NsatStatus solver_error(NsatSolver *s, NsatStopReason reason)
{
    nsat_stop(&s->budget, NSAT_ERROR, reason);
    return s->budget.status;
}

static NsatStatus bind_checker(NsatSolver *s, NsatVerifier *v)
{
    if (!v || s->started) return solver_error(s, NSAT_STOP_STATE);
    if (s->checker) return s->checker == v ? NSAT_OK : solver_error(s, NSAT_STOP_STATE);
    if (nsat_verify_variables(v) != s->variables || nsat_verify_original_count(v) != s->count)
        return solver_error(s, NSAT_STOP_STATE);
    NsatStatus ready = nsat_verify_begin(v);
    if (ready != NSAT_OK) return checker_status(s, v, ready);
    s->checker = v;
    s->last_id = s->count;
    return NSAT_OK;
}

NsatStatus nsat_add_root_rup(NsatSolver *s, NsatVerifier *v,
                             const int *lits, size_t count,
                             const NsatClauseId *hints, size_t hint_count,
                             NsatClauseId *id)
{
    if (id) *id = 0;
    if (!s) return NSAT_ERROR;
    if (!nsat_poll(&s->budget)) return s->budget.status;
    NsatStatus ready = bind_checker(s, v);
    if (ready != NSAT_OK) return ready;
    if (s->last_id == UINT64_MAX) return solver_error(s, NSAT_STOP_STATE);
    NsatClauseId next = s->last_id + 1;
    NsatStatus status = nsat_verify_derive(v, next, lits, count, hints, hint_count);
    if (status != NSAT_OK) return checker_status(s, v, status);
    status = store_clause(s, next, lits, count);
    if (status != NSAT_OK) return status;
    s->last_id = next;
    if (id) *id = next;
    return NSAT_OK;
}

static void link_watch(NsatSolver *s, size_t index, size_t side, int lit)
{
    size_t *head = &s->heads[lit_index(lit)];
    SolverClause *clause = &s->clauses[index];
    bool binary = clause->count == 2;
    if (binary) clause->watch[side].mate = s->words[clause->offset + (side ^ 1)];
    clause->watch[side].next = *head;
    *head = 2 * (2 * index + side + 1) + binary;
}

static void watch_clause(NsatSolver *s, size_t index)
{
    const SolverClause *clause = &s->clauses[index];
    link_watch(s, index, 0, s->words[clause->offset]);
    link_watch(s, index, 1, s->words[clause->offset + 1]);
}

static NsatStatus enqueue(NsatSolver *s, NsatVerifier *v, int lit, size_t ref)
{
    if (!s->level) return root_event(s, v, ref - 1, lit);
    Variable *var = &s->vars[variable(lit)];
    var->level = s->level;
    var->reason = ref;
    var->phase = lit > 0 ? 1 : -1;
    set_value(s, lit, 1);
    s->trail[s->trail_size++] = lit;
    return NSAT_OK;
}

/* Store conflicts by clause reference. Binary propagation reads the
 * other literal from the watch and leaves the clause arena untouched. */
static NsatStatus propagate_false(NsatSolver *s, NsatVerifier *v, int lit,
                                  size_t *conflict)
{
    NsatBudget *b = &s->budget;
    size_t *link = &s->heads[lit_index(lit)];
    while (*link) {
        if (!nsat_work(b, 1)) return b->status;
        size_t encoded = *link;
        size_t watch = (encoded >> 1) - 1;
        size_t index = watch / 2, side = watch % 2;
        SolverClause *clause = &s->clauses[index];
        Watch *entry = &clause->watch[side];
        const bool binary = encoded & 1;
        int *lits = binary ? NULL : s->words + clause->offset;
        int other = binary ? entry->mate : lits[side ^ 1];
        int value = lit_value(s, other);
        if (value > 0) {
            link = &entry->next;
            continue;
        }
        if (!binary) {
            size_t i = 2;
            while (i < clause->count) {
                if (!nsat_work(b, 1)) return b->status;
                if (lit_value(s, lits[i]) >= 0) break;
                ++i;
            }
            if (i < clause->count) {
                *link = entry->next;
                lits[side] = lits[i];
                lits[i] = lit;
                link_watch(s, index, side, lits[side]);
                continue;
            }
        }
        if (value < 0) {
            *conflict = index + 1;
            return NSAT_OK;
        }
        if (binary) ++s->search.binary_propagations;
        NsatStatus status = enqueue(s, v, other, index + 1);
        if (status != NSAT_OK) return status;
        link = &entry->next;
    }
    return NSAT_OK;
}

static NsatStatus propagate(NsatSolver *s, NsatVerifier *v, size_t *conflict)
{
    *conflict = 0;
    while (s->propagated < s->trail_size) {
        if (!nsat_work(&s->budget, 1)) return s->budget.status;
        NsatStatus status = propagate_false(s, v, -s->trail[s->propagated++],
                                            conflict);
        if (status != NSAT_OK) return status;
        if (*conflict) return s->level ? NSAT_OK : root_event(s, v, *conflict - 1, 0);
    }
    return NSAT_OK;
}

static int tier(const NsatSolver *s, int variable)
{
    return s->priority_count && variable >= s->priority_first &&
        variable - s->priority_first < s->priority_count;
}

static NsatStatus prepare_search(NsatSolver *s, NsatSearchOptions options)
{
    NsatBudget *b = &s->budget;
    size_t slots = (size_t)s->variables + 1;
    s->vars = nsat_alloc(b, slots, sizeof(*s->vars));
    if (!s->vars || !nsat_copy(b, s->vars, NULL, slots * sizeof(*s->vars)))
        return b->status;
    s->level_start = nsat_alloc(b, slots, sizeof(*s->level_start));
    s->grade = nsat_alloc(b, slots, sizeof(*s->grade));
    s->learned = nsat_alloc(b, slots, sizeof(*s->learned));
    s->touched = nsat_alloc(b, slots, sizeof(*s->touched));
    s->frames = nsat_alloc(b, slots, sizeof(*s->frames));
    if (b->status != NSAT_OK ||
        !nsat_copy(b, s->grade, NULL, slots * sizeof(*s->grade))) return b->status;

    if (options.priority_count < 0 ||
        (options.priority_count && (options.priority_first < 1 ||
         options.priority_first > s->variables ||
         options.priority_count - 1 > s->variables - options.priority_first)))
        return solver_error(s, NSAT_STOP_INPUT);
    s->priority_first = options.priority_first;
    s->priority_count = options.priority_count;
    int tail[2] = {0};
    for (size_t i = 1; i < slots; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int group = tier(s, (int)i);
        s->vars[i].prev = tail[group];
        if (tail[group]) s->vars[tail[group]].next = (int)i;
        else s->queue_front[group] = s->queue_scan[group] = (int)i;
        tail[group] = (int)i;
        s->vars[i].stamp = slots - i;
        s->vars[i].phase = group ? 1 : -1;
    }
    s->stamp = (uint64_t)s->variables;
    return NSAT_OK;
}

static NsatStatus start_root(NsatSolver *s, NsatVerifier *v, bool search)
{
    if (!s) return NSAT_ERROR;
    NsatBudget *b = &s->budget;
    if (!nsat_poll(b)) return b->status;
    NsatStatus ready = bind_checker(s, v);
    if (ready != NSAT_OK) return ready;
    s->started = true;
    s->pinned_count = s->count;
    size_t slots = (size_t)s->variables + 1;
    size_t heads = 2 * slots; /* Checked by nsat_init. */
    s->heads = nsat_alloc(b, heads, sizeof(*s->heads));
    if (!s->heads || !nsat_copy(b, s->heads, NULL, heads * sizeof(*s->heads)))
        return b->status;
    s->units = nsat_alloc(b, slots, sizeof(*s->units));
    s->trail = nsat_alloc(b, slots, sizeof(*s->trail));
    s->hints = nsat_alloc(b, search ? slots : s->max_width + 1, sizeof(*s->hints));
    if (b->status != NSAT_OK) return b->status;

    for (size_t i = 0; i < s->count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        if (s->clauses[i].count >= 2) watch_clause(s, i);
    }
    for (size_t i = 0; i < s->count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        const SolverClause *clause = &s->clauses[i];
        if (clause->count <= 1) {
            int forced = clause->count ? s->words[clause->offset] : 0;
            NsatStatus status = root_event(s, v, i, forced);
            if (status != NSAT_OK) return status;
        }
    }
    size_t conflict;
    NsatStatus status = propagate(s, v, &conflict);
    if (status != NSAT_OK) return status;
    return nsat_poll(b) ? NSAT_ROOT_FIXED_POINT : b->status;
}

NsatStatus nsat_propagate_root(NsatSolver *s, NsatVerifier *v)
{
    return start_root(s, v, false);
}

static bool bump(NsatSolver *s, int i)
{
    Variable *var = &s->vars[i];
    if (!var->prev) return true;
    if (s->stamp == UINT64_MAX)
        return nsat_stop(&s->budget, NSAT_ERROR, NSAT_STOP_STATE);
    s->vars[var->prev].next = var->next;
    if (var->next) s->vars[var->next].prev = var->prev;
    var->prev = 0;
    int group = tier(s, i);
    var->next = s->queue_front[group];
    s->vars[var->next].prev = i;
    s->queue_front[group] = s->queue_scan[group] = i;
    var->stamp = ++s->stamp;
    return true;
}

static void mark(NsatSolver *s, int i, unsigned char flag)
{
    Variable *var = &s->vars[i];
    if (!var->marks) s->touched[s->touched_count++] = i;
    var->marks |= flag;
}

static bool clear_marks(NsatSolver *s)
{
    while (s->touched_count) {
        if (!nsat_work(&s->budget, 1)) return false;
        s->vars[s->touched[--s->touched_count]].marks = 0;
    }
    return true;
}

static bool first_uip(NsatSolver *s, size_t ref)
{
    NsatBudget *b = &s->budget;
    size_t pending = 0, cursor = s->trail_size;
    int pivot = 0;
    s->learned_count = 1;

    for (;;) {
        SolverClause *clause = &s->clauses[ref - 1];
        clause->used = true;
        for (size_t i = 0; i < clause->count; ++i) {
            if (!nsat_work(b, 1)) return false;
            int lit = s->words[clause->offset + i];
            int var = variable(lit);
            Variable *state = &s->vars[var];
            if (lit == pivot || !state->level || (state->marks & ANALYSED)) continue;
            mark(s, var, ANALYSED);
            if (!bump(s, var)) return false;
            if (state->level == s->level) {
                ++pending;
            } else {
                mark(s, var, SOURCE);
                s->learned[s->learned_count++] = lit;
            }
        }

        /* Resolve the last marked assignment at this level. Clauses only
         * refer to earlier assignments, so the trail cursor moves one way. */
        if (!pending) return nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
        for (;;) {
            if (!nsat_work(b, 1)) return false;
            if (cursor == s->level_start[s->level])
                return nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
            pivot = s->trail[--cursor];
            if (s->vars[variable(pivot)].marks & ANALYSED) break;
        }
        if (--pending == 0) {
            s->learned[0] = -pivot;
            mark(s, variable(pivot), SOURCE);
            return true;
        }
        ref = s->vars[variable(pivot)].reason;
        if (!ref) return nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
    }
}

/* A reason can justify removing a literal if each antecedent is at root,
 * remains in the learned clause, or is itself recursively justified.
 * The iterative traversal leaves successful dependency marks in place
 * when a sibling fails. A later pass reconstructs the actual proof. */
static bool removable(NsatSolver *s, int candidate)
{
    size_t depth = 1;
    s->frames[0] = (ReasonFrame){.variable = candidate};
    while (depth) {
        if (!nsat_work(&s->budget, 1)) return false;
        ReasonFrame *frame = &s->frames[depth - 1];
        int var = frame->variable;
        size_t ref = s->vars[var].reason;
        if (!ref) break;
        SolverClause *reason = &s->clauses[ref - 1];
        reason->used = true;
        if (frame->next_literal == reason->count) {
            mark(s, var, REMOVABLE);
            --depth;
            continue;
        }
        int lit = s->words[reason->offset + frame->next_literal++];
        int dep = variable(lit);
        if (dep == var || !s->vars[dep].level) continue;
        unsigned char flags = s->vars[dep].marks;
        if (flags & (SOURCE | REMOVABLE)) continue;
        if ((flags & BLOCKED) || !s->vars[dep].reason) break;
        s->frames[depth++] = (ReasonFrame){.variable = dep};
    }
    if (!depth) return true;
    while (depth) {
        if (!nsat_work(&s->budget, 1)) return false;
        mark(s, s->frames[--depth].variable, BLOCKED);
    }
    return false;
}

static bool minimise(NsatSolver *s)
{
    size_t kept = 1;
    for (size_t i = 1; i < s->learned_count; ++i) {
        if (!nsat_work(&s->budget, 1)) return false;
        int lit = s->learned[i], var = variable(lit);
        if (removable(s, var)) {
            s->vars[var].marks &= (unsigned char)~SOURCE;
            ++s->search.minimised_literals;
        } else {
            if (s->budget.status != NSAT_OK) return false;
            s->learned[kept++] = lit;
        }
    }
    s->learned_count = kept;
    return true;
}

/* Replay starts with the negation of the final lemma (SOURCE).
 * Walk the reason graph in postorder. Each antecedent precedes its
 * dependent on the assignment trail, so this graph is acyclic.
 * The conflict is the synthetic first frame and the final hint. */
static bool proof_chain(NsatSolver *s, size_t conflict, size_t *length)
{
    NsatBudget *b = &s->budget;
    size_t depth = 1, count = 0;
    s->frames[0] = (ReasonFrame){0};
    while (depth) {
        if (!nsat_work(b, 1)) return false;
        ReasonFrame *frame = &s->frames[depth - 1];
        int var = frame->variable;
        if (var && !s->vars[var].level) {
            if (!s->units[var]) return nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
            s->hints[count++] = s->units[var];
            --depth;
            continue;
        }
        size_t ref = var ? s->vars[var].reason : conflict;
        if (!ref) return nsat_stop(b, NSAT_ERROR, NSAT_STOP_STATE);
        const SolverClause *reason = &s->clauses[ref - 1];
        if (frame->next_literal == reason->count) {
            s->hints[count++] = reason->id;
            --depth;
            continue;
        }
        int dep = variable(s->words[reason->offset + frame->next_literal++]);
        if (dep == var || (s->vars[dep].marks & (SOURCE | NEEDED))) continue;
        mark(s, dep, NEEDED);
        s->frames[depth++] = (ReasonFrame){.variable = dep};
    }
    *length = count;
    return true;
}

static bool backjump(NsatSolver *s, int to)
{
    size_t retained = s->level_start[to + 1];
    while (s->trail_size > retained) {
        if (!nsat_work(&s->budget, 1)) return false;
        int lit = s->trail[--s->trail_size], var = variable(lit);
        set_value(s, lit, 0);
        s->vars[var].level = 0;
        s->vars[var].reason = 0;
        int *scan = &s->queue_scan[tier(s, var)];
        if (!*scan || s->vars[var].stamp > s->vars[*scan].stamp) *scan = var;
    }
    s->level = to;
    s->propagated = retained;
    return true;
}

static NsatStatus learn(NsatSolver *s, NsatVerifier *v, size_t conflict,
                         bool minimise_literals)
{
    NsatBudget *b = &s->budget;
    size_t hints;
    if (!first_uip(s, conflict) ||
        (minimise_literals && !minimise(s)) ||
        !proof_chain(s, conflict, &hints) || !clear_marks(s))
        return b->status;

    int jump_to = 0;
    size_t second_watch = 1;
    for (size_t i = 1; i < s->learned_count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int level = s->vars[variable(s->learned[i])].level;
        if (level > jump_to) {
            jump_to = level;
            second_watch = i;
        }
    }
    if (s->learned_count > 1) {
        int other = s->learned[second_watch];
        s->learned[second_watch] = s->learned[1];
        s->learned[1] = other;
    }
    if (s->last_id == UINT64_MAX) return solver_error(s, NSAT_STOP_STATE);
    NsatClauseId id = ++s->last_id;
    unsigned glue = 0;
    for (size_t i = 0; i < s->learned_count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int level = s->vars[variable(s->learned[i])].level;
        if (s->grade[level] != id) {
            s->grade[level] = id;
            ++glue;
        }
    }
    NsatStatus status = nsat_verify_derive(v, id, s->learned, s->learned_count,
                                          s->hints, hints);
    if (status != NSAT_OK) return checker_status(s, v, status);
    status = store_clause(s, id, s->learned, s->learned_count);
    if (status != NSAT_OK) return status;
    s->clauses[s->count - 1].glue = glue;
    bool nonchronological = jump_to < s->level - 1;
    if (!backjump(s, jump_to)) return b->status;
    if (s->learned_count > 1) watch_clause(s, s->count - 1);
    ++s->search.learned;
    if (nonchronological) ++s->search.nonchronological_backjumps;
    return enqueue(s, v, s->learned[0], s->count);
}

/* At root there are no live arena reason references. Pinned root units
 * use proof IDs, and a unit in the arena is always retained. Move survivors
 * down in place, then rebuild links from their new clause references. */
static NsatStatus reduce_at_root(NsatSolver *s, NsatVerifier *v)
{
    NsatBudget *b = &s->budget;
    if (s->level || s->propagated != s->trail_size)
        return solver_error(s, NSAT_STOP_STATE);
    for (size_t i = 0; i < s->trail_size; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        int var = variable(s->trail[i]);
        if (s->vars[var].reason || s->vars[var].level || !s->units[var])
            return solver_error(s, NSAT_STOP_STATE);
    }

    size_t kept = s->pinned_count;
    size_t words = kept < s->count ? s->clauses[kept].offset : s->word_count;
    uint64_t deleted = 0, moved = 0;
    for (size_t i = kept; i < s->count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        SolverClause clause = s->clauses[i];
        if (clause.count > 2 && clause.glue > 2 && !clause.used) {
            NsatStatus status = nsat_verify_delete(v, clause.id);
            if (status != NSAT_OK) return checker_status(s, v, status);
            ++deleted;
            continue;
        }
        clause.used = false;
        if (words != clause.offset) {
            /* Destination precedes source, so copying forward in chunks
             * preserves overlapping ranges and polls during long moves. */
            for (size_t j = 0; j < clause.count;) {
                size_t chunk = clause.count - j;
                if (chunk > 4096 / sizeof(int)) chunk = 4096 / sizeof(int);
                if (!nsat_work(b, chunk * sizeof(int))) return b->status;
                memmove(s->words + words + j, s->words + clause.offset + j,
                        chunk * sizeof(int));
                j += chunk;
            }
            ++moved;
        }
        clause.offset = words;
        s->clauses[kept++] = clause;
        words += clause.count;
    }
    ++s->search.reductions;
    if (!deleted) return NSAT_OK;
    ++s->search.compactions;
    s->search.deleted += deleted;
    s->search.moved += moved;
    s->search.words_reclaimed += s->word_count - words;
    s->count = kept;
    s->word_count = words;

    size_t heads = 2 * ((size_t)s->variables + 1);
    if (!nsat_copy(b, s->heads, NULL, heads * sizeof(*s->heads))) return b->status;
    for (size_t i = 0; i < s->count; ++i) {
        if (!nsat_work(b, 1)) return b->status;
        if (s->clauses[i].count >= 2) watch_clause(s, i);
    }
    return NSAT_OK;
}

typedef struct {
    uint64_t block, scale, next;
} RestartSchedule;

/* The low bit of block is the final term of that block of the Luby sequence. */
static void advance_restart(RestartSchedule *r, uint64_t base, uint64_t conflicts)
{
    if ((r->block & (~r->block + 1)) == r->scale) {
        if (r->block != UINT64_MAX) ++r->block;
        r->scale = 1;
    } else {
        r->scale *= 2;
    }
    r->next = r->scale > (UINT64_MAX - conflicts) / base
        ? UINT64_MAX : conflicts + base * r->scale;
}

static NsatStatus search(NsatSolver *s, NsatVerifier *v)
{
    NsatBudget *b = &s->budget;
    NsatSearchOptions options = s->options;
    uint64_t base = options.restart_base ? options.restart_base : 64;
    uint64_t start = s->search.conflicts;
    RestartSchedule restart = {.block = 1, .scale = 1,
        .next = base > UINT64_MAX - start ? UINT64_MAX : start + base};
    uint64_t interval = options.reduction_interval ? options.reduction_interval : 2000;
    uint64_t last_reduction = start;
    for (;;) {
        size_t conflict;
        NsatStatus status = propagate(s, v, &conflict);
        if (status != NSAT_OK) return status;
        if (conflict) {
            if (b->limits.conflicts && s->search.conflicts >= b->limits.conflicts) {
                nsat_stop(b, NSAT_UNKNOWN, NSAT_STOP_CONFLICTS);
                return b->status;
            }
            ++s->search.conflicts;
            status = learn(s, v, conflict, !options.disable_minimisation);
            if (status != NSAT_OK) return status;
            continue;
        }

        /* A complete assignment is checked before considering a restart. */
        if (s->trail_size == (size_t)s->variables) {
            if (!nsat_poll(b)) return b->status;
            status = nsat_verify_conclude_sat(v, s->trail, s->trail_size);
            s->model_valid = status == NSAT_SAT;
            return checker_status(s, v, status);
        }
        if (s->search.conflicts >= restart.next) {
            if (s->level) {
                if (!backjump(s, 0)) return b->status;
                ++s->search.restarts;
            }
            advance_restart(&restart, base, s->search.conflicts);
        }
        if (!options.disable_reduction && !s->level &&
            s->search.conflicts - last_reduction >= interval) {
            status = reduce_at_root(s, v);
            if (status != NSAT_OK) return status;
            last_reduction = s->search.conflicts;
        }
        int next = 0;
        for (int group = 1; group >= 0 && !next; --group) {
            int *scan = &s->queue_scan[group];
            while (*scan && lit_value(s, *scan)) {
                if (!nsat_work(b, 1)) return b->status;
                *scan = s->vars[*scan].next;
            }
            next = *scan;
        }
        if (!next) return solver_error(s, NSAT_STOP_STATE);
        if (!nsat_work(b, 1)) return b->status;
        s->level_start[++s->level] = s->trail_size;
        ++s->search.decisions;
        enqueue(s, v, s->vars[next].phase * next, 0);
    }
}

NsatStatus nsat_solve(NsatSolver *s, NsatVerifier *v, NsatSearchOptions options)
{
    NsatStatus status = start_root(s, v, true);
    if (status != NSAT_ROOT_FIXED_POINT) return status;
    status = prepare_search(s, options);
    if (status != NSAT_OK) return status;
    s->options = options;
    return search(s, v);
}

NsatStatus nsat_solve_tighter_cover(NsatSolver *s, NsatVerifier *v, size_t target)
{
    if (!s) return NSAT_ERROR;
    NsatBudget *b = &s->budget;
    if (!nsat_poll(b)) return b->status;
    if (!s->model_valid || !v || s->checker != v)
        return solver_error(s, NSAT_STOP_STATE);
    s->model_valid = false;
    NsatClauseId id;
    int literal;
    NsatStatus status = nsat_verify_tighten_cover(v, target, &id, &literal);
    if (status != NSAT_OK) return checker_status(s, v, status);
    if (id <= s->last_id) return solver_error(s, NSAT_STOP_STATE);
    s->last_id = id;
    if (s->level && !backjump(s, 0)) return b->status;
    /* Units survive reduction and need no watches. Existing root
     * assignments either satisfy this unit, conflict, or propagate it. */
    status = store_clause(s, id, &literal, 1);
    if (status != NSAT_OK) return status;
    status = root_event(s, v, s->count - 1, literal);
    return status == NSAT_OK ? search(s, v) : status;
}

int nsat_model_value(const NsatSolver *s, int variable)
{
    return s->model_valid && s->budget.status == NSAT_OK &&
        variable > 0 && variable <= s->variables ? lit_value(s, variable) : 0;
}

int nsat_root_value(const NsatSolver *s, int variable)
{
    if (s->budget.status != NSAT_OK || variable <= 0 || variable > s->variables ||
        (s->vars && s->vars[variable].level))
        return 0;
    return lit_value(s, variable);
}

NsatClauseId nsat_root_unit_id(const NsatSolver *s, int variable)
{
    return nsat_root_value(s, variable) ? s->units[variable] : 0;
}

NsatStats nsat_stats(const NsatSolver *s) { return s->budget.stats; }
NsatSearchStats nsat_search_stats(const NsatSolver *s) { return s->search; }

void nsat_destroy(NsatSolver *s)
{
    if (!s) return;
    free(s->clauses);
    free(s->words);
    free(s->values);
    free(s->units);
    free(s->hints);
    free(s->trail);
    free(s->heads);
    free(s->vars);
    free(s->level_start);
    free(s->grade);
    free(s->learned);
    free(s->touched);
    free(s->frames);
    nsat_release_budget(&s->budget);
    free(s);
}
