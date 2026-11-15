/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "cover_encoding.h"
#include "native_sat_support.h"

enum { MAX_ROWS = 96, MAX_OFF = 128, MAX_LEMMAS = 32768 };
static const uint64_t max_intersections = UINT64_C(2000000);

typedef struct {
    unsigned char rows[3], width, off;
} ForbiddenGroup;

typedef struct {
    NsatBudget budget;
    NsatCoverDomainStats *stats;
    ForbiddenGroup *groups;
    size_t capacity, pairs, triples, triples_seen;
    unsigned random;
    uint64_t intersections;
    size_t order[MAX_ROWS];
} Domain;

static bool prepare_masks(const NsatCoverInput *q, Domain *d, uint64_t *masks, size_t words)
{
    for (size_t r = 0; r < q->positives; ++r) {
        for (size_t z = 0; z < q->negatives; ++z) {
            const int *on = q->on + r * q->inputs, *off = q->off + z * q->inputs;
            for (size_t word = 0, i = 0; word < words; ++word) {
                uint64_t mask = 0;
                for (unsigned bit = 0; bit < 64 && i < q->inputs; ++bit, ++i) {
                    if (!nsat_work(&d->budget, 1)) return false;
                    if (on[i] != off[i]) mask |= UINT64_C(1) << bit;
                }
                *masks++ = mask;
            }
        }
    }
    return true;
}

static bool charge_intersection(Domain *d, size_t words)
{
    if (words > max_intersections - d->intersections) {
        d->stats->discovery_capped = true;
        return false;
    }
    d->intersections += words;
    return nsat_work(&d->budget, words);
}

static bool has_common_difference(const uint64_t *a, const uint64_t *b,
                                  const uint64_t *c, size_t words)
{
    for (size_t word = 0; word < words; ++word) {
        uint64_t common = a[word] & b[word];
        if (c) common &= c[word];
        if (common) return true;
    }
    return false;
}

static void save_triple(Domain *d, ForbiddenGroup group)
{
    size_t capacity = d->capacity - d->pairs;
    ++d->triples_seen;
    if (d->triples < capacity) {
        d->groups[d->pairs + d->triples++] = group;
    } else {
        d->stats->discovery_capped = true;
        /* Deterministic reservoir to cover rows beyond the lexicographic prefix. */
        d->random ^= d->random << 13;
        d->random ^= d->random >> 17;
        d->random ^= d->random << 5;
        size_t i = d->random % d->triples_seen;
        if (i < capacity) d->groups[d->pairs + i] = group;
    }
}

static bool precedes(int a, int b, const int *degree, const int *triangles)
{
    if (b < 0) return true;
    if (degree[a] != degree[b]) return degree[a] > degree[b];
    return triangles[a] > triangles[b];
}

static void choose_order(Domain *d, size_t p, const uint64_t graph[][2],
                         const int *degree, const int *triangles)
{
    bool taken[MAX_ROWS] = {0};
    uint64_t remaining[2] = {UINT64_MAX, UINT64_MAX};
    size_t anchors = 0;
    for (;;) {
        int best = -1;
        for (size_t r = 0; r < p; ++r) {
            if (!nsat_work(&d->budget, 1)) return;
            if ((remaining[r / 64] & (UINT64_C(1) << (r % 64))) &&
                precedes((int)r, best, degree, triangles))
                best = (int)r;
        }
        if (best < 0) break;
        taken[best] = true;
        d->order[anchors++] = (size_t)best;
        remaining[0] &= graph[best][0];
        remaining[1] &= graph[best][1];
    }
    d->stats->anchors = anchors;
    for (size_t at = anchors; at < p; ++at) {
        int best = -1;
        for (size_t r = 0; r < p; ++r) {
            if (!nsat_work(&d->budget, 1)) return;
            if (!taken[r] && precedes((int)r, best, degree, triangles)) best = (int)r;
        }
        taken[best] = true;
        d->order[at] = (size_t)best;
    }
}

static void discover(const NsatCoverInput *q, Domain *d, bool groups, bool ordering)
{
    const size_t p = q->positives, off = q->negatives;
    const size_t words = (q->inputs - 1) / 64 + 1;
    size_t cells;
    if (!nsat_mul(p * off, words, &cells)) {
        nsat_stop(&d->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        return;
    }
    uint64_t *masks = nsat_alloc(&d->budget, cells, sizeof(*masks));
    if (!masks) return;
    uint64_t graph[MAX_ROWS][2] = {{0}};
    int degree[MAX_ROWS] = {0}, triangles[MAX_ROWS] = {0};
    if (!prepare_masks(q, d, masks, words)) goto done;
    d->random = 20261114u;
    for (size_t a = 0; a < p; ++a) {
        for (size_t b = a + 1; b < p; ++b) {
            for (size_t z = 0; z < off; ++z) {
                if (!charge_intersection(d, words)) goto ordered;
                const uint64_t *first = masks + (a * off + z) * words;
                const uint64_t *second = masks + (b * off + z) * words;
                if (has_common_difference(first, second, NULL, words)) continue;
                graph[a][b / 64] |= UINT64_C(1) << (b % 64);
                graph[b][a / 64] |= UINT64_C(1) << (a % 64);
                ++degree[a];
                ++degree[b];
                if (groups) {
                    if (d->pairs < d->capacity)
                        d->groups[d->pairs++] = (ForbiddenGroup){{a, b}, 2, z};
                    else d->stats->discovery_capped = true;
                }
                break;
            }
        }
    }
    for (size_t a = 0; a < p; ++a) {
        for (size_t b = a + 1; b < p; ++b) {
            if (graph[a][b / 64] & (UINT64_C(1) << (b % 64))) continue;
            for (size_t c = b + 1; c < p; ++c) {
                if (!nsat_work(&d->budget, 1)) goto done;
                if ((graph[a][c / 64] | graph[b][c / 64]) & (UINT64_C(1) << (c % 64)))
                    continue; /* A pair already excludes this triple. */
                for (size_t z = 0; z < off; ++z) {
                    if (!charge_intersection(d, words)) goto ordered;
                    const uint64_t *first = masks + (a * off + z) * words;
                    const uint64_t *second = masks + (b * off + z) * words;
                    const uint64_t *third = masks + (c * off + z) * words;
                    if (has_common_difference(first, second, third, words)) continue;
                    ++triangles[a];
                    ++triangles[b];
                    ++triangles[c];
                    if (groups) save_triple(d, (ForbiddenGroup){{a, b, c}, 3, z});
                    break;
                }
            }
        }
    }
ordered:
    if (ordering && d->budget.status == NSAT_OK) choose_order(d, p, graph, degree, triangles);
done:
    nsat_free(&d->budget, masks, cells * sizeof(*masks));
}

static size_t source_row(const NsatCoverInput *q, size_t r)
{
    return q->row_order ? q->row_order[r] : r;
}

static int assigned(const NsatCoverInput *q, size_t r, size_t t)
{
    return (int)(2 * q->inputs * q->target + r * q->target + t + 1);
}

static int used(const NsatCoverInput *q, size_t r, size_t t)
{
    return assigned(q, r, t) + (int)(q->positives * q->target);
}

/* Position of the first binary selector clause for this row/slot in v1. */
static NsatClauseId row_term_first(const NsatCoverInput *q, size_t r, size_t t)
{
    size_t first = r ? 1 + q->target * (q->inputs + 3 + (r - 1) * (q->inputs + 4)) : 1;
    return first + 1 + t * (q->inputs + 3 + (r != 0)) - (t != 0);
}

static NsatClauseId excludes_off(const NsatCoverInput *q, size_t t, size_t z)
{
    size_t on_clauses = q->target * (q->inputs + 3 + (q->positives - 1) * (q->inputs + 4));
    return on_clauses + t * (q->inputs + q->negatives) + q->inputs + z + 1;
}

static NsatStatus add_groups(const NsatCoverInput *q, Domain *d, NsatSolver *s, NsatVerifier *v)
{
    if (!d->pairs && !d->triples) return NSAT_OK;
    const size_t hint_count = q->inputs + 1;
    NsatClauseId *hints = nsat_alloc(&d->budget, hint_count, sizeof(*hints));
    if (!hints) return d->budget.status;
    NsatStatus status = NSAT_OK;
    size_t position[MAX_ROWS];
    for (size_t r = 0; r < q->positives; ++r) position[source_row(q, r)] = r;
    for (size_t g = 0; g < d->pairs + d->triples; ++g) {
        ForbiddenGroup group = d->groups[g];
        const int *off = q->off + group.off * q->inputs;
        size_t members[3];
        for (size_t j = 0; j < group.width; ++j) members[j] = position[group.rows[j]];
        for (size_t t = 0; t < q->target; ++t) {
            int clause[3];
            for (size_t j = 0; j < group.width; ++j) clause[j] = -assigned(q, members[j], t);
            for (size_t i = 0; i < q->inputs; ++i) {
                if (!nsat_work(&d->budget, 1)) {
                    status = d->budget.status;
                    goto done;
                }
                size_t j = 0;
                while (j < group.width && q->on[group.rows[j] * q->inputs + i] != off[i]) ++j;
                if (j == group.width) {
                    nsat_stop(&d->budget, NSAT_ERROR, NSAT_STOP_PROOF);
                    status = d->budget.status;
                    goto done;
                }
                hints[i] = row_term_first(q, members[j], t) + i;
            }
            hints[q->inputs] = excludes_off(q, t, group.off);
            status = nsat_add_root_rup(s, v, clause, group.width, hints, hint_count, NULL);
            if (status != NSAT_OK) goto done;
            if (group.width == 2) ++d->stats->pairs;
            else ++d->stats->triples;
        }
    }
done:
    nsat_free(&d->budget, hints, hint_count * sizeof(*hints));
    return status;
}

static NsatStatus add_staircase(const NsatCoverInput *q, Domain *d, NsatSolver *s, NsatVerifier *v)
{
    for (size_t t = 1; t < q->target; ++t) {
        NsatClauseId prior = 0;
        for (size_t r = 1; r < q->positives; ++r) {
            if (!nsat_work(&d->budget, 1)) return d->budget.status;
            NsatClauseId slot = row_term_first(q, r, t);
            NsatClauseId hints[] = {
                !prior ? row_term_first(q, 0, t) + q->inputs + 2
                       : row_term_first(q, r - 1, t - 1) + q->inputs + 1,
                !prior ? row_term_first(q, 0, t) + q->inputs + 1 : prior,
                slot + q->inputs + 2,
                slot + q->inputs + 3
            };
            int vertical[] = {-used(q, r, t), used(q, r - 1, t - 1)};
            NsatStatus status = nsat_add_root_rup(s, v, vertical, 2, hints, 4, &prior);
            if (status != NSAT_OK) return status;
            ++d->stats->staircase;
            int horizontal[] = {-used(q, r, t), used(q, r, t - 1)};
            NsatClauseId horizontal_hints[] = {row_term_first(q, r, t - 1) + q->inputs + 1, prior};
            status = nsat_add_root_rup(s, v, horizontal, 2, horizontal_hints, 2, NULL);
            if (status != NSAT_OK) return status;
            ++d->stats->staircase;
        }
    }
    return NSAT_OK;
}

NsatStatus nsat_cover_solve(const NsatCoverInput *input, NsatCoverStrengthening mode,
                            NsatLimits limits, NsatSolver **s, NsatVerifier **v,
                            NsatStopReason *stop, NsatCoverDomainStats *stats)
{
    *s = NULL;
    *v = NULL;
    *stop = NSAT_STOP_NONE;
    *stats = (NsatCoverDomainStats){0};
    NsatCoverInput q = *input;
    Domain d = {.stats = stats, .budget = {.limits = limits}};
    bool groups = mode == NSAT_COVER_GROUPS || mode == NSAT_COVER_COMBINED;
    bool ordering = mode == NSAT_COVER_ORDERING || mode == NSAT_COVER_COMBINED;
    const size_t max = (size_t)INT_MAX / 2;
    bool eligible = q.on && q.off && q.inputs &&
        q.positives && q.positives <= MAX_ROWS && q.negatives && q.negatives <= MAX_OFF &&
        q.target && q.target <= q.positives && !q.row_order && !q.order_count &&
        q.inputs <= max - q.positives && q.target <= max / (q.inputs + q.positives) &&
        q.inputs <= SIZE_MAX / sizeof(int) / q.positives &&
        q.inputs <= SIZE_MAX / sizeof(int) / q.negatives;
    if ((groups || ordering) && eligible) {
        if (groups) {
            size_t p = q.positives;
            size_t max_groups = p * (p - 1) / 2 + (p >= 3 ? p * (p - 1) * (p - 2) / 6 : 0);
            d.capacity = MAX_LEMMAS / q.target;
            if (d.capacity > max_groups) d.capacity = max_groups;
            if (d.capacity)
                d.groups = nsat_alloc(&d.budget, d.capacity, sizeof(*d.groups));
        }
        if (d.budget.status == NSAT_OK) discover(&q, &d, groups, ordering);
        if (ordering) {
            q.row_order = d.order;
            q.order_count = q.positives;
        }
    }
    NsatStatus status = d.budget.status;
    if (status == NSAT_OK) status = nsat_verify_init_cover(v, q.on, q.off, q.inputs,
        q.positives, q.negatives, q.target, q.row_order, q.order_count, limits);
    if (status == NSAT_OK) status = nsat_cover_load_solver(s, &q, *v, limits, stop);
    if (status == NSAT_OK && ordering && eligible) status = add_staircase(&q, &d, *s, *v);
    if (status == NSAT_OK && groups && eligible) status = add_groups(&q, &d, *s, *v);
    nsat_free(&d.budget, d.groups, d.capacity * sizeof(*d.groups));
    stats->prep_work = d.budget.stats.work;
    if (status == NSAT_OK) {
        NsatSearchOptions options = {0};
        if (ordering && eligible) {
            options.priority_first = assigned(&q, 0, 0);
            options.priority_count = (int)(q.positives * q.target);
        }
        status = nsat_solve(*s, *v, options);
    }
    if (!*stop) {
        if (d.budget.stats.stop) *stop = d.budget.stats.stop;
        else if (*s && nsat_stats(*s).stop) *stop = nsat_stats(*s).stop;
        else if (*v) *stop = nsat_verify_stats(*v).stop;
    }
    return status;
}
