/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "cover_encoding.h"
#include "native_sat_support.h"

enum { MAX_ROWS = 96, MAX_OFF = 128, MAX_LEMMAS = 32768 };
static const uint64_t max_intersections = UINT64_C(2000000);
static const uint64_t max_distance_comparisons = UINT64_C(2000000);

typedef struct {
    unsigned char rows[2], off;
} IncompatiblePair;

typedef struct {
    NsatBudget budget;
    NsatCoverDomainStats *stats;
    IncompatiblePair *pairs;
    size_t capacity, pair_count;
    uint64_t intersections;
    size_t order[MAX_ROWS];
    uint64_t (*graph)[2];
    bool active_pairs;
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

static bool precedes(int a, int b, const int *degree, const int *triangles)
{
    if (b < 0) return true;
    if (degree[a] != degree[b]) return degree[a] > degree[b];
    return triangles[a] > triangles[b];
}

static void choose_order(Domain *d, const NsatCoverInput *q, const uint64_t graph[][2],
                         const int *degree, const int *triangles)
{
    const size_t p = q->positives;
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
    /* Keep the clique unchanged. On a partial graph, missing edges leave
     * slots eligible rather than asserting compatibility or incompatibility. */
    bool rank_slots = anchors > 1 && anchors <= q->target && anchors < p;
    int compatible[MAX_ROWS];
    if (rank_slots) {
        for (size_t r = 0; r < p; ++r) {
            if (taken[r]) continue;
            compatible[r] = (int)anchors;
            for (size_t a = 0; a < anchors; ++a) {
                if (!nsat_work(&d->budget, 1)) return;
                if (graph[d->order[a]][r / 64] & (UINT64_C(1) << (r % 64)))
                    --compatible[r];
            }
        }
    }
    for (size_t at = anchors; at < p; ++at) {
        int best = -1;
        for (size_t r = 0; r < p; ++r) {
            if (!nsat_work(&d->budget, 1)) return;
            if (taken[r]) continue;
            if (best < 0) {
                best = (int)r;
            } else if (rank_slots && compatible[r] != compatible[best]) {
                if (compatible[r] < compatible[best]) best = (int)r;
            } else if (precedes((int)r, best, degree, triangles)) {
                best = (int)r;
            }
        }
        taken[best] = true;
        d->order[at] = (size_t)best;
    }
}

/* Earlier anchors occupy their own slots by induction. Only exclude a later
 * anchor from each earlier anchor's slot; clauses for other slots are not
 * needed. The first target+1 anchors suffice when the clique exceeds target.
 * Otherwise retain all known pairs for active-slot filtering. Reconstruct
 * witnesses after ordering. At most 96*95/2 witnesses are retained. */
static void retain_ordered_pairs(const NsatCoverInput *q, Domain *d,
                                const uint64_t *masks, size_t words,
                                const uint64_t graph[][2])
{
    size_t anchors = d->stats->anchors;
    /* The clique alone already refutes this target. Preserve that small proof. */
    if (anchors > q->target) d->active_pairs = false;
    if (anchors > q->target + 1) anchors = q->target + 1;
    if (anchors < 2) return;
    d->capacity = anchors * (anchors - 1) / 2;
    if (d->active_pairs) {
        d->capacity = q->positives * (q->positives - 1) / 2;
        d->graph = nsat_alloc(&d->budget, q->positives, sizeof(*d->graph));
        if (!d->graph || !nsat_copy(&d->budget, d->graph, graph,
            q->positives * sizeof(*d->graph))) return;
    }
    d->pairs = nsat_alloc(&d->budget, d->capacity, sizeof(*d->pairs));
    if (!d->pairs) return;
    size_t later = d->active_pairs ? q->positives : anchors;
    size_t earlier = d->active_pairs ? q->positives - 1 : anchors;
    for (size_t a = 0; a < earlier; ++a) {
        size_t first_row = d->order[a];
        for (size_t b = a + 1; b < later; ++b) {
            size_t second_row = d->order[b], z = 0;
            if (b >= anchors) {
                if (!nsat_work(&d->budget, 1)) return;
                if (!(graph[first_row][second_row / 64] &
                    (UINT64_C(1) << (second_row % 64)))) continue;
            }
            for (; z < q->negatives; ++z) {
                if (!nsat_work(&d->budget, words)) return;
                const uint64_t *first = masks + (first_row * q->negatives + z) * words;
                const uint64_t *second = masks + (second_row * q->negatives + z) * words;
                if (!has_common_difference(first, second, NULL, words)) break;
            }
            if (z == q->negatives) {
                nsat_stop(&d->budget, NSAT_ERROR, NSAT_STOP_PROOF);
                return;
            }
            d->pairs[d->pair_count++] = (IncompatiblePair){{first_row, second_row}, z};
        }
    }
}

static void discover(const NsatCoverInput *q, Domain *d)
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
                    break;
                }
            }
        }
    }
ordered:
    if (d->budget.status == NSAT_OK) choose_order(d, q, graph, degree, triangles);
    if (d->budget.status == NSAT_OK) retain_ordered_pairs(q, d, masks, words, graph);
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

NsatStatus nsat_cover_distance_one(const NsatCoverInput *q,
                                    NsatSolver *s, NsatVerifier *v, NsatLimits limits,
                                    NsatCoverDistanceStats *stats)
{
    if (!stats) return NSAT_ERROR;
    *stats = (NsatCoverDistanceStats){0};
    if (!q || !s || !v || !q->on || !q->off || !q->inputs || !q->positives ||
        !q->negatives || !q->target) return NSAT_ERROR;
    NsatBudget budget = {.limits = limits};
    NsatClauseId *hints = nsat_alloc(&budget, q->inputs, sizeof(*hints));
    unsigned char *seen = nsat_alloc(&budget, q->inputs, sizeof(*seen));
    NsatStatus status = budget.status;
    if (status != NSAT_OK) goto done;
    for (size_t r = 0; r < q->positives; ++r) {
        if (!nsat_copy(&budget, seen, NULL, q->inputs)) goto stopped;
        const int *on = q->on + r * q->inputs;
        for (size_t z = 0; z < q->negatives; ++z) {
            const int *off = q->off + z * q->inputs;
            size_t difference = q->inputs, i = 0;
            for (; i < q->inputs; ++i) {
                if (stats->comparisons == max_distance_comparisons) {
                    stats->capped = true;
                    goto done;
                }
                if (!nsat_work(&budget, 1)) goto stopped;
                ++stats->comparisons;
                if (on[i] == off[i]) continue;
                if (difference != q->inputs) break;
                difference = i;
            }
            if (i != q->inputs || difference == q->inputs || seen[difference]) continue;
            seen[difference] = 1;
            for (size_t t = 0; t < q->target; ++t) {
                if (stats->clauses == MAX_LEMMAS) {
                    stats->capped = true;
                    goto done;
                }
                NsatClauseId first = q->target * (q->inputs + q->negatives) +
                    r * (q->target * q->inputs + 1) + t * q->inputs + 1;
                /* Negating the shortcut assigns the ON row to this slot
                 * and denies its required selector. The other ON binaries
                 * falsify all remaining literals of the OFF clause. */
                size_t count = 0;
                for (size_t c = 0; c < q->inputs; ++c) {
                    if (!nsat_work(&budget, 1)) goto stopped;
                    if (c != difference) hints[count++] = first + c;
                }
                hints[count++] = t * (q->inputs + q->negatives) + q->inputs + z + 1;
                int selector = nsat_cover_selector((int)q->inputs, (int)t,
                    (int)difference, on[difference]);
                int clause[] = {-assigned(q, r, t), selector};
                status = nsat_add_root_rup(s, v, clause, 2, hints, count, NULL);
                if (status != NSAT_OK) goto done;
                ++stats->clauses;
            }
        }
    }
    goto done;
stopped:
    status = budget.status;
done:
    nsat_free(&budget, hints, q->inputs * sizeof(*hints));
    nsat_free(&budget, seen, q->inputs * sizeof(*seen));
    stats->work = budget.stats.work;
    stats->stop = budget.stats.stop;
    return status;
}

static NsatStatus add_incompatible_pairs(const NsatCoverInput *q, Domain *d,
                                         NsatSolver *s, NsatVerifier *v)
{
    if (!d->pair_count) return NSAT_OK;
    const size_t hint_count = q->inputs + 1;
    NsatClauseId *hints = nsat_alloc(&d->budget, hint_count, sizeof(*hints));
    if (!hints) return d->budget.status;
    NsatStatus status = NSAT_OK;
    size_t position[MAX_ROWS];
    for (size_t r = 0; r < q->positives; ++r) position[source_row(q, r)] = r;
    for (size_t p = 0; p < d->pair_count; ++p) {
        IncompatiblePair pair = d->pairs[p];
        const int *off = q->off + pair.off * q->inputs;
        const int *first = q->on + pair.rows[0] * q->inputs;
        const int *second = q->on + pair.rows[1] * q->inputs;
        size_t first_row = position[pair.rows[0]], second_row = position[pair.rows[1]];
        size_t first_slot = first_row;
        size_t slots = first_slot + 1;
        bool non_anchor_pair = d->active_pairs && first_row >= d->stats->anchors;
        if (non_anchor_pair) {
            first_slot = 0;
            if (slots > q->target) slots = q->target;
        }
        for (size_t t = first_slot; t < slots; ++t) {
            if (non_anchor_pair) {
                if (!nsat_work(&d->budget, 1)) {
                    status = d->budget.status;
                    goto done;
                }
                /* Earlier pairs pin each anchor and exclude both of its
                 * incompatible rows. Missing graph edges remain eligible. */
                if (t < d->stats->anchors) {
                    const uint64_t *adjacent = d->graph[source_row(q, t)];
                    if ((adjacent[pair.rows[0] / 64] & (UINT64_C(1) << (pair.rows[0] % 64))) ||
                        (adjacent[pair.rows[1] / 64] & (UINT64_C(1) << (pair.rows[1] % 64)))) continue;
                }
                if (d->stats->pairs == MAX_LEMMAS) {
                    d->stats->discovery_capped = true;
                    goto done;
                }
            }
            int clause[] = {-assigned(q, first_row, t), -assigned(q, second_row, t)};
            for (size_t i = 0; i < q->inputs; ++i) {
                if (!nsat_work(&d->budget, 1)) {
                    status = d->budget.status;
                    goto done;
                }
                size_t row;
                if (first[i] == off[i]) {
                    row = first_row;
                } else if (second[i] == off[i]) {
                    row = second_row;
                } else {
                    nsat_stop(&d->budget, NSAT_ERROR, NSAT_STOP_PROOF);
                    status = d->budget.status;
                    goto done;
                }
                hints[i] = row_term_first(q, row, t) + i;
            }
            hints[q->inputs] = excludes_off(q, t, pair.off);
            status = nsat_add_root_rup(s, v, clause, 2, hints, hint_count, NULL);
            if (status != NSAT_OK) goto done;
            ++d->stats->pairs;
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

NsatStatus nsat_cover_solve_default(const NsatCoverInput *input, NsatLimits limits,
                                    NsatSolver **s, NsatVerifier **v,
                                    NsatStopReason *stop, NsatCoverDomainStats *stats)
{
    *s = NULL;
    *v = NULL;
    *stop = NSAT_STOP_NONE;
    *stats = (NsatCoverDomainStats){0};
    NsatCoverInput q = *input;
    Domain d = {.stats = stats, .budget = {.limits = limits}, .active_pairs = true};
    const size_t max = (size_t)INT_MAX / 2;
    bool eligible = q.on && q.off && q.inputs &&
        q.positives && q.positives <= MAX_ROWS && q.negatives && q.negatives <= MAX_OFF &&
        q.target && q.target <= q.positives && !q.row_order && !q.order_count &&
        q.inputs <= max - q.positives && q.target <= max / (q.inputs + q.positives) &&
        q.inputs <= SIZE_MAX / sizeof(int) / q.positives &&
        q.inputs <= SIZE_MAX / sizeof(int) / q.negatives;
    if (eligible) {
        discover(&q, &d);
        q.row_order = d.order;
        q.order_count = q.positives;
    }
    NsatStatus status = d.budget.status;
    if (status == NSAT_OK) status = nsat_verify_init_cover(v, q.on, q.off, q.inputs,
        q.positives, q.negatives, q.target, q.row_order, q.order_count, limits);
    if (status == NSAT_OK) status = nsat_cover_load_solver(s, &q, *v, limits, stop);
    if (status == NSAT_OK && eligible) status = add_staircase(&q, &d, *s, *v);
    if (status == NSAT_OK && eligible) status = add_incompatible_pairs(&q, &d, *s, *v);
    nsat_free(&d.budget, d.pairs, d.capacity * sizeof(*d.pairs));
    nsat_free(&d.budget, d.graph, q.positives * sizeof(*d.graph));
    stats->prep_work = d.budget.stats.work;
    if (status == NSAT_OK) {
        NsatSearchOptions options = {.sticky_priority = true, .cache_root_units = true};
        if (eligible) {
            options.priority_first = assigned(&q, 0, 0);
            options.priority_count = (int)(q.positives * q.target);
            options.row_width = (int)q.target;
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
