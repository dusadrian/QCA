/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "ipc_weighted_verify.h"
#include "native_sat_support.h"

typedef struct {
    NsatBudget budget;
    const NsatCoverInput *input;
    const uint64_t *weights;
    const uint64_t *on, *off;
    size_t input_words, row_words;
    uint64_t *columns;
    size_t used, capacity;
    const IpcWeightHints *hints;
} WeightedCheck;

static void *matrix(NsatBudget *budget, size_t rows, size_t cols, size_t width)
{
    size_t cells;
    if (nsat_mul(rows, cols, &cells)) return nsat_alloc(budget, cells, width);
    nsat_stop(budget, NSAT_ERROR, NSAT_STOP_INPUT);
    return NULL;
}

static bool pack_points(WeightedCheck *check, const int *points, size_t count, uint64_t *bits)
{
    size_t width = check->input->inputs, words = check->input_words;
    for (size_t r = 0; r < count; ++r) {
        if (!nsat_work(&check->budget, width)) return false;
        memset(bits + r * words, 0, words * sizeof(*bits));
        for (size_t i = 0; i < width; ++i) {
            int value = points[r * width + i];
            if (value == 2) bits[r * words + i / 64] |= UINT64_C(1) << (i % 64);
            else if (value != 1) return nsat_stop(&check->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        }
    }
    return true;
}

/* Geometry uses two input masks, for the literals fixed to one and zero.
 * A hull is legal precisely when it excludes each original OFF point. */
static bool extend_hull(WeightedCheck *check, const uint64_t *parent, size_t row,
                        uint64_t *child)
{
    size_t words = check->input_words, negatives = check->input->negatives;
    if (!nsat_work(&check->budget, words)) return false;
    const uint64_t *point = check->on + row * words;
    for (size_t w = 0; w < words; ++w) {
        child[w] = parent[w] & point[w];
        child[words + w] = parent[words + w] & ~point[w];
    }
    if (words == 1) {
        size_t n = 0;
        while (n < negatives && ((child[0] & ~check->off[n]) || (child[1] & check->off[n])))
            ++n;
        return nsat_work(&check->budget, n + (n < negatives)) && n == negatives;
    }
    for (size_t n = 0; n < negatives; ++n) {
        const uint64_t *off = check->off + n * words;
        if (!nsat_work(&check->budget, words)) return false;
        size_t w = 0;
        while (w < words && !(child[w] & ~off[w]) && !(child[words + w] & off[w])) ++w;
        if (w == words) return false;
    }
    return true;
}

static bool retain(WeightedCheck *check, const uint64_t *selected)
{
    size_t stride = check->row_words * sizeof(*selected);
    if (check->used == check->capacity) {
        if (check->capacity > SIZE_MAX / 2)
            return nsat_stop(&check->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        size_t capacity = check->capacity ? 2 * check->capacity : 64;
        uint64_t *columns = nsat_alloc(&check->budget, capacity, stride);
        if (!columns) return false;
        if (!nsat_copy(&check->budget, columns, check->columns, check->used * stride)) {
            nsat_free(&check->budget, columns, capacity * stride);
            return false;
        }
        nsat_free(&check->budget, check->columns, check->capacity * stride);
        check->columns = columns;
        check->capacity = capacity;
    }
    if (!nsat_work(&check->budget, check->row_words)) return false;
    memcpy(check->columns + check->used++ * check->row_words, selected, stride);
    return true;
}

typedef struct {
    size_t last, next;
    uint64_t weight;
    bool extended;
    bool primed;
    uint64_t rest;
} HullFrame;

/* Visit legal ON subsets in increasing row order. The sum of weights of
 * compatible later rows bounds all descendants of the current hull.
 * Only globally maximal heavy subsets can be needed in a cover. */
static bool check_hulls(WeightedCheck *check, const uint64_t *suffix, uint64_t maximum,
                        bool by_weight, uint64_t minimum)
{
    size_t rows = check->input->positives, stride = 2 * check->input_words;
    size_t row_bytes = check->row_words * sizeof(uint64_t);
    uint64_t cutoff = by_weight ? maximum + 1 : minimum;
    HullFrame *frames = nsat_alloc(&check->budget, rows + 1, sizeof(*frames));
    uint64_t *hulls = matrix(&check->budget, rows + 2, stride, sizeof(*hulls));
    uint64_t *selected = nsat_alloc(&check->budget, check->row_words, sizeof(*selected));
    bool ok = false;
    uint64_t *compatible = NULL;
    if (!frames || !hulls || !selected) goto done;
    if (rows + 1 <= nsat_memory_available(&check->budget) / row_bytes) {
        compatible = matrix(&check->budget, rows + 1, check->row_words, sizeof(*compatible));
        if (!compatible) goto done;
    }
    if (!nsat_copy(&check->budget, selected, NULL, row_bytes) ||
        !nsat_work(&check->budget, stride)) goto done;
    memset(hulls, 0xff, stride * sizeof(*hulls));
    frames[0] = (HullFrame){0};
    for (size_t depth = 0;;) {
        if (!nsat_work(&check->budget, 1)) goto done;
        HullFrame *frame = frames + depth;
        uint64_t *hull = hulls + depth * stride, *child = hull + stride;
        if (frame->next < rows && frame->weight + suffix[frame->next] >= cutoff) {
            size_t r = frame->next;
            if (compatible) {
                uint64_t *bits = compatible + depth * check->row_words;
                const uint64_t *parent = depth ? bits - check->row_words : NULL;
                if (!frame->primed) {
                    if (!nsat_work(&check->budget, check->row_words)) goto done;
                    memset(bits, 0, row_bytes);
                    for (size_t next = r; next < rows; ++next) {
                        if (!nsat_work(&check->budget, 1)) goto done;
                        uint64_t bit = UINT64_C(1) << (next % 64);
                        if (parent && !(parent[next / 64] & bit)) continue;
                        if (extend_hull(check, hull, next, child)) {
                            bits[next / 64] |= bit;
                            frame->rest += check->weights[next];
                        }
                        if (check->budget.status != NSAT_OK) goto done;
                    }
                    frame->primed = true;
                }
                if (frame->weight + frame->rest >= cutoff) {
                    for (; r < rows; ++r) {
                        if (!nsat_work(&check->budget, 1)) goto done;
                        if (bits[r / 64] & (UINT64_C(1) << (r % 64))) break;
                    }
                } else {
                    r = rows;
                }
                if (r < rows) {
                    /* Keep the original compatibility mask for siblings. */
                    if (!nsat_work(&check->budget, check->input_words)) goto done;
                    for (size_t w = 0; w < check->input_words; ++w) {
                        uint64_t point = check->on[r * check->input_words + w];
                        child[w] = hull[w] & point;
                        child[check->input_words + w] = hull[check->input_words + w] & ~point;
                    }
                    frame->next = r + 1;
                    frame->rest -= check->weights[r];
                }
            } else {
                ++frame->next;
                if (!extend_hull(check, hull, r, child)) continue;
            }
            if (r < rows) {
                uint64_t weight = frame->weight + check->weights[r];
                if (weight > maximum) {
                    nsat_stop(&check->budget, NSAT_ERROR, NSAT_STOP_PROOF);
                    goto done;
                }
                frame->extended = true;
                selected[r / 64] |= UINT64_C(1) << (r % 64);
                frames[++depth] = (HullFrame){.last = r, .next = r + 1, .weight = weight};
                continue;
            }
        }
        if (!depth) break;
        if (!by_weight && !frame->extended && frame->weight >= minimum) {
            size_t r = 0;
            for (; r < frame->last; ++r) {
                if (!nsat_work(&check->budget, 1)) goto done;
                if ((selected[r / 64] >> (r % 64)) & 1u) continue;
                if (extend_hull(check, hull, r, child)) break;
            }
            if (r == frame->last && !retain(check, selected)) goto done;
        }
        size_t last = frame->last;
        selected[last / 64] &= ~(UINT64_C(1) << (last % 64));
        --depth;
    }
    ok = nsat_poll(&check->budget);
done:
    nsat_free(&check->budget, compatible, (rows + 1) * row_bytes);
    nsat_free(&check->budget, selected, row_bytes);
    nsat_free(&check->budget, hulls, (rows + 2) * stride * sizeof(*hulls));
    nsat_free(&check->budget, frames, (rows + 1) * sizeof(*frames));
    return ok;
}

static uint64_t row_weight(const WeightedCheck *check, const uint64_t *rows,
                           const uint64_t *cube)
{
    uint64_t sum = 0;
    for (size_t w = 0; w < check->row_words; ++w)
        for (uint64_t bits = rows[w] & cube[w]; bits; bits &= bits - 1)
            sum += check->weights[w * 64 + (size_t)__builtin_ctzll(bits)];
    return sum;
}

typedef struct {
    size_t next, anchor;
    uint64_t minimum;
} CoverFrame;

/* Independently branch on a row with the fewest eligible extensions. Every
 * completion contains a cube covering it. If U bounds each remaining trace,
 * every cube of a k-cube completion carries at least W - (k - 1) U. */
static NsatStatus check_cover(WeightedCheck *check, size_t target)
{
    size_t words = check->row_words, count = check->input->positives;
    if (!check->used) return NSAT_UNSAT;
    CoverFrame *frames = nsat_alloc(&check->budget, target, sizeof(*frames));
    uint64_t *remaining = matrix(&check->budget, target, words, sizeof(*remaining));
    uint64_t *prices = nsat_alloc(&check->budget, check->used, sizeof(*prices));
    size_t *frequency = nsat_alloc(&check->budget, count, sizeof(*frequency));
    NsatStatus status = NSAT_UNKNOWN;
    if (!frames || !remaining || !prices || !frequency ||
        !nsat_work(&check->budget, words)) goto done;
    for (size_t w = 0; w < words; ++w) remaining[w] = UINT64_MAX;
    if (count % 64) remaining[words - 1] = (UINT64_C(1) << (count % 64)) - 1;
    bool entering = true;
    for (size_t depth = 0;;) {
        uint64_t *rows = remaining + depth * words;
        CoverFrame *frame = frames + depth;
        if (entering) {
            if (!nsat_work(&check->budget, count)) goto done;
            uint64_t need = row_weight(check, rows, rows), largest = 0;
            size_t slots = target - depth;
            for (size_t c = 0; c < check->used; ++c) {
                if (!nsat_work(&check->budget, words)) goto done;
                const uint64_t *cube = check->columns + c * words;
                size_t w = 0;
                while (w < words && !(rows[w] & ~cube[w])) ++w;
                if (w == words) {
                    status = NSAT_SAT;
                    goto done;
                }
                prices[c] = row_weight(check, rows, cube);
                if (prices[c] > largest) largest = prices[c];
            }
            *frame = (CoverFrame){.next = check->used};
            if (slots > 1 && (!need || (largest && (need - 1) / largest < slots))) {
                uint64_t others = largest <= need / (slots - 1) ? largest * (slots - 1) : need;
                frame->minimum = need - others;
                memset(frequency, 0, count * sizeof(*frequency));
                for (size_t c = 0; c < check->used; ++c) {
                    if (!nsat_work(&check->budget, words)) goto done;
                    if (prices[c] < frame->minimum) continue;
                    const uint64_t *cube = check->columns + c * words;
                    for (size_t w = 0; w < words; ++w)
                        for (uint64_t bits = rows[w] & cube[w]; bits; bits &= bits - 1)
                            ++frequency[w * 64 + (size_t)__builtin_ctzll(bits)];
                }
                size_t anchor = SIZE_MAX;
                for (size_t w = 0; w < words; ++w) {
                    for (uint64_t bits = rows[w]; bits; bits &= bits - 1) {
                        size_t r = w * 64 + (size_t)__builtin_ctzll(bits);
                        if (anchor == SIZE_MAX || frequency[r] < frequency[anchor]) anchor = r;
                    }
                }
                frame->anchor = anchor;
                if (anchor != SIZE_MAX && frequency[anchor]) frame->next = 0;
            }
        }
        entering = false;
        while (frame->next < check->used) {
            if (!nsat_work(&check->budget, words)) goto done;
            const uint64_t *cube = check->columns + frame->next++ * words;
            size_t r = frame->anchor;
            if (!((cube[r / 64] >> (r % 64)) & 1u) ||
                row_weight(check, rows, cube) < frame->minimum) continue;
            for (size_t w = 0; w < words; ++w) rows[words + w] = rows[w] & ~cube[w];
            ++depth;
            entering = true;
            break;
        }
        if (entering) continue;
        if (!depth) {
            if (nsat_poll(&check->budget)) status = NSAT_UNSAT;
            break;
        }
        --depth;
    }
done:
    nsat_free(&check->budget, frequency, count * sizeof(*frequency));
    nsat_free(&check->budget, prices, check->used * sizeof(*prices));
    nsat_free(&check->budget, remaining, target * words * sizeof(*remaining));
    nsat_free(&check->budget, frames, target * sizeof(*frames));
    return status;
}

static uint64_t mask_weight(uint64_t bits, const uint64_t *weights)
{
    uint64_t sum = 0;
    for (; bits; bits &= bits - 1) sum += weights[__builtin_ctzll(bits)];
    return sum;
}

static const uint64_t *find_hint(WeightedCheck *check, uint64_t residual, size_t slots)
{
    const IpcWeightHints *hints = check->hints;
    if (!hints || !hints->entries || !hints->weights ||
        hints->rows != check->input->positives || hints->count > IPC_WEIGHT_HINT_BUCKETS) return NULL;
    size_t bucket = ipc_weight_hint_bucket(residual, (int)slots);
    for (size_t probe = 0; probe < IPC_WEIGHT_HINT_BUCKETS; ++probe) {
        if (!nsat_work(&check->budget, 1)) return NULL;
        IpcWeightHint h = hints->entries[bucket];
        if (!h.index) return NULL;
        if (h.residual == residual && h.slots == (int)slots) {
            if (h.index > hints->count) return NULL;
            const uint64_t *weights = hints->weights + (h.index - 1) * hints->rows;
            uint64_t sum = 0;
            for (uint64_t bits = residual; bits; bits &= bits - 1) {
                if (!nsat_work(&check->budget, 1)) return NULL;
                uint64_t w = weights[__builtin_ctzll(bits)];
                if (w > UINT64_MAX - sum) return NULL;
                sum += w;
            }
            return weights;
        }
        bucket = (bucket + 1) & (IPC_WEIGHT_HINT_BUCKETS - 1);
    }
    return NULL;
}

/* Candidates dropped at a parent cannot appear in a completion at any
 * child. Derive every filter using this checker's own chart and arithmetic. */
static size_t filter_chart(WeightedCheck *check, uint64_t residual, size_t slots,
                           uint64_t *choices, size_t count, const uint64_t *weights,
                           uint64_t *need_out, uint64_t *max_out,
                           uint64_t previous_need, uint64_t previous_max, bool *used)
{
    uint64_t need = mask_weight(residual, weights), maximum = 0;
    for (size_t c = 0; c < count; ++c) {
        if (!nsat_work(&check->budget, 1)) return 0;
        uint64_t value = mask_weight(choices[c] & residual, weights);
        if (value > maximum) maximum = value;
    }
    if (need && (!maximum || (need - 1) / maximum >= slots)) return 0;
    if (need_out) *need_out = need;
    if (max_out) *max_out = maximum;
    if (used && (!maximum || (previous_max &&
        (long double)need / maximum < (long double)previous_need / previous_max)))
        return count;
    if (used) *used = true;
    uint64_t other = maximum <= need / (slots - 1) ? maximum * (slots - 1) : need;
    uint64_t minimum = need - other;
    size_t keep = 0;
    for (size_t c = 0; c < count; ++c) {
        if (!nsat_work(&check->budget, 1)) return 0;
        if (mask_weight(choices[c] & residual, weights) >= minimum)
            choices[keep++] = choices[c];
    }
    return keep;
}

static NsatStatus visit_chart(WeightedCheck *check, uint64_t residual, size_t slots,
                              const uint64_t *parent, size_t count, uint64_t *choices,
                              const uint64_t *weights)
{
    if (!nsat_work(&check->budget, 1)) return check->budget.status;
    size_t keep = 0;
    for (size_t c = 0; c < count; ++c) {
        if (!nsat_work(&check->budget, 1)) return check->budget.status;
        uint64_t mask = parent[c];
        if (!(residual & ~mask)) return NSAT_SAT;
        if ((residual & mask) && slots > 1) choices[keep++] = mask;
    }
    if (!keep) return NSAT_UNSAT;
    uint64_t need = 0, maximum = 0;
    keep = filter_chart(check, residual, slots, choices, keep, weights,
                        &need, &maximum, 0, 0, NULL);
    if (check->budget.status != NSAT_OK) return check->budget.status;
    if (!keep) return NSAT_UNSAT;

    const uint64_t *local = find_hint(check, residual, slots);
    if (check->budget.status != NSAT_OK) return check->budget.status;
    if (local) {
        bool used = false;
        keep = filter_chart(check, residual, slots, choices, keep, local,
                            NULL, NULL, need, maximum, &used);
        if (check->budget.status != NSAT_OK) return check->budget.status;
        if (!keep) return NSAT_UNSAT;
        /* The checked vector can also bound residuals of this subproblem. */
        if (used) weights = local;
    }
    size_t frequency[64] = {0};
    if (!nsat_work(&check->budget, check->input->positives)) return check->budget.status;
    for (size_t c = 0; c < keep; ++c) {
        if (!nsat_work(&check->budget, 1)) return check->budget.status;
        for (uint64_t bits = choices[c] & residual; bits; bits &= bits - 1)
            ++frequency[__builtin_ctzll(bits)];
    }
    unsigned anchor = __builtin_ctzll(residual);
    for (uint64_t bits = residual & (residual - 1); bits; bits &= bits - 1) {
        unsigned row = __builtin_ctzll(bits);
        if (frequency[row] < frequency[anchor]) anchor = row;
    }
    uint64_t bit = UINT64_C(1) << anchor;
    for (size_t c = 0; c < keep;) {
        if (!nsat_work(&check->budget, 1)) return check->budget.status;
        if (!(choices[c] & bit)) {
            ++c;
            continue;
        }
        uint64_t trace = choices[c] & residual;
        bool dominated = false;
        for (size_t d = 0; d < keep; ++d) {
            if (!nsat_work(&check->budget, 1)) return check->budget.status;
            if (d == c || (trace & ~choices[d])) continue;
            if (d < c || trace != (choices[d] & residual)) {
                dominated = true;
                break;
            }
        }
        if (dominated) {
            ++c;
            continue;
        }
        /* Partition covers by the first chosen column covering the anchor.
         * It covers no residual row in this child and is excluded in later siblings. */
        size_t following = keep - c - 1;
        if (!nsat_work(&check->budget, following)) return check->budget.status;
        memmove(choices + c, choices + c + 1, following * sizeof(*choices));
        --keep;
        NsatStatus status = visit_chart(check, residual & ~trace, slots - 1,
                                        choices, keep, choices + check->used, weights);
        if (status != NSAT_UNSAT) return status;
    }
    return NSAT_UNSAT;
}

static NsatStatus check_cover_inherited(WeightedCheck *check, size_t target)
{
    if (!check->used) return NSAT_UNSAT;
    /* A compact scalar chart is optional under a tight memory limit. */
    if (check->used > nsat_memory_available(&check->budget) / (target + 1) / sizeof(uint64_t))
        return check_cover(check, target);
    NsatStatus status = NSAT_UNKNOWN;
    uint64_t *chart = matrix(&check->budget, target + 1, check->used, sizeof(*chart));
    if (!chart) return check->budget.status;
    /* Stable wide-first order of this independently enumerated chart. */
    size_t offsets[65] = {0}, rows = check->input->positives;
    for (size_t c = 0; c < check->used; ++c) {
        if (!nsat_work(&check->budget, 1)) goto done;
        ++offsets[__builtin_popcountll(check->columns[c])];
    }
    for (size_t size = rows, next = 0; size; --size) {
        size_t count = offsets[size];
        offsets[size] = next;
        next += count;
    }
    for (size_t c = 0; c < check->used; ++c) {
        if (!nsat_work(&check->budget, 1)) goto done;
        uint64_t mask = check->columns[c];
        chart[offsets[__builtin_popcountll(mask)]++] = mask;
    }
    uint64_t residual = rows == 64 ? UINT64_MAX : (UINT64_C(1) << rows) - 1;
    status = visit_chart(check, residual, target, chart, check->used,
                         chart + check->used, check->weights);
    if (!nsat_poll(&check->budget)) status = check->budget.status;
done:
    nsat_free(&check->budget, chart, (target + 1) * check->used * sizeof(*chart));
    return status;
}

NsatStatus ipc_verify_weighted_hint_bound(const NsatCoverInput *q, const uint64_t *weights,
                                     size_t count, uint64_t maximum, int lower,
                                     NsatLimits limits, NsatStats *stats,
                                     const IpcWeightHints *hints)
{
    WeightedCheck check = {.budget = {.limits = limits}, .input = q, .weights = weights,
                           .hints = hints};
    NsatStatus status = NSAT_UNKNOWN;
    uint64_t *on = NULL, *off = NULL, *suffix = NULL;
    if (!q || !weights || !q->on || !q->off || !q->inputs || !count ||
        count != q->positives || count > INT_MAX || !q->negatives ||
        lower < 2 || (size_t)lower > q->target || q->target > count ||
        q->inputs > SIZE_MAX / sizeof(int) / count ||
        q->inputs > SIZE_MAX / sizeof(int) / q->negatives) {
        nsat_stop(&check.budget, NSAT_ERROR, NSAT_STOP_INPUT);
        goto done;
    }
    check.input_words = (q->inputs - 1) / 64 + 1;
    check.row_words = (count - 1) / 64 + 1;
    suffix = nsat_alloc(&check.budget, count + 1, sizeof(*suffix));
    if (!suffix || !nsat_work(&check.budget, count)) goto done;
    suffix[count] = 0;
    for (size_t r = count; r-- > 0;) {
        if (weights[r] > UINT64_MAX - suffix[r + 1]) {
            nsat_stop(&check.budget, NSAT_ERROR, NSAT_STOP_INPUT);
            goto done;
        }
        suffix[r] = weights[r] + suffix[r + 1];
    }
    uint64_t total = suffix[0], target = (uint64_t)lower - 1;
    if (maximum > total) {
        nsat_stop(&check.budget, NSAT_ERROR, NSAT_STOP_PROOF);
        goto done;
    }
    bool by_weight = total && maximum <= (total - 1) / target;
    uint64_t others = target == 1 ? 0 :
        maximum <= total / (target - 1) ? maximum * (target - 1) : total;
    on = matrix(&check.budget, count, check.input_words, sizeof(*on));
    off = matrix(&check.budget, q->negatives, check.input_words, sizeof(*off));
    if (!on || !off || !pack_points(&check, q->on, count, on) ||
        !pack_points(&check, q->off, q->negatives, off)) goto done;
    check.on = on;
    check.off = off;
    bool hulls_ok = check_hulls(&check, suffix, maximum, by_weight, total - others);
    if (hulls_ok) {
        if (by_weight) status = NSAT_UNSAT;
        else if (check.row_words == 1) status = check_cover_inherited(&check, (size_t)target);
        else status = check_cover(&check, (size_t)target);
    }
    if (status == NSAT_SAT) nsat_stop(&check.budget, NSAT_ERROR, NSAT_STOP_PROOF);
done:
    nsat_free(&check.budget, suffix, (count + 1) * sizeof(*suffix));
    if (on) nsat_free(&check.budget, on, count * check.input_words * sizeof(*on));
    if (off) nsat_free(&check.budget, off, q->negatives * check.input_words * sizeof(*off));
    nsat_free(&check.budget, check.columns, check.capacity * check.row_words * sizeof(uint64_t));
    if (stats) *stats = check.budget.stats;
    return check.budget.status == NSAT_OK ? status : check.budget.status;
}

NsatStatus ipc_verify_weighted_bound(const NsatCoverInput *q,
    const uint64_t *weights, size_t count, uint64_t maximum, int lower,
    NsatLimits limits, NsatStats *stats)
{
    return ipc_verify_weighted_hint_bound(q, weights, count, maximum, lower, limits, stats, NULL);
}
