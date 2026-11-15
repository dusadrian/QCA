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
} HullFrame;

/* Visit legal ON subsets in increasing row order. A suffix sum bounds the
 * weight of every descendant, even if its rows cannot share a legal hull.
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
    if (!frames || !hulls || !selected || !nsat_copy(&check->budget, selected, NULL, row_bytes) ||
        !nsat_work(&check->budget, stride)) goto done;
    memset(hulls, 0xff, stride * sizeof(*hulls));
    frames[0] = (HullFrame){0};
    for (size_t depth = 0;;) {
        if (!nsat_work(&check->budget, 1)) goto done;
        HullFrame *frame = frames + depth;
        uint64_t *hull = hulls + depth * stride, *child = hull + stride;
        if (frame->next < rows && frame->weight + suffix[frame->next] >= cutoff) {
            size_t r = frame->next++;
            if (!extend_hull(check, hull, r, child)) continue;
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

NsatStatus ipc_verify_weighted_bound(const NsatCoverInput *q, const uint64_t *weights,
                                     size_t count, uint64_t maximum, int lower,
                                     NsatLimits limits, NsatStats *stats)
{
    WeightedCheck check = {.budget = {.limits = limits}, .input = q, .weights = weights};
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
    if (check_hulls(&check, suffix, maximum, by_weight, total - others))
        status = by_weight ? NSAT_UNSAT : check_cover(&check, (size_t)target);
    if (status == NSAT_SAT) nsat_stop(&check.budget, NSAT_ERROR, NSAT_STOP_PROOF);
done:
    nsat_free(&check.budget, suffix, (count + 1) * sizeof(*suffix));
    if (on) nsat_free(&check.budget, on, count * check.input_words * sizeof(*on));
    if (off) nsat_free(&check.budget, off, q->negatives * check.input_words * sizeof(*off));
    nsat_free(&check.budget, check.columns, check.capacity * check.row_words * sizeof(uint64_t));
    if (stats) *stats = check.budget.stats;
    return check.budget.status == NSAT_OK ? status : check.budget.status;
}
