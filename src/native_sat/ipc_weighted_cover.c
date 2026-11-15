/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "ipc_weighted_cover.h"
#include "ipc_weighted_verify.h"
#include "native_sat_support.h"

#include <math.h>

typedef struct {
    NsatBudget budget;
    IpcWeightedCoverResult result;
    const NsatCoverInput *input;
    uint64_t *diff, *columns;
    size_t capacity, input_words, row_words, mask_words;
    uint64_t (*weight_bytes)[256];
} WeightedCover;

static uint64_t row_bit(unsigned row)
{
    return UINT64_C(1) << row;
}

static void *allocate_matrix(NsatBudget *budget, size_t rows, size_t columns, size_t width)
{
    size_t cells;
    if (!nsat_mul(rows, columns, &cells)) {
        nsat_stop(budget, NSAT_ERROR, NSAT_STOP_INPUT);
        return NULL;
    }
    return nsat_alloc(budget, cells, width);
}

static size_t next_row(const uint64_t *mask, size_t words, size_t first)
{
    size_t word = first / 64;
    if (word >= words) return SIZE_MAX;
    uint64_t bits = mask[word] & (UINT64_MAX << (first % 64));
    while (!bits && ++word < words) bits = mask[word];
    return bits ? 64 * word + (size_t)__builtin_ctzll(bits) : SIZE_MAX;
}

static void select_row(uint64_t *mask, size_t row)
{
    mask[row / 64] |= row_bit((unsigned)(row % 64));
}

static void clear_row(uint64_t *mask, size_t row)
{
    mask[row / 64] &= ~row_bit((unsigned)(row % 64));
}

static bool contains_row(const uint64_t *mask, size_t row)
{
    return (mask[row / 64] & row_bit((unsigned)(row % 64))) != 0;
}

static void all_rows(const WeightedCover *w, uint64_t *mask)
{
    for (size_t i = 0; i < w->mask_words; ++i) mask[i] = UINT64_MAX;
    size_t tail = w->input->positives % 64;
    if (tail) mask[w->mask_words - 1] = (UINT64_C(1) << tail) - 1;
}

static const uint64_t *column(const WeightedCover *w, size_t index)
{
    return w->columns + index * w->mask_words;
}

static bool prepare_differences(WeightedCover *w)
{
    const NsatCoverInput *q = w->input;
    for (size_t r = 0; r < q->positives; ++r) {
        for (size_t z = 0; z < q->negatives; ++z) {
            uint64_t *diff = w->diff + r * w->row_words + z * w->input_words;
            const int *on = q->on + r * q->inputs, *off = q->off + z * q->inputs;
            bool distinct = false;
            for (size_t word = 0, i = 0; word < w->input_words; ++word) {
                size_t width = q->inputs - i < 64 ? q->inputs - i : 64;
                if (!nsat_work(&w->budget, width)) return false;
                uint64_t bits = 0;
                for (unsigned bit = 0; bit < width; ++bit, ++i) {
                    int a = on[i], b = off[i];
                    if ((a != 1 && a != 2) || (b != 1 && b != 2))
                        return nsat_stop(&w->budget, NSAT_ERROR, NSAT_STOP_INPUT);
                    if (a != b) bits |= row_bit(bit);
                }
                diff[word] = bits;
                distinct |= bits != 0;
            }
            if (!distinct) return nsat_stop(&w->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        }
    }
    return true;
}

static bool keep_column(WeightedCover *w, const uint64_t *rows)
{
    size_t used = w->result.maximal;
    size_t stride = w->mask_words * sizeof(*rows);
    if (used == w->capacity) {
        if (w->capacity > SIZE_MAX / 2)
            return nsat_stop(&w->budget, NSAT_ERROR, NSAT_STOP_INPUT);
        size_t capacity = w->capacity ? 2 * w->capacity : 64;
        uint64_t *columns = nsat_alloc(&w->budget, capacity, stride);
        if (!columns) return false;
        if (!nsat_copy(&w->budget, columns, w->columns, used * stride)) {
            nsat_free(&w->budget, columns, capacity * stride);
            return false;
        }
        nsat_free(&w->budget, w->columns, w->capacity * stride);
        w->columns = columns;
        w->capacity = capacity;
    }
    if (!nsat_work(&w->budget, w->mask_words)) return false;
    memcpy(w->columns + used * w->mask_words, rows, stride);
    ++w->result.maximal;
    return true;
}

static bool intersects_each_off(WeightedCover *w, const uint64_t *common,
                                const uint64_t *diff)
{
    size_t z = 0;
    if (w->input_words == 1) {
        for (; z < w->row_words && (common[z] & diff[z]); ++z) {}
        return nsat_work(&w->budget, z) && z == w->row_words;
    }
    for (; z < w->row_words; z += w->input_words) {
        if (!nsat_work(&w->budget, w->input_words)) return false;
        size_t word = 0;
        while (word < w->input_words && !(common[z + word] & diff[z + word])) ++word;
        if (word == w->input_words) return false;
    }
    return true;
}

typedef struct {
    size_t row, cursor;
} EnumerationFrame;

/* A legal hull must differ from each OFF row on a common input.
 * Keep earlier compatible rows in the candidate mask so that maximality
 * is checked globally, even though branches add rows in ascending order. */
static bool enumerate(WeightedCover *w)
{
    size_t levels = w->input->positives + 1;
    size_t stride = w->mask_words * sizeof(uint64_t);
    const size_t common_stride = w->row_words * sizeof(uint64_t);
    EnumerationFrame *frames = nsat_alloc(&w->budget, levels, sizeof(*frames));
    uint64_t *common = nsat_alloc(&w->budget, levels, common_stride);
    uint64_t *candidates = nsat_alloc(&w->budget, levels, stride);
    uint64_t *selected = nsat_alloc(&w->budget, w->mask_words, sizeof(*selected));
    bool complete = false;
    if (!frames || !common || !candidates || !selected ||
        !nsat_copy(&w->budget, selected, NULL, stride) ||
        !nsat_work(&w->budget, w->row_words + w->mask_words))
        goto done;
    memset(common, 0xff, common_stride);
    all_rows(w, candidates);
    frames[0] = (EnumerationFrame){0};
    for (size_t depth = 0;;) {
        EnumerationFrame *frame = frames + depth;
        const uint64_t *available = candidates + depth * w->mask_words;
        if (!nsat_work(&w->budget, w->mask_words)) goto done;
        size_t row = next_row(available, w->mask_words, frame->cursor);
        if (row == SIZE_MAX) {
            if (!depth) break;
            clear_row(selected, frame->row);
            --depth;
            continue;
        }
        frame->cursor = row + 1;

        const uint64_t *prior = common + depth * w->row_words;
        uint64_t *next = common + (depth + 1) * w->row_words;
        const uint64_t *diff = w->diff + row * w->row_words;
        if (!nsat_work(&w->budget, w->row_words + w->mask_words)) goto done;
        for (size_t word = 0; word < w->row_words; ++word) next[word] = prior[word] & diff[word];
        uint64_t *viable = candidates + (depth + 1) * w->mask_words;
        memset(viable, 0, stride);
        bool maximal = true;
        for (size_t other = next_row(available, w->mask_words, 0); other != SIZE_MAX;
             other = next_row(available, w->mask_words, other + 1)) {
            if (other == row) continue;
            if (!nsat_work(&w->budget, 1)) goto done;
            if (intersects_each_off(w, next, w->diff + other * w->row_words)) {
                select_row(viable, other);
                maximal = false;
            }
            if (w->budget.status != NSAT_OK) goto done;
        }
        ++w->result.compatible;
        select_row(selected, row);
        if (maximal) {
            if (!keep_column(w, selected)) goto done;
            clear_row(selected, row);
        } else {
            frames[++depth] = (EnumerationFrame){.row = row, .cursor = row + 1};
        }
    }
    complete = true;
done:
    nsat_free(&w->budget, selected, stride);
    nsat_free(&w->budget, candidates, levels * stride);
    nsat_free(&w->budget, common, levels * common_stride);
    nsat_free(&w->budget, frames, levels * sizeof(*frames));
    return complete;
}

static const double pivot_tolerance = 1e-9;

typedef struct {
    size_t rows, stride;
    double *inverse, *rhs;
    bool *surplus;
    size_t *coordinates;
} PricingBasis;

/* Resolve equal ratios in lexicographic order to handle degenerate bases. */
static bool leaves_first(const PricingBasis *basis, const double *direction, size_t a, size_t b)
{
    double difference = basis->rhs[a] / direction[a] - basis->rhs[b] / direction[b];
    if (fabs(difference) > pivot_tolerance) return difference < 0;
    for (size_t r = 0; r < basis->rows; ++r) {
        difference = basis->inverse[a * basis->stride + r] / direction[a] -
            basis->inverse[b * basis->stride + r] / direction[b];
        if (fabs(difference) > pivot_tolerance) return difference < 0;
    }
    return a < b;
}

static bool pivot(PricingBasis *basis, double *direction, const uint64_t *coverage,
                  size_t words, size_t surplus, const size_t *positions)
{
    size_t rows = basis->rows, leave = rows;
    for (size_t b = 0; b < rows; ++b) {
        double *inverse = basis->inverse + b * basis->stride;
        double value = 0;
        if (surplus < rows) {
            value = -inverse[surplus];
        } else {
            for (size_t r = next_row(coverage, words, 0); r != SIZE_MAX;
                 r = next_row(coverage, words, r + 1)) {
                size_t index = positions ? positions[r] : r;
                if (index < rows) value += inverse[index];
            }
        }
        direction[b] = value;
        if (value > pivot_tolerance &&
            (leave == rows || leaves_first(basis, direction, b, leave))) leave = b;
    }
    if (leave == rows) return false;
    double *leaving = basis->inverse + leave * basis->stride;
    basis->rhs[leave] /= direction[leave];
    for (size_t r = 0; r < rows; ++r) leaving[r] /= direction[leave];
    for (size_t b = 0; b < rows; ++b) {
        if (b == leave) continue;
        basis->rhs[b] -= direction[b] * basis->rhs[leave];
        double *inverse = basis->inverse + b * basis->stride;
        for (size_t r = 0; r < rows; ++r) inverse[r] -= direction[b] * leaving[r];
    }
    basis->surplus[leave] = surplus < rows;
    return true;
}

static void free_basis(NsatBudget *budget, PricingBasis basis)
{
    nsat_free(budget, basis.surplus, basis.rows * sizeof(*basis.surplus));
    nsat_free(budget, basis.rhs, basis.rows * 3 * sizeof(*basis.rhs));
    nsat_free(budget, basis.inverse, basis.rows * basis.rows * sizeof(*basis.inverse));
}

/* The complete family is available for the integer search. The revised
 * simplex proposes weights only. Every bound below uses integer arithmetic. */
static bool price_weights(WeightedCover *w, uint64_t *weights, PricingBasis *saved)
{
    size_t rows = w->input->positives;
    PricingBasis basis = {.rows = rows, .stride = rows};
    basis.inverse = allocate_matrix(&w->budget, rows, rows, sizeof(*basis.inverse));
    basis.rhs = nsat_alloc(&w->budget, rows, 3 * sizeof(*basis.rhs));
    basis.surplus = nsat_alloc(&w->budget, rows, sizeof(*basis.surplus));
    if (!basis.inverse || !basis.rhs || !basis.surplus) goto done;
    if (!nsat_copy(&w->budget, basis.surplus, NULL, rows * sizeof(*basis.surplus))) goto done;
    /* Charge the first row-square pass for setting up the initial basis. */
    if (!nsat_work(&w->budget, (uint64_t)rows * rows)) goto done;
    memset(basis.inverse, 0, rows * rows * sizeof(*basis.inverse));
    double *prices = basis.rhs + rows, *direction = prices + rows;
    for (size_t r = 0; r < rows; ++r) {
        basis.inverse[r * rows + r] = 1;
        basis.rhs[r] = 1;
        weights[r] = 1;
    }

    double best = 0;
    for (unsigned pass = 0; pass < 1024; ++pass) {
        if (!nsat_work(&w->budget, (uint64_t)rows * rows)) break;
        memset(prices, 0, rows * sizeof(*prices));
        for (size_t b = 0; b < rows; ++b) {
            if (basis.surplus[b]) continue;
            for (size_t r = 0; r < rows; ++r) prices[r] += basis.inverse[b * rows + r];
        }
        size_t surplus = rows;
        double total = 0;
        for (size_t r = 0; r < rows; ++r) {
            if (!isfinite(prices[r]) || fabs(prices[r]) > (double)rows) goto done;
            if (prices[r] < -pivot_tolerance && (surplus == rows || prices[r] < prices[surplus]))
                surplus = r;
            if (prices[r] < 0 && prices[r] >= -pivot_tolerance) prices[r] = 0;
            total += prices[r];
        }
        const uint64_t *entering = NULL;
        if (surplus == rows) {
            double maximum = 0;
            for (size_t c = 0; c < w->result.maximal; ++c) {
                if (!nsat_work(&w->budget, w->mask_words)) goto done;
                const uint64_t *coverage = column(w, c);
                double value = 0;
                for (size_t r = next_row(coverage, w->mask_words, 0); r != SIZE_MAX;
                     r = next_row(coverage, w->mask_words, r + 1))
                    value += prices[r];
                if (value > maximum) {
                    maximum = value;
                    entering = coverage;
                }
            }
            if (maximum > 0 && total / maximum > best) {
                best = total / maximum;
                /* Scaling larger charts keeps k times any integer weight
                 * sum in uint64_t, including the per-row rounding error. */
                double scale = rows <= 64 ? 1e9 : 1e9 / fmax(1.0, total);
                for (size_t r = 0; r < rows; ++r)
                    weights[r] = (uint64_t)(prices[r] * scale + 0.5);
            }
            if (maximum <= 1 + pivot_tolerance) break;
        }
        if (!nsat_work(&w->budget, (uint64_t)rows * rows)) break;
        if (!pivot(&basis, direction, entering, w->mask_words, surplus, NULL)) break;
        ++w->result.pivots;
    }
done:
    if (saved) *saved = basis;
    else free_basis(&w->budget, basis);
    return w->budget.status == NSAT_OK;
}

static uint64_t weight(const WeightedCover *w, const uint64_t *a, const uint64_t *b,
                       const uint64_t *local)
{
    uint64_t total = 0;
    for (size_t word = 0; word < w->mask_words; ++word) {
        uint64_t bits = a[word] & b[word];
        if (!bits) continue;
        if (local) {
            for (; bits; bits &= bits - 1)
                total += local[64 * word + (size_t)__builtin_ctzll(bits)];
            continue;
        }
        for (size_t byte = 0; byte < 8; ++byte, bits >>= 8)
            total += w->weight_bytes[word * 8 + byte][bits & 255u];
    }
    return total;
}

static size_t coverage_size(const WeightedCover *w, size_t index)
{
    size_t size = 0;
    const uint64_t *mask = column(w, index);
    for (size_t word = 0; word < w->mask_words; ++word) size += __builtin_popcountll(mask[word]);
    return size;
}

/* Prefer wider columns, with a stable tie break. Counting sort keeps
 * preparation budgeted even when many columns fit in memory. */
static bool order_columns(WeightedCover *w, size_t *order)
{
    size_t rows = w->input->positives, cols = w->result.maximal;
    size_t *offsets = nsat_alloc(&w->budget, rows + 1, sizeof(*offsets));
    bool ok = false;
    if (!offsets || !nsat_copy(&w->budget, offsets, NULL, (rows + 1) * sizeof(*offsets)))
        goto done;
    for (size_t c = 0; c < cols; ++c) {
        if (!nsat_work(&w->budget, w->mask_words)) goto done;
        ++offsets[coverage_size(w, c)];
    }
    for (size_t r = rows, next = 0; r; --r) {
        if (!nsat_work(&w->budget, 1)) goto done;
        size_t count = offsets[r];
        offsets[r] = next;
        next += count;
    }
    for (size_t c = 0; c < cols; ++c) {
        if (!nsat_work(&w->budget, w->mask_words)) goto done;
        order[offsets[coverage_size(w, c)]++] = c;
    }
    ok = true;
done:
    nsat_free(&w->budget, offsets, (rows + 1) * sizeof(*offsets));
    return ok;
}

typedef struct {
    size_t count, cursor, anchor;
    const uint64_t *weights;
    PricingBasis basis;
} SearchFrame;

typedef struct {
    size_t min_branches, levels, square;
    unsigned pivots;
    PricingBasis root;
    double *inverse, *rhs;
    bool *surplus;
    size_t *coordinates, *positions;
    uint64_t *weights;
} Repricing;

static bool allocate_repricing(WeightedCover *w, Repricing *pricing, size_t levels)
{
    size_t rows = w->input->positives;
    pricing->levels = levels;
    if (!nsat_mul(rows, rows, &pricing->square))
        return nsat_stop(&w->budget, NSAT_ERROR, NSAT_STOP_INPUT);
    pricing->inverse = allocate_matrix(&w->budget, levels, pricing->square, sizeof(double));
    pricing->rhs = allocate_matrix(&w->budget, levels, rows, sizeof(double));
    pricing->surplus = allocate_matrix(&w->budget, levels, rows, sizeof(bool));
    pricing->coordinates = allocate_matrix(&w->budget, levels, rows, sizeof(size_t));
    pricing->positions = nsat_alloc(&w->budget, rows, sizeof(size_t));
    pricing->weights = allocate_matrix(&w->budget, levels, rows, sizeof(uint64_t));
    return pricing->inverse && pricing->rhs && pricing->surplus && pricing->coordinates &&
        pricing->positions && pricing->weights;
}

static void free_repricing(WeightedCover *w, Repricing pricing)
{
    size_t cells = pricing.levels * w->input->positives;
    nsat_free(&w->budget, pricing.weights, cells * sizeof(*pricing.weights));
    nsat_free(&w->budget, pricing.positions, w->input->positives * sizeof(*pricing.positions));
    nsat_free(&w->budget, pricing.coordinates, cells * sizeof(*pricing.coordinates));
    nsat_free(&w->budget, pricing.surplus, cells * sizeof(*pricing.surplus));
    nsat_free(&w->budget, pricing.rhs, cells * sizeof(*pricing.rhs));
    nsat_free(&w->budget, pricing.inverse, pricing.levels * pricing.square * sizeof(*pricing.inverse));
    free_basis(&w->budget, pricing.root);
}

/* Removing a row frees one direction in the old primal basis. Move in
 * that direction up to the first zero basic variable. Removing that variable
 * and the covered row gives a smaller feasible basis by elimination. */
static bool project_basis(WeightedCover *w, PricingBasis *basis, size_t coordinate,
                           double *direction)
{
    size_t rows = basis->rows, stride = basis->stride, leave = rows;
    if (!nsat_work(&w->budget, (uint64_t)rows * rows)) return false;
    bool positive = false;
    for (size_t b = 0; b < rows; ++b) {
        direction[b] = basis->inverse[b * stride + coordinate];
        positive |= direction[b] > pivot_tolerance;
    }
    if (!positive)
        for (size_t b = 0; b < rows; ++b) direction[b] = -direction[b];
    for (size_t b = 0; b < rows; ++b) {
        if (direction[b] > pivot_tolerance &&
            (leave == rows || leaves_first(basis, direction, b, leave))) leave = b;
    }
    if (leave == rows) return false;
    const double *leaving = basis->inverse + leave * stride;
    size_t last = rows - 1;
    for (size_t b = 0; b < rows; ++b) {
        if (b == leave) continue;
        double *row = basis->inverse + b * stride;
        double factor = direction[b] / direction[leave];
        basis->rhs[b] -= factor * basis->rhs[leave];
        for (size_t c = 0; c < rows; ++c) row[c] -= factor * leaving[c];
        row[coordinate] = row[last];
    }
    if (leave != last) {
        memcpy(basis->inverse + leave * stride, basis->inverse + last * stride,
               last * sizeof(*basis->inverse));
        basis->rhs[leave] = basis->rhs[last];
        basis->surplus[leave] = basis->surplus[last];
    }
    basis->coordinates[coordinate] = basis->coordinates[last];
    basis->rows = last;
    return true;
}

static double cube_price(const WeightedCover *w, size_t cube, const double *prices,
                          const size_t *positions, size_t rows)
{
    const uint64_t *mask = column(w, cube);
    double value = 0;
    for (size_t r = next_row(mask, w->mask_words, 0); r != SIZE_MAX;
         r = next_row(mask, w->mask_words, r + 1)) {
        size_t index = positions[r];
        if (index < rows) value += prices[index];
    }
    return value;
}

/* A basic column filtered out of this node can leave the projected LP
 * but cannot re-enter. We test the suggested weights with integer arithmetic
 * on the actual surviving columns before they can affect the search. */
static uint64_t *reprice_node(WeightedCover *w, Repricing *pricing, SearchFrame *frame,
                             const uint64_t *remaining, const size_t *choices, size_t count,
                             size_t depth, int slots)
{
    size_t stride = w->input->positives, offset = (depth - 1) * stride;
    PricingBasis parent = frame->basis;
    size_t rows = parent.rows;
    PricingBasis basis = {.rows = rows, .stride = stride,
        .inverse = pricing->inverse + (depth - 1) * pricing->square,
        .rhs = pricing->rhs + offset, .surplus = pricing->surplus + offset,
        .coordinates = pricing->coordinates + offset};
    if (!nsat_work(&w->budget, (uint64_t)rows * rows + 3 * rows)) return NULL;
    for (size_t b = 0; b < rows; ++b) {
        memcpy(basis.inverse + b * stride, parent.inverse + b * parent.stride,
               rows * sizeof(*basis.inverse));
        basis.coordinates[b] = parent.coordinates ? parent.coordinates[b] : b;
    }
    memcpy(basis.rhs, parent.rhs, rows * sizeof(*basis.rhs));
    memcpy(basis.surplus, parent.surplus, rows * sizeof(*basis.surplus));

    double *prices = pricing->root.rhs + stride, *direction = prices + stride;
    ++w->result.repriced;
    for (size_t i = 0; i < basis.rows;) {
        if (contains_row(remaining, basis.coordinates[i])) {
            ++i;
        } else if (!project_basis(w, &basis, i, direction)) {
            return NULL;
        }
    }
    rows = basis.rows;
    for (size_t i = 0; i < stride; ++i) pricing->positions[i] = SIZE_MAX;
    for (size_t i = 0; i < rows; ++i) pricing->positions[basis.coordinates[i]] = i;
    uint64_t *weights = pricing->weights + offset;
    double best = 0;
    for (unsigned pass = 0;; ++pass) {
        if (!nsat_work(&w->budget, (uint64_t)rows * rows)) return NULL;
        memset(prices, 0, rows * sizeof(*prices));
        for (size_t b = 0; b < rows; ++b) {
            if (basis.surplus[b]) continue;
            for (size_t r = 0; r < rows; ++r) prices[r] += basis.inverse[b * stride + r];
        }
        size_t surplus = rows;
        double total = 0;
        for (size_t r = 0; r < rows; ++r) {
            if (!isfinite(prices[r]) || fabs(prices[r]) > (double)stride) goto done;
            if (prices[r] < -pivot_tolerance && (surplus == rows || prices[r] < prices[surplus]))
                surplus = r;
            if (prices[r] < 0 && prices[r] >= -pivot_tolerance) prices[r] = 0;
            total += prices[r];
        }
        size_t entering = 0;
        if (surplus == rows) {
            double largest = 0;
            for (size_t c = 0; c < count; ++c) {
                if (!nsat_work(&w->budget, w->mask_words)) return NULL;
                double value = cube_price(w, choices[c], prices, pricing->positions, rows);
                if (value <= largest) continue;
                largest = value;
                entering = choices[c];
            }
            if (largest > 0 && total / largest > best) {
                best = total / largest;
                double scale = stride <= 64 ? 1e9 : 1e9 / fmax(1.0, total);
                for (size_t i = 0; i < rows; ++i)
                    weights[basis.coordinates[i]] = (uint64_t)(prices[i] * scale + 0.5);
            }
            if (largest <= 1 + pivot_tolerance || best > slots + 1e-7) break;
        }
        if (pass == pricing->pivots) break;
        if (!nsat_work(&w->budget, (uint64_t)rows * rows)) return NULL;
        const uint64_t *coverage = surplus < rows ? NULL : column(w, entering);
        if (!pivot(&basis, direction, coverage, w->mask_words, surplus, pricing->positions)) break;
        ++w->result.reprice_pivots;
    }
done:
    frame->basis = basis;
    return best > 0 ? weights : NULL;
}

/* Every cube of a k-cube completion must cover at least W-(k-1)U. */
static bool filter_choices(WeightedCover *w, const uint64_t *remaining, int slots,
                           size_t *child, size_t count, SearchFrame *frame, size_t *frequency,
                           uint64_t need, uint64_t maximum)
{
    uint64_t other = (uint64_t)(slots - 1) * maximum;
    uint64_t minimum = need > other ? need - other : 0;
    size_t rows = w->input->positives;
    if (!nsat_work(&w->budget, rows)) return false;
    memset(frequency, 0, rows * sizeof(*frequency));
    size_t eligible = 0;
    for (size_t c = 0; c < count; ++c) {
        if (!nsat_work(&w->budget, w->mask_words)) return false;
        const uint64_t *mask = column(w, child[c]);
        if (weight(w, remaining, mask, frame->weights) < minimum) continue;
        child[eligible++] = child[c];
        for (size_t word = 0; word < w->mask_words; ++word) {
            for (uint64_t trace = mask[word] & remaining[word]; trace; trace &= trace - 1)
                ++frequency[64 * word + (size_t)__builtin_ctzll(trace)];
        }
    }
    size_t anchor = next_row(remaining, w->mask_words, 0);
    for (size_t r = anchor; r != SIZE_MAX; r = next_row(remaining, w->mask_words, r + 1))
        if (frequency[r] < frequency[anchor]) anchor = r;
    frame->count = eligible;
    frame->anchor = anchor;
    return true;
}

static NsatStatus prepare_search(WeightedCover *w, const uint64_t *remaining, int slots,
                                 const size_t *parent, size_t count, size_t *child,
                                 SearchFrame *frame, size_t *frequency, size_t *chosen,
                                 Repricing *pricing, size_t depth)
{
    if (!nsat_work(&w->budget, 1)) return w->budget.status;
    ++w->result.nodes;
    *frame = (SearchFrame){0};
    if (pricing->pivots) {
        frame->weights = depth ? frame[-1].weights : NULL;
        frame->basis = depth ? frame[-1].basis : pricing->root;
    }
    size_t nonempty = 0;
    uint64_t maximum = 0;
    for (size_t c = 0; c < count; ++c) {
        if (!nsat_work(&w->budget, w->mask_words)) return w->budget.status;
        const uint64_t *mask = column(w, parent[c]);
        bool complete = true;
        uint64_t any = 0;
        for (size_t word = 0; word < w->mask_words; ++word) {
            complete &= !(remaining[word] & ~mask[word]);
            any |= remaining[word] & mask[word];
        }
        if (complete) {
            chosen[depth] = parent[c];
            return NSAT_SAT;
        }
        if (!any) continue;
        uint64_t value = weight(w, remaining, mask, frame->weights);
        if (value > maximum) maximum = value;
        child[nonempty++] = parent[c];
    }
    uint64_t need = weight(w, remaining, remaining, frame->weights);
    if (!nonempty || slots == 1 || need > (uint64_t)slots * maximum) return NSAT_UNSAT;
    if (!filter_choices(w, remaining, slots, child, nonempty, frame, frequency, need, maximum))
        return w->budget.status;
    if (!pricing->pivots || !depth || frequency[frame->anchor] < pricing->min_branches)
        return NSAT_OK;
    uint64_t *local = reprice_node(w, pricing, frame, remaining, child, frame->count,
        depth, slots);
    if (w->budget.status != NSAT_OK) return w->budget.status;
    if (!local) return NSAT_OK;
    uint64_t local_total = weight(w, remaining, remaining, local), local_maximum = 0;
    for (size_t c = 0; c < frame->count; ++c) {
        if (!nsat_work(&w->budget, w->mask_words)) return w->budget.status;
        uint64_t value = weight(w, remaining, column(w, child[c]), local);
        if (value > local_maximum) local_maximum = value;
    }
    if (local_total > (uint64_t)slots * local_maximum) {
        ++w->result.reprice_pruned;
        return NSAT_UNSAT;
    }
    if (local_maximum && (!maximum ||
        (long double)local_total / local_maximum >= (long double)need / maximum)) {
        frame->weights = local;
        if (!filter_choices(w, remaining, slots, child, frame->count, frame, frequency,
                            local_total, local_maximum)) return w->budget.status;
        if (!frequency[frame->anchor]) ++w->result.reprice_pruned;
    }
    return NSAT_OK;
}

static bool dominated(WeightedCover *w, const uint64_t *remaining, const size_t *choices,
                       size_t count, size_t c)
{
    const uint64_t *mask = column(w, choices[c]);
    for (size_t d = 0; d < count; ++d) {
        if (!nsat_work(&w->budget, w->mask_words)) return false;
        if (d == c) continue;
        const uint64_t *other = column(w, choices[d]);
        bool equal = true;
        size_t word = 0;
        for (; word < w->mask_words; ++word) {
            uint64_t trace = mask[word] & remaining[word], all = other[word] & remaining[word];
            if (trace & ~all) break;
            equal &= trace == all;
        }
        if (word == w->mask_words && (!equal || d < c)) return true;
    }
    return false;
}

/* Heap-backed frames let resource budgets govern depth. Each child sees
 * only candidates that could appear in a solution of its parent's target. */
static NsatStatus search(WeightedCover *w, int slots, const size_t *order,
                          size_t *scratch, uint64_t *remaining, SearchFrame *frames,
                          size_t *frequency, size_t *chosen, int *terms, Repricing *pricing)
{
    const size_t cols = w->result.maximal;
    all_rows(w, remaining);
    bool entering = true;
    for (size_t depth = 0;;) {
        uint64_t *rows = remaining + depth * w->mask_words;
        size_t *choices = scratch + depth * cols;
        SearchFrame *frame = frames + depth;
        if (entering) {
            const size_t *parent = depth ? choices - cols : order;
            size_t count = depth ? frames[depth - 1].count : cols;
            NsatStatus status = prepare_search(w, rows, slots - (int)depth, parent, count,
                choices, frame, frequency, chosen, pricing, depth);
            if (status == NSAT_SAT) {
                *terms = (int)depth + 1;
                return status;
            }
            if (status == NSAT_UNSAT) frame->count = frame->cursor = 0;
            else if (status != NSAT_OK) return status;
        }
        entering = false;
        while (frame->cursor < frame->count) {
            if (!nsat_work(&w->budget, 1)) return w->budget.status;
            size_t c = frame->cursor++;
            const uint64_t *mask = column(w, choices[c]);
            if (!contains_row(mask, frame->anchor) || dominated(w, rows, choices, frame->count, c))
                continue;
            if (!nsat_work(&w->budget, w->mask_words)) return w->budget.status;
            uint64_t *next = rows + w->mask_words;
            for (size_t word = 0; word < w->mask_words; ++word) next[word] = rows[word] & ~mask[word];
            chosen[depth++] = choices[c];
            entering = true;
            break;
        }
        if (entering) continue;
        if (w->budget.status != NSAT_OK) return w->budget.status;
        if (!depth) return NSAT_UNSAT;
        --depth;
    }
}

static int *materialize(WeightedCover *w, const size_t *chosen, int terms)
{
    const NsatCoverInput *q = w->input;
    size_t cells = q->inputs * (size_t)terms;
    int *cover = nsat_alloc(&w->budget, cells, sizeof(*cover));
    if (!cover) return NULL;
    for (int c = 0; c < terms; ++c) {
        const uint64_t *mask = column(w, chosen[c]);
        size_t first = next_row(mask, w->mask_words, 0);
        for (size_t i = 0; i < q->inputs; ++i) {
            int value = q->on[first * q->inputs + i];
            for (size_t row = first; row != SIZE_MAX; row = next_row(mask, w->mask_words, row + 1)) {
                if (!nsat_work(&w->budget, 1)) {
                    nsat_free(&w->budget, cover, cells * sizeof(*cover));
                    return NULL;
                }
                if (q->on[row * q->inputs + i] == value) continue;
                value = 0;
                break;
            }
            cover[i * (size_t)terms + (size_t)c] = value;
        }
    }
    return cover;
}

static void verify_bound(WeightedCover *w, const uint64_t *weights,
                         uint64_t maximum, int lower)
{
    if (w->result.lower == lower) return;
    NsatLimits limits = w->budget.limits;
    NsatStatus status = NSAT_UNKNOWN;
    if (limits.work && w->budget.stats.work >= limits.work) {
        w->result.check.stop = NSAT_STOP_WORK;
    } else if (limits.memory_bytes && w->budget.stats.memory_bytes >= limits.memory_bytes) {
        w->result.check.stop = NSAT_STOP_MEMORY;
    } else {
        if (limits.work) limits.work -= w->budget.stats.work;
        if (limits.memory_bytes) limits.memory_bytes -= w->budget.stats.memory_bytes;
        status = ipc_verify_weighted_bound(w->input, weights, w->input->positives,
            maximum, w->result.lower, limits, &w->result.check);
    }
    w->result.check_status = status;
    if (status == NSAT_UNSAT) return;
    w->result.lower = lower;
    w->result.status = status;
    w->result.stop = w->result.check.stop;
}

static IpcWeightedCoverResult weighted_cover(const NsatCoverInput *q, int lower,
                                              NsatLimits limits, size_t min_branches,
                                              unsigned pivots, bool seed_first)
{
    WeightedCover w = {.budget = {.limits = limits}, .input = q,
                       .result = {.lower = lower, .status = NSAT_UNKNOWN}};
    if (!q || !q->on || !q->off || !q->inputs ||
        !q->positives || q->positives > INT_MAX || !q->negatives ||
        lower < 1 || q->target <= (size_t)lower || q->target > q->positives)
        return w.result;
    w.result.supported = true;
    size_t *order = NULL, *scratch = NULL, *frequency = NULL, *chosen = NULL;
    uint64_t *remaining = NULL, *weights = NULL;
    uint64_t maximum = 0;
    SearchFrame *frames = NULL;
    Repricing pricing = {.min_branches = min_branches, .pivots = min_branches ? pivots : 0};
    Repricing fixed = {0};
    size_t diff_count = 0, levels = q->target - 1;
    w.input_words = (q->inputs - 1) / 64 + 1;
    w.mask_words = (q->positives - 1) / 64 + 1;
    if (q->inputs > SIZE_MAX / sizeof(int) / q->positives ||
        q->inputs > SIZE_MAX / sizeof(int) / q->negatives ||
        !nsat_mul(q->negatives, w.input_words, &w.row_words) ||
        w.row_words > SIZE_MAX / sizeof(uint64_t) ||
        !nsat_mul(q->positives, w.row_words, &diff_count)) {
        nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_INPUT);
        goto done;
    }
    w.diff = nsat_alloc(&w.budget, diff_count, sizeof(*w.diff));
    if (!w.diff || !prepare_differences(&w) || !enumerate(&w)) goto done;
    nsat_free(&w.budget, w.diff, diff_count * sizeof(*w.diff));
    w.diff = NULL;

    const size_t padded_rows = w.mask_words * 64;
    weights = nsat_alloc(&w.budget, padded_rows, sizeof(*weights));
    w.weight_bytes = nsat_alloc(&w.budget, w.mask_words * 8, sizeof(*w.weight_bytes));
    if (!weights || !w.weight_bytes ||
        !nsat_copy(&w.budget, weights + q->positives, NULL,
                   (padded_rows - q->positives) * sizeof(*weights))) goto done;
    if (!price_weights(&w, weights, pricing.pivots ? &pricing.root : NULL)) goto done;
    uint64_t total = 0;
    for (size_t byte = 0; byte < w.mask_words * 8; ++byte) {
        if (!nsat_work(&w.budget, 256)) goto done;
        w.weight_bytes[byte][0] = 0;
        for (unsigned bits = 1; bits < 256; ++bits) {
            w.weight_bytes[byte][bits] = w.weight_bytes[byte][bits & (bits - 1)] +
                weights[byte * 8 + (unsigned)__builtin_ctz(bits)];
        }
        total += w.weight_bytes[byte][255];
    }
    for (size_t c = 0; c < w.result.maximal; ++c) {
        if (!nsat_work(&w.budget, w.mask_words)) goto done;
        const uint64_t *mask = column(&w, c);
        uint64_t value = weight(&w, mask, mask, NULL);
        if (value > maximum) maximum = value;
    }
    int incumbent = (int)q->target;
    if (maximum) {
        uint64_t bound = total / maximum + (total % maximum != 0);
        if (bound > q->target) {
            nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_PROOF);
            goto done;
        }
        if ((int)bound > w.result.lower) w.result.lower = (int)bound;
    }
    if (w.result.lower < incumbent) {
        order = nsat_alloc(&w.budget, w.result.maximal, sizeof(*order));
        scratch = allocate_matrix(&w.budget, levels, w.result.maximal, sizeof(*scratch));
        remaining = allocate_matrix(&w.budget, levels, w.mask_words, sizeof(*remaining));
        frames = nsat_alloc(&w.budget, levels, sizeof(*frames));
        chosen = nsat_alloc(&w.budget, levels, sizeof(*chosen));
        frequency = nsat_alloc(&w.budget, q->positives, sizeof(*frequency));
        if (!order || !scratch || !remaining || !frames || !chosen || !frequency ||
            !order_columns(&w, order)) goto done;
    }
    while (w.result.lower < incumbent) {
        Repricing *active = seed_first && !w.result.first_terms ? &fixed : &pricing;
        if (active->pivots && incumbent > 2 && !active->levels &&
            !allocate_repricing(&w, active, (size_t)incumbent - 2)) break;
        int terms = 0;
        NsatStatus decision = search(&w, incumbent - 1, order, scratch, remaining, frames,
                                     frequency, chosen, &terms, active);
        if (!nsat_poll(&w.budget)) break;
        if (decision == NSAT_UNSAT) {
            w.result.lower = incumbent;
            break;
        }
        if (decision != NSAT_SAT || terms < w.result.lower || terms >= incumbent) {
            nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_PROOF);
            break;
        }
        int *cover = materialize(&w, chosen, terms);
        if (!cover) break;
        if (!w.result.first_terms) {
            w.result.first_terms = terms;
            w.result.first_cover_work = w.budget.stats.work;
            w.result.first_cover_nodes = w.result.nodes;
        }
        nsat_free(&w.budget, w.result.cover, q->inputs * (size_t)w.result.terms * sizeof(*cover));
        w.result.cover = cover;
        w.result.terms = incumbent = terms;
    }
    if (w.result.lower == incumbent && nsat_poll(&w.budget)) w.result.status = NSAT_OK;
done:
    if (w.budget.status == NSAT_ERROR) {
        nsat_free(&w.budget, w.result.cover, q->inputs * (size_t)w.result.terms * sizeof(int));
        w.result.cover = NULL;
        w.result.terms = 0;
        w.result.lower = lower;
    }
    if (w.budget.status != NSAT_OK) w.result.status = w.budget.status;
    w.result.stop = w.budget.stats.stop;
    nsat_free(&w.budget, order, w.result.maximal * sizeof(*order));
    nsat_free(&w.budget, scratch, levels * w.result.maximal * sizeof(*scratch));
    nsat_free(&w.budget, remaining, levels * w.mask_words * sizeof(*remaining));
    nsat_free(&w.budget, frames, levels * sizeof(*frames));
    nsat_free(&w.budget, chosen, levels * sizeof(*chosen));
    nsat_free(&w.budget, frequency, q->positives * sizeof(*frequency));
    nsat_free(&w.budget, w.columns, w.capacity * w.mask_words * sizeof(*w.columns));
    nsat_free(&w.budget, w.diff, diff_count * sizeof(*w.diff));
    nsat_free(&w.budget, w.weight_bytes, w.mask_words * 8 * sizeof(*w.weight_bytes));
    free_repricing(&w, pricing);
    verify_bound(&w, weights, maximum, lower);
    nsat_free(&w.budget, weights, w.mask_words * 64 * sizeof(*weights));
    nsat_release_budget(&w.budget);
    return w.result;
}

IpcWeightedCoverResult ipc_weighted_cover(const NsatCoverInput *q, int lower,
                                          NsatLimits limits)
{
    return weighted_cover(q, lower, limits, 0, 0, false);
}

IpcWeightedCoverResult ipc_weighted_cover_reprice(const NsatCoverInput *q, int lower,
                                                  NsatLimits limits, size_t min_branches,
                                                  unsigned pivots)
{
    return weighted_cover(q, lower, limits, min_branches, pivots, false);
}

IpcWeightedCoverResult ipc_weighted_cover_seed_then_reprice(const NsatCoverInput *q,
                                                            int lower, NsatLimits limits,
                                                            size_t min_branches, unsigned pivots)
{
    return weighted_cover(q, lower, limits, min_branches, pivots, true);
}
