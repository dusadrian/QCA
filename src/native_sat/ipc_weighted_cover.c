/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#include "ipc_weighted_cover.h"
#include "ipc_weighted_verify.h"
#include "native_sat_support.h"

#include <math.h>
#include "qca_threads.h"
#if defined(HAVE_PTHREAD)
#include <stdatomic.h>
#include <time.h>
#endif

typedef struct ChartWorker ChartWorker;

typedef struct {
    NsatBudget budget;
    ChartWorker *worker;
    IpcWeightedCoverResult result;
    const NsatCoverInput *input;
    uint64_t *diff, *columns;
    size_t capacity, input_words, row_words, mask_words;
    uint64_t (*weight_bytes)[256];
    NsatLimits caller_limits;
    uint64_t probe_work, next_probe, progress, checkpoint_progress, extensions;
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
    size_t next = 0;
    for (size_t r = rows;; --r) {
        if (!nsat_work(&w->budget, 1)) goto done;
        size_t count = offsets[r];
        offsets[r] = next;
        next += count;
        if (!r) break;
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

static bool chart_claim(ChartWorker *worker, size_t depth);
static void chart_closed(ChartWorker *worker, size_t depth);

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
            NsatStatus status = chart_claim(w->worker, depth)
                ? prepare_search(w, rows, slots - (int)depth, parent, count,
                    choices, frame, frequency, chosen, pricing, depth)
                : NSAT_UNSAT;
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
            /* This branch selects the column. Later siblings exclude it,
             * so a cover containing several alternatives for the anchor
             * row is not explored again in different selection orders.
             * It covers no remaining row in its own child either. */
            size_t following = frame->count - c - 1;
            if (!nsat_work(&w->budget, following)) return w->budget.status;
            memmove(choices + c, choices + c + 1, following * sizeof(*choices));
            --frame->count;
            --frame->cursor;
            entering = true;
            break;
        }
        if (entering) continue;
        if (w->budget.status != NSAT_OK) return w->budget.status;
        chart_closed(w->worker, depth);
        if (!depth) return NSAT_UNSAT;
        --depth;
        if (!depth) ++w->progress;
    }
}

#if defined(HAVE_PTHREAD)
typedef struct {
    pthread_mutex_t mutex;
    pthread_cond_t changed;
    atomic_bool stop;
    atomic_size_t next_task;
    atomic_uint_fast64_t completed;
    uint64_t available, future, probe, checkpoint, extensions;
    unsigned waiting, finished, count;
    NsatStopReason reason;
} ChartPool;

struct ChartWorker {
    WeightedCover cover;
    Repricing pricing;
    ChartPool *pool;
    const size_t *order;
    size_t *scratch, *frequency, *chosen;
    uint64_t *remaining;
    SearchFrame *frames;
    size_t task, split;
    uint64_t margin;
    bool owns_task;
    int slots, terms;
    NsatStatus status;
};

/* Each worker walks the same small, deterministic prefix. A task number is
 * claimed once, immediately before searching its subtree. Sibling exclusions
 * make these subtrees disjoint. Walking the prefix again needs no task snapshots
 * and lets an idle worker take the next branch, regardless of its depth-first
 * predecessor's cost. No bound or incumbent depends on another worker. */
static bool chart_claim(ChartWorker *worker, size_t depth)
{
    if (!worker || depth != worker->split) return true;
    size_t task = worker->task++;
    worker->owns_task = atomic_compare_exchange_strong(&worker->pool->next_task,
        &task, worker->task);
    return worker->owns_task;
}

static void chart_closed(ChartWorker *worker, size_t depth)
{
    if (worker && depth == worker->split && worker->owns_task) {
        atomic_fetch_add(&worker->pool->completed, 1);
        worker->owns_task = false;
    }
}

/* Reserve work in chunks rather than locking every scan. The sum of grants
 * never exceeds the one caller allowance. Unused work is returned on exit.
 * Workers wait at a probe with their search state intact; only the coordinator
 * may extend it. A small unspent tail can cause an early UNKNOWN, never a proof. */
static bool chart_worker_poll(void *context)
{
    ChartWorker *worker = context;
    ChartPool *pool = worker->pool;
    NsatBudget *budget = &worker->cover.budget;
    if (atomic_load(&pool->stop)) return true;
    if (budget->limits.work - budget->stats.work > worker->margin) return false;
    pthread_mutex_lock(&pool->mutex);
    while (!pool->available && !atomic_load(&pool->stop)) {
        ++pool->waiting;
        pthread_cond_broadcast(&pool->changed);
        pthread_cond_wait(&pool->changed, &pool->mutex);
        --pool->waiting;
    }
    if (!atomic_load(&pool->stop)) {
        uint64_t grant = UINT64_C(1048576);
        if (grant < worker->margin * 2) grant = worker->margin * 2;
        if (grant > pool->available) grant = pool->available;
        pool->available -= grant;
        budget->limits.work += grant;
    }
    bool stopped = atomic_load(&pool->stop);
    pthread_mutex_unlock(&pool->mutex);
    return stopped;
}

static void *chart_worker_run(void *context)
{
    ChartWorker *worker = context;
    ChartPool *pool = worker->pool;
    /* The coordinator holds the mutex until every successfully started worker
     * is counted. A failed pthread_create simply reduces the worker count. */
    pthread_mutex_lock(&pool->mutex);
    pthread_mutex_unlock(&pool->mutex);
    if (nsat_poll(&worker->cover.budget)) {
        worker->status = search(&worker->cover, worker->slots, worker->order,
            worker->scratch, worker->remaining, worker->frames, worker->frequency,
            worker->chosen, &worker->terms, &worker->pricing);
    } else {
        worker->status = worker->cover.budget.status;
    }
    pthread_mutex_lock(&pool->mutex);
    pool->available += worker->cover.budget.limits.work - worker->cover.budget.stats.work;
    ++pool->finished;
    if (worker->status == NSAT_SAT) atomic_store(&pool->stop, true);
    else if (worker->status != NSAT_UNSAT && !atomic_load(&pool->stop)) {
        pool->reason = worker->cover.budget.stats.stop;
        atomic_store(&pool->stop, true);
    }
    pthread_cond_broadcast(&pool->changed);
    pthread_mutex_unlock(&pool->mutex);
    return NULL;
}

static void free_chart_worker(WeightedCover *owner, ChartWorker *worker, size_t levels)
{
    size_t rows = owner->input->positives, cols = owner->result.maximal;
    free_repricing(owner, worker->pricing);
    nsat_free(&owner->budget, worker->frames, levels * sizeof(SearchFrame));
    nsat_free(&owner->budget, worker->chosen, levels * sizeof(size_t));
    nsat_free(&owner->budget, worker->frequency, rows * sizeof(size_t));
    nsat_free(&owner->budget, worker->scratch, levels * cols * sizeof(size_t));
    nsat_free(&owner->budget, worker->remaining, levels * owner->mask_words * sizeof(uint64_t));
}

static bool prepare_chart_worker(WeightedCover *owner, ChartWorker *worker,
    const Repricing *pricing, size_t levels)
{
    size_t rows = owner->input->positives, cols = owner->result.maximal;
    worker->cover = *owner;
    worker->cover.worker = worker;
    worker->cover.result.nodes = worker->cover.result.repriced = 0;
    worker->cover.budget = (NsatBudget){.limits = {
        .poll = chart_worker_poll, .context = worker}};
    worker->pricing.min_branches = pricing->min_branches;
    worker->pricing.pivots = pricing->pivots;
    worker->frames = nsat_alloc(&owner->budget, levels, sizeof(SearchFrame));
    worker->chosen = nsat_alloc(&owner->budget, levels, sizeof(size_t));
    worker->frequency = nsat_alloc(&owner->budget, rows, sizeof(size_t));
    worker->scratch = allocate_matrix(&owner->budget, levels, cols, sizeof(size_t));
    worker->remaining = allocate_matrix(&owner->budget, levels, owner->mask_words, sizeof(uint64_t));
    if (!worker->frames || !worker->chosen || !worker->frequency ||
        !worker->scratch || !worker->remaining) return false;
    if (pricing->pivots) {
        PricingBasis *root = &worker->pricing.root;
        root->rows = root->stride = rows;
        root->inverse = allocate_matrix(&owner->budget, rows, rows, sizeof(double));
        root->rhs = nsat_alloc(&owner->budget, rows, 3 * sizeof(double));
        root->surplus = nsat_alloc(&owner->budget, rows, sizeof(bool));
        if (!root->inverse || !root->rhs || !root->surplus ||
            !nsat_copy(&owner->budget, root->inverse, pricing->root.inverse, rows * rows * sizeof(double)) ||
            !nsat_copy(&owner->budget, root->rhs, pricing->root.rhs, rows * sizeof(double)) ||
            !nsat_copy(&owner->budget, root->surplus, pricing->root.surplus, rows * sizeof(bool))) return false;
        if (levels > 2 && !allocate_repricing(owner, &worker->pricing, levels - 2)) return false;
    }
    /* Largest indivisible search charge, plus the polling interval. */
    worker->margin = rows * rows + cols + 8192;
    return true;
}

static NsatStatus parallel_chart_search(WeightedCover *w, int slots, const size_t *order,
    size_t *scratch, uint64_t *remaining, SearchFrame *frames, size_t *frequency,
    size_t *chosen, int *terms, Repricing *pricing, unsigned requested,
    unsigned *used, uint64_t *tasks)
{
    const size_t levels = (size_t)slots + 1, rows = w->input->positives;
    const size_t cols = w->result.maximal;
    /* The production adapter restricts charts to 128 x 1024. Keep arithmetic
     * and optional workspace sizing bounded for other callers too. */
    if (requested < 2 || rows > 128 || cols > 1024 || levels > 1024)
        return search(w, slots, order, scratch, remaining, frames, frequency, chosen, terms, pricing);
    if (requested > 8) requested = 8;
    size_t bytes = sizeof(ChartWorker) + levels * (sizeof(SearchFrame) +
        sizeof(size_t) + cols * sizeof(size_t) + w->mask_words * sizeof(uint64_t)) + rows * sizeof(size_t);
    if (pricing->pivots) {
        bytes += rows * rows * sizeof(double) + rows * (3 * sizeof(double) + sizeof(bool));
        if (levels > 2) bytes += (levels - 2) * (rows * rows * sizeof(double) +
            rows * (sizeof(double) + sizeof(bool) + sizeof(size_t) + sizeof(uint64_t))) + rows * sizeof(size_t);
    }
    size_t available = nsat_memory_available(&w->budget);
    while (requested > 1 && bytes > available / requested) --requested;
    if (requested < 2)
        return search(w, slots, order, scratch, remaining, frames, frequency, chosen, terms, pricing);

    ChartWorker *workers = nsat_alloc(&w->budget, requested, sizeof(*workers));
    if (!workers) return w->budget.status;
    memset(workers, 0, requested * sizeof(*workers));
    ChartPool pool = {0};
    atomic_init(&pool.stop, false);
    atomic_init(&pool.next_task, 0);
    atomic_init(&pool.completed, 0);
    NsatStatus status = NSAT_UNKNOWN;
    unsigned prepared = 0;
    for (; prepared < requested; ++prepared) {
        ChartWorker *worker = workers + prepared;
        worker->pool = &pool;
        worker->slots = slots;
        worker->split = slots > 2 ? 2 : 1;
        worker->order = order;
        if (!prepare_chart_worker(w, worker, pricing, levels)) {
            ++prepared;
            goto cleanup;
        }
    }
    if (pthread_mutex_init(&pool.mutex, NULL) != 0) goto serial;
    if (pthread_cond_init(&pool.changed, NULL) != 0) {
        pthread_mutex_destroy(&pool.mutex);
        goto serial;
    }
    uint64_t hard = w->budget.limits.work ? w->budget.limits.work : UINT64_MAX;
    if (w->budget.limits.shared && w->budget.limits.shared->work_limit) {
        uint64_t left = w->budget.limits.shared->work_limit - w->budget.limits.shared->work;
        if (left < hard - w->budget.stats.work) hard = w->budget.stats.work + left;
    }
    pool.available = hard - w->budget.stats.work;
    pool.probe = w->probe_work;
    if (pool.probe && w->next_probe < hard) {
        pool.available = w->next_probe > w->budget.stats.work
            ? w->next_probe - w->budget.stats.work : 0;
        pool.future = hard - w->budget.stats.work - pool.available;
    }
    pthread_t threads[8];
    pthread_mutex_lock(&pool.mutex);
    for (unsigned i = 0; i < requested; ++i) {
        if (pthread_create(threads + pool.count, NULL, chart_worker_run, workers + pool.count) != 0) break;
        ++pool.count;
    }
    if (pool.count > *used) *used = pool.count;
    while (pool.finished < pool.count) {
        if (!atomic_load(&pool.stop) && pool.waiting + pool.finished == pool.count && !pool.available) {
            uint64_t completed = atomic_load(&pool.completed);
            if (pool.future && completed > pool.checkpoint) {
                uint64_t grant = pool.future < pool.probe ? pool.future : pool.probe;
                pool.available += grant;
                pool.future -= grant;
                pool.checkpoint = completed;
                ++pool.extensions;
            } else {
                pool.reason = NSAT_STOP_WORK;
                atomic_store(&pool.stop, true);
            }
            pthread_cond_broadcast(&pool.changed);
        }
        pthread_mutex_unlock(&pool.mutex);
        /* Caller callbacks (including R interruption checks) run only here. */
        bool cancelled = w->caller_limits.poll && w->caller_limits.poll(w->caller_limits.context);
        pthread_mutex_lock(&pool.mutex);
        if (cancelled) {
            pool.reason = NSAT_STOP_CANCELLED;
            atomic_store(&pool.stop, true);
            pthread_cond_broadcast(&pool.changed);
        }
        struct timespec deadline;
        timespec_get(&deadline, TIME_UTC);
        deadline.tv_nsec += 10000000;
        if (deadline.tv_nsec >= 1000000000) {
            ++deadline.tv_sec;
            deadline.tv_nsec -= 1000000000;
        }
        if (pool.finished < pool.count) pthread_cond_timedwait(&pool.changed, &pool.mutex, &deadline);
    }
    pthread_mutex_unlock(&pool.mutex);
    for (unsigned i = 0; i < pool.count; ++i) pthread_join(threads[i], NULL);
    pthread_cond_destroy(&pool.changed);
    pthread_mutex_destroy(&pool.mutex);
    if (!pool.count) goto serial;
    status = NSAT_UNSAT;
    int witness = -1;
    uint64_t work = 0;
    for (unsigned i = 0; i < pool.count; ++i) {
        ChartWorker *worker = workers + i;
        work += worker->cover.budget.stats.work;
        w->result.nodes += worker->cover.result.nodes;
        w->result.repriced += worker->cover.result.repriced;
        if (worker->status == NSAT_SAT && witness < 0) witness = (int)i;
        else if (worker->status != NSAT_UNSAT) status = NSAT_UNKNOWN;
    }
    /* Reservations already checked both limits; account actual work once. */
    w->budget.stats.work += work;
    if (w->budget.limits.shared) w->budget.limits.shared->work += work;
    w->extensions += pool.extensions;
    uint64_t completed = atomic_load(&pool.completed);
    w->progress += completed;
    w->checkpoint_progress = w->progress;
    *tasks += completed;
    w->next_probe = hard - pool.future;
    if (witness >= 0) {
        *terms = workers[witness].terms;
        memcpy(chosen, workers[witness].chosen, (size_t)*terms * sizeof(*chosen));
        status = NSAT_SAT;
    } else if (status != NSAT_UNSAT) {
        nsat_stop(&w->budget, NSAT_UNKNOWN, pool.reason ? pool.reason : NSAT_STOP_CANCELLED);
    }
    goto cleanup;
serial:
    /* Thread creation/support failure leaves the exact serial path available. */
    status = NSAT_OK;
cleanup:
    for (unsigned i = 0; i < prepared; ++i) free_chart_worker(w, workers + i, levels);
    nsat_free(&w->budget, workers, requested * sizeof(*workers));
    if (w->budget.status != NSAT_OK) return w->budget.status;
    if (status == NSAT_OK)
        return search(w, slots, order, scratch, remaining, frames, frequency, chosen, terms, pricing);
    return status;
}
#else
static bool chart_claim(ChartWorker *worker, size_t depth)
{
    (void)worker;
    (void)depth;
    return true;
}

static void chart_closed(ChartWorker *worker, size_t depth)
{
    (void)worker;
    (void)depth;
}
#endif

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

static bool chart_progress_poll(void *context)
{
    WeightedCover *w = context;
    if (w->caller_limits.poll && w->caller_limits.poll(w->caller_limits.context))
        return true;
    if (w->budget.stats.work < w->next_probe ||
        w->next_probe >= w->budget.limits.work) return false;
    if (w->progress == w->checkpoint_progress) {
        nsat_stop(&w->budget, NSAT_UNKNOWN, NSAT_STOP_WORK);
        return true;
    }
    w->checkpoint_progress = w->progress;
    ++w->extensions;
    uint64_t left = w->budget.limits.work - w->next_probe;
    w->next_probe += left < w->probe_work ? left : w->probe_work;
    return false;
}

static IpcChartCoverResult chart_cover(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots, uint64_t probe_work, unsigned threads)
{
    IpcChartCoverResult result = {.status = NSAT_ERROR};
    if (rows <= 0 || columns <= 0 || !masks || !solution) return result;
    if (probe_work && (!limits.work || probe_work > limits.work)) return result;
    int incumbent = 0;
    for (int c = 0; c < columns; ++c) {
        if (solution[c] != 0 && solution[c] != 1) return result;
        incumbent += solution[c];
    }
    if (!incumbent) return result;
    result.terms = incumbent;
    NsatCoverInput input = {.positives = (size_t)rows};
    WeightedCover w = {.budget = {.limits = limits}, .input = &input,
        .mask_words = ((size_t)rows - 1) / 64 + 1,
        .result = {.maximal = (size_t)columns}, .caller_limits = limits,
        .probe_work = probe_work, .next_probe = probe_work};
    if (probe_work && probe_work < limits.work) {
        w.budget.limits.poll = chart_progress_poll;
        w.budget.limits.context = &w;
    }
    Repricing pricing = {.min_branches = min_branches, .pivots = pivots};
    size_t levels = (size_t)incumbent;
    size_t *order = NULL, *scratch = NULL, *frequency = NULL, *chosen = NULL;
    uint64_t *remaining = NULL, *weights = NULL;
    SearchFrame *frames = NULL;
    result.status = NSAT_UNKNOWN;
    w.columns = allocate_matrix(&w.budget, columns, w.mask_words, sizeof(*w.columns));
    remaining = allocate_matrix(&w.budget, levels, w.mask_words, sizeof(*remaining));
    if (!w.columns || !remaining) goto done;
    for (size_t c = 0; c < (size_t)columns; ++c) {
        if (!nsat_work(&w.budget, w.mask_words)) goto done;
        for (size_t word = 0; word < w.mask_words; ++word)
            w.columns[c * w.mask_words + word] = masks[c * w.mask_words + word];
        if (rows % 64)
            w.columns[(c + 1) * w.mask_words - 1] &= (UINT64_C(1) << (rows % 64)) - 1;
    }
    all_rows(&w, remaining);
    for (size_t c = 0; c < (size_t)columns; ++c) {
        if (!nsat_work(&w.budget, w.mask_words)) goto done;
        if (!solution[c]) continue;
        for (size_t word = 0; word < w.mask_words; ++word)
            remaining[word] &= ~column(&w, c)[word];
    }
    for (size_t word = 0; word < w.mask_words; ++word) {
        if (remaining[word]) {
            nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_INPUT);
            goto done;
        }
    }
    if (incumbent == 1) {
        result.status = NSAT_UNSAT;
        goto done;
    }

    weights = nsat_alloc(&w.budget, w.mask_words * 64, sizeof(*weights));
    w.weight_bytes = nsat_alloc(&w.budget, w.mask_words * 8, sizeof(*w.weight_bytes));
    if (!weights || !w.weight_bytes ||
        !nsat_copy(&w.budget, weights, NULL, w.mask_words * 64 * sizeof(*weights)) ||
        !price_weights(&w, weights, pivots ? &pricing.root : NULL)) goto done;
    for (size_t byte = 0; byte < w.mask_words * 8; ++byte) {
        if (!nsat_work(&w.budget, 256)) goto done;
        w.weight_bytes[byte][0] = 0;
        for (unsigned bits = 1; bits < 256; ++bits)
            w.weight_bytes[byte][bits] = w.weight_bytes[byte][bits & (bits - 1)] +
                weights[byte * 8 + (unsigned)__builtin_ctz(bits)];
    }
    order = nsat_alloc(&w.budget, columns, sizeof(*order));
    scratch = allocate_matrix(&w.budget, levels, columns, sizeof(*scratch));
    frequency = nsat_alloc(&w.budget, rows, sizeof(*frequency));
    chosen = nsat_alloc(&w.budget, levels, sizeof(*chosen));
    frames = nsat_alloc(&w.budget, levels, sizeof(*frames));
    if (!order || !scratch || !frequency || !chosen || !frames ||
        !order_columns(&w, order)) goto done;
    if (pivots && incumbent > 2 && !allocate_repricing(&w, &pricing, levels - 2)) goto done;

    while (incumbent > 1) {
        int terms = 0;
        NsatStatus status;
#if defined(HAVE_PTHREAD)
        status = parallel_chart_search(&w, incumbent - 1, order, scratch, remaining,
            frames, frequency, chosen, &terms, &pricing, threads,
            &result.threads, &result.tasks);
#else
        (void)threads;
        status = search(&w, incumbent - 1, order, scratch, remaining,
            frames, frequency, chosen, &terms, &pricing);
#endif
        if (status == NSAT_UNSAT) {
            result.status = status;
            break;
        }
        if (status != NSAT_SAT) goto done;
        if (terms < 1 || terms >= incumbent) {
            nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_PROOF);
            goto done;
        }
        /* Check the witness against every original row before publishing it. */
        all_rows(&w, remaining);
        for (int i = 0; i < terms; ++i) {
            if (!nsat_work(&w.budget, w.mask_words)) goto done;
            if (chosen[i] >= (size_t)columns) {
                nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_PROOF);
                goto done;
            }
            for (size_t word = 0; word < w.mask_words; ++word)
                remaining[word] &= ~column(&w, chosen[i])[word];
        }
        for (size_t word = 0; word < w.mask_words; ++word) {
            if (remaining[word]) {
                nsat_stop(&w.budget, NSAT_ERROR, NSAT_STOP_PROOF);
                goto done;
            }
        }
        memset(solution, 0, (size_t)columns * sizeof(*solution));
        incumbent = 0;
        for (int i = 0; i < terms; ++i) {
            if (!solution[chosen[i]]) ++incumbent;
            solution[chosen[i]] = 1;
        }
        result.terms = incumbent;
        ++w.progress;
    }
    if (incumbent == 1) result.status = NSAT_UNSAT;
done:
    if (w.budget.status != NSAT_OK) result.status = w.budget.status;
    result.nodes = w.result.nodes;
    result.repriced = w.result.repriced;
    result.extensions = w.extensions;
    nsat_free(&w.budget, frames, levels * sizeof(*frames));
    nsat_free(&w.budget, chosen, levels * sizeof(*chosen));
    nsat_free(&w.budget, frequency, (size_t)rows * sizeof(*frequency));
    nsat_free(&w.budget, scratch, levels * (size_t)columns * sizeof(*scratch));
    nsat_free(&w.budget, order, (size_t)columns * sizeof(*order));
    free_repricing(&w, pricing);
    nsat_free(&w.budget, w.weight_bytes, w.mask_words * 8 * sizeof(*w.weight_bytes));
    nsat_free(&w.budget, weights, w.mask_words * 64 * sizeof(*weights));
    nsat_free(&w.budget, remaining, levels * w.mask_words * sizeof(*remaining));
    nsat_free(&w.budget, w.columns, (size_t)columns * w.mask_words * sizeof(*w.columns));
    result.stats = w.budget.stats;
    nsat_release_budget(&w.budget);
    return result;
}

IpcChartCoverResult ipc_chart_cover(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots)
{
    return chart_cover(rows, columns, masks, solution, limits, min_branches, pivots, 0, 1);
}

IpcChartCoverResult ipc_chart_cover_adaptive(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots, uint64_t probe_work)
{
    return chart_cover(rows, columns, masks, solution, limits, min_branches, pivots, probe_work, 1);
}

IpcChartCoverResult ipc_chart_cover_parallel(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots, uint64_t probe_work, unsigned threads)
{
    return chart_cover(rows, columns, masks, solution, limits,
        min_branches, pivots, probe_work, threads);
}
