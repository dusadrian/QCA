#ifndef QCA_ENUMERATE_MODELS_H
#define QCA_ENUMERATE_MODELS_H

#include "qca_rinternals.h"
#include <limits.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

/* Enumeration keeps every PI identity, including equal and dominated columns.
   States partition covers by the first available PI covering a chosen row. */
typedef struct {
    int ready, pivot, next;
} QCAEnumFrame;

typedef struct {
    const int *chart;
    int rows, cols, k, words;
    QCAEnumFrame *frames;
    unsigned char *active, *packed;
    uint64_t *masks, *uncovered;
    int *chosen, *output;
    size_t count, capacity;
    int complete;
} QCAEnum;

static int qca_enum_popcount(uint64_t x) {
#if defined(__GNUC__) || defined(__clang__)
    return __builtin_popcountll(x);
#else
    int count = 0;
    while (x) { x &= x - 1; ++count; }
    return count;
#endif
}

static int qca_enum_int_compare(const void *a, const void *b) {
    int x = *(const int *)a, y = *(const int *)b;
    return (x > y) - (x < y);
}

static void qca_enum_emit(QCAEnum *s) {
    size_t limit = (size_t)INT_MAX / s->k;
    if (s->count == limit) error("Too many enumerated models for the solution matrix.");
    if (s->count == s->capacity) {
        size_t capacity = s->capacity ? s->capacity * 2 : 64;
        if (capacity > limit) capacity = limit;
        if (capacity > SIZE_MAX / sizeof(int) / s->k)
            error("Enumerated solution matrix is too large.");
        if (s->output) s->output = R_Realloc(s->output, capacity * s->k, int);
        else s->output = R_Calloc(capacity * s->k, int);
        s->capacity = capacity;
    }
    int *record = s->output + s->count * s->k;
    memcpy(record, s->chosen, (size_t)s->k * sizeof(int));
    qsort(record, (size_t)s->k, sizeof(int), qca_enum_int_compare);
    ++s->count;
}

static int qca_enum_compare(QCAEnum *s, size_t a, size_t b) {
    const int *x = s->output + a * s->k, *y = s->output + b * s->k;
    for (int i = 0; i < s->k; ++i) {
        if (x[i] != y[i]) return x[i] < y[i] ? -1 : 1;
    }
    return 0;
}

static void qca_enum_swap(QCAEnum *s, size_t a, size_t b) {
    int *x = s->output + a * s->k, *y = s->output + b * s->k;
    for (int i = 0; i < s->k; ++i) {
        int temp = x[i]; x[i] = y[i]; y[i] = temp;
    }
}

static void qca_enum_sift(QCAEnum *s, size_t root, size_t end) {
    while (root < end / 2) {
        size_t child = root * 2 + 1;
        if (child + 1 < end && qca_enum_compare(s, child, child + 1) < 0) ++child;
        if (qca_enum_compare(s, root, child) >= 0) break;
        qca_enum_swap(s, root, child);
        root = child;
    }
}

/* In-place heapsort preserves the previous combination order without global
   comparator state or another full copy of a potentially large result. */
static void qca_enum_sort(QCAEnum *s) {
    for (size_t start = s->count / 2; start > 0;) {
        --start;
        if ((start & 1023U) == 0) R_CheckUserInterrupt();
        qca_enum_sift(s, start, s->count);
    }
    for (size_t end = s->count; end > 1;) {
        --end;
        if ((end & 1023U) == 0) R_CheckUserInterrupt();
        qca_enum_swap(s, 0, end);
        qca_enum_sift(s, 0, end);
    }
}

/* Returns -2 for infeasibility and -1 when coverage is already complete.
   Bounds use only necessary coverage conditions, never PI dominance. */
static int qca_enum_pivot(QCAEnum *s, unsigned char *active,
                          uint64_t *uncovered, int slots) {
    int available = 0, remaining = 0, max_gain = 0;
    for (int w = 0; w < s->words; ++w) remaining += qca_enum_popcount(uncovered[w]);
    for (int c = 0; c < s->cols; ++c) {
        if (!active[c]) continue;
        ++available;
        int gain = 0;
        for (int w = 0; w < s->words; ++w)
            gain += qca_enum_popcount(uncovered[w] & s->masks[(size_t)c * s->words + w]);
        if (gain > max_gain) max_gain = gain;
    }
    if (available < slots || (int64_t)slots * max_gain < remaining) return -2;
    if (!remaining) return -1;

    int pivot = -1, least = INT_MAX, packing = 0;
    memset(s->packed, 0, (size_t)s->cols);
    for (int r = 0; r < s->rows; ++r) {
        if (!(uncovered[r >> 6] & (UINT64_C(1) << (r & 63)))) continue;
        int support = 0, disjoint = 1;
        for (int c = 0; c < s->cols; ++c) {
            if (active[c] && s->chart[(size_t)c * s->rows + r]) {
                ++support;
                if (s->packed[c]) disjoint = 0;
            }
        }
        if (!support) return -2;
        if (support < least) { least = support; pivot = r; }
        /* Rows with disjoint candidate supports require distinct PIs. */
        if (disjoint) {
            if (++packing > slots) return -2;
            for (int c = 0; c < s->cols; ++c)
                if (active[c] && s->chart[(size_t)c * s->rows + r]) s->packed[c] = 1;
        }
    }
    return pivot;
}

static SEXP qca_enum_run(void *data) {
    QCAEnum *s = (QCAEnum *)data;
    size_t levels = (size_t)s->k + 1;
    if (levels > SIZE_MAX / (size_t)s->cols ||
        levels > SIZE_MAX / sizeof(uint64_t) / (size_t)s->words ||
        (size_t)s->cols > SIZE_MAX / sizeof(uint64_t) / (size_t)s->words)
        error("Model enumeration workspace is too large.");
    s->frames = R_Calloc(levels, QCAEnumFrame);
    s->active = R_Calloc(levels * s->cols, unsigned char);
    s->packed = R_Calloc(s->cols, unsigned char);
    s->masks = R_Calloc((size_t)s->cols * s->words, uint64_t);
    s->uncovered = R_Calloc(levels * s->words, uint64_t);
    s->chosen = R_Calloc(s->k, int);
    for (int c = 0; c < s->cols; ++c) {
        if ((c & 1023) == 0) R_CheckUserInterrupt();
        for (int r = 0; r < s->rows; ++r)
            if (s->chart[(size_t)c * s->rows + r])
                s->masks[(size_t)c * s->words + (r >> 6)] |= UINT64_C(1) << (r & 63);
    }
    for (int r = 0; r < s->rows; ++r)
        s->uncovered[r >> 6] |= UINT64_C(1) << (r & 63);
    memset(s->active, 1, (size_t)s->cols);

    int depth = 0;
    unsigned long long ticks = 0;
    while (depth >= 0) {
        if ((++ticks & 1023ULL) == 0) R_CheckUserInterrupt();
        QCAEnumFrame *frame = s->frames + depth;
        unsigned char *active = s->active + (size_t)depth * s->cols;
        uint64_t *uncovered = s->uncovered + (size_t)depth * s->words;
        if (!frame->ready) {
            if (depth == s->k) {
                int empty = 1;
                for (int w = 0; w < s->words; ++w) if (uncovered[w]) { empty = 0; break; }
                if (empty) qca_enum_emit(s);
                --depth;
                continue;
            }
            frame->pivot = qca_enum_pivot(s, active, uncovered, s->k - depth);
            if (frame->pivot == -2) { --depth; continue; }
            frame->next = 0;
            frame->ready = 1;
        }
        int c = frame->next;
        while (c < s->cols && (!active[c] ||
            (frame->pivot >= 0 && !s->chart[(size_t)c * s->rows + frame->pivot]))) ++c;
        if (c == s->cols) { --depth; continue; }
        frame->next = c + 1;
        /* Later siblings exclude this PI; the child includes it. Hence a
           cover containing multiple pivot candidates is emitted only once. */
        active[c] = 0;
        s->chosen[depth] = c + 1;
        memcpy(active + s->cols, active, (size_t)s->cols);
        uint64_t *child = uncovered + s->words;
        for (int w = 0; w < s->words; ++w)
            child[w] = uncovered[w] & ~s->masks[(size_t)c * s->words + w];
        s->frames[++depth].ready = 0;
    }
    qca_enum_sort(s);
    if (!s->output) s->output = R_Calloc(1, int);
    s->complete = 1;
    return R_NilValue;
}

static void qca_enum_cleanup(void *data) {
    QCAEnum *s = (QCAEnum *)data;
    R_Free(s->frames); R_Free(s->active); R_Free(s->packed);
    R_Free(s->masks); R_Free(s->uncovered); R_Free(s->chosen);
    if (!s->complete) R_Free(s->output);
}

static void qca_enumerate_fixed_models(const int *chart, int rows, unsigned int cols,
                                      int k, int **solutions, int *nr, int *nc) {
    if (k <= 0 || cols > INT_MAX || (unsigned int)k > cols) {
        *nr = 0; *nc = 0;
        return;
    }
    QCAEnum s = {.chart = chart, .rows = rows, .cols = (int)cols, .k = k,
                 .words = rows / 64 + (rows % 64 != 0)};
    R_ExecWithCleanup(qca_enum_run, &s, qca_enum_cleanup, &s);
    R_Free(*solutions);
    *solutions = s.output;
    *nr = k;
    *nc = (int)s.count;
}

#endif
