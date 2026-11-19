#include "qca_consistency_verify.h"
#include "native_sat/native_sat_support.h"
#include <float.h>
#include <math.h>

/* This intentionally small checker uses exhaustive finite-domain reasoning.
 * There are no floating-point pruning bounds: every rejected model is scored
 * in case order using the package's double arithmetic and tolerance rules.
 * Exhausting a resource allowance is UNKNOWN, never an optimality proof. */
typedef struct {
    const QcaConsistencyQuery *q;
    NsatBudget budget;
    int *levels, *term, *coordinates, *primes;
    double *memberships;
    size_t count, capacity, cached_primes, largest_singleton_group;
} ConsistencyCheck;

typedef struct {
    uint64_t known, covered;
} CoverageWord;

typedef struct {
    CoverageWord *masks;
    size_t words, bytes;
} CoverageCache;

static bool invalid(ConsistencyCheck *v)
{
    return nsat_stop(&v->budget, NSAT_ERROR, NSAT_STOP_INPUT);
}

static bool validate(ConsistencyCheck *v)
{
    const QcaConsistencyQuery *q = v->q;
    size_t cells;
    if (!q || !q->tt || !q->data || !q->fuzzy || q->conditions < 1 || q->conditions == INT_MAX ||
        q->rows < 1 || q->cases < 1 || q->depth < 1 || q->depth > q->conditions ||
        !isfinite(q->pi_cons) || q->pi_cons < 0 || q->pi_cons > 1 ||
        !isfinite(q->sol_cons) || q->sol_cons < 0 || q->sol_cons > 1 ||
        !isfinite(q->sol_cov) || q->sol_cov < 0 || q->sol_cov > 1 ||
        !nsat_mul((size_t)q->conditions + 1, (size_t)q->rows, &cells) ||
        cells > SIZE_MAX / sizeof(int) ||
        !nsat_mul((size_t)q->conditions + 1, (size_t)q->cases, &cells) ||
        cells > SIZE_MAX / sizeof(double)) return invalid(v);
    size_t width = (size_t)q->conditions;
    v->levels = nsat_alloc(&v->budget, width, sizeof(int));
    v->term = nsat_alloc(&v->budget, width, sizeof(int));
    v->coordinates = nsat_alloc(&v->budget, width, sizeof(int));
    if (!v->levels || !v->term || !v->coordinates) return false;
    if (!nsat_copy(&v->budget, v->term, NULL, width * sizeof(int))) return false;
    bool positive = false;
    for (int c = 0; c <= q->conditions; ++c) {
        if (c < q->conditions) {
            if (q->fuzzy[c] != 0 && q->fuzzy[c] != 1) return invalid(v);
            v->levels[c] = 1;
        }
        for (int r = 0; r < q->rows; ++r) {
            if (!nsat_work(&v->budget, 1)) return false;
            int value = q->tt[(size_t)c * q->rows + r];
            if (value < 0 || value == INT_MAX) return invalid(v);
            if (c == q->conditions) {
                if (value > 1) return invalid(v);
                positive |= value == 1;
            } else {
                if (q->fuzzy[c] && value > 1) return invalid(v);
                if (value + 1 > v->levels[c]) v->levels[c] = value + 1;
            }
        }
        for (int r = 0; r < q->cases; ++r) {
            if (!nsat_work(&v->budget, 1)) return false;
            double value = q->data[(size_t)c * q->cases + r];
            if (!isfinite(value) || value < 0) return invalid(v);
            if (c == q->conditions || q->fuzzy[c]) {
                if (value > 1) return invalid(v);
            } else if (value != floor(value)) return invalid(v);
        }
    }
    return positive || invalid(v);
}

/* Cache the end of a wildcard run in the verifier's private PI copy. */
static int next_fixed(ConsistencyCheck *v, int *term, int start) {
    if (term[start] < 0) return -term[start];
    int c = start + 1;
    while (c < v->q->conditions && term[c] == 0) {
        if (!nsat_work(&v->budget, 1)) return v->q->conditions;
        ++c;
    }
    if (c < v->q->conditions && term[c] < 0) c = -term[c];
    term[start] = -c;
    return c;
}

static bool matches(ConsistencyCheck *v, int *term, int row) {
    for (int c = 0; c < v->q->conditions; ) {
        if (!nsat_work(&v->budget, 1)) return false;
        if (term[c] <= 0) {
            c = next_fixed(v, term, c);
        } else {
            if (term[c] - 1 != v->q->tt[(size_t)c * v->q->rows + row]) return false;
            ++c;
        }
    }
    return v->budget.status == NSAT_OK;
}

static double membership(ConsistencyCheck *v, int *term, int row) {
    double result = 1.0;
    for (int c = 0; c < v->q->conditions; ) {
        if (!nsat_work(&v->budget, 1)) return 0;
        if (term[c] <= 0) {
            c = next_fixed(v, term, c);
            continue;
        }
        double value = v->q->data[(size_t)c * v->q->cases + row];
        if (v->q->fuzzy[c]) {
            if (term[c] == 1) value = 1.0 - value;
        } else {
            value = value == term[c] - 1 ? 1.0 : 0.0;
        }
        if (value < result) result = value;
        ++c;
    }
    return result;
}

static bool candidate_matches(ConsistencyCheck *v, int row, int depth) {
    const QcaConsistencyQuery *q = v->q;
    for (int i = 0; i < depth; ++i) {
        if (!nsat_work(&v->budget, 1)) return false;
        int c = v->coordinates[i];
        if (v->term[c] - 1 != q->tt[(size_t)c * q->rows + row]) return false;
    }
    return true;
}

static double candidate_membership(ConsistencyCheck *v, int row, int depth) {
    double result = 1.0;
    /* The enumerator lists fixed conditions in increasing column order. */
    for (int i = 0; i < depth; ++i) {
        if (!nsat_work(&v->budget, 1)) return 0;
        int c = v->coordinates[i];
        int literal = v->term[c];
        double value = v->q->data[(size_t)c * v->q->cases + row];
        if (v->q->fuzzy[c]) {
            if (literal == 1) value = 1.0 - value;
        } else {
            value = value == literal - 1 ? 1.0 : 0.0;
        }
        if (value < result) result = value;
    }
    return result;
}

static bool eligible(ConsistencyCheck *v, int depth) {
    const QcaConsistencyQuery *q = v->q;
    bool covers = false;
    for (int r = 0; r < q->rows; ++r) {
        if (!nsat_work(&v->budget, 1)) return false;
        bool on = q->tt[(size_t)q->conditions * q->rows + r] != 0;
        if ((on && !covers) || (!on && q->pi_cons == 0)) {
            bool match = candidate_matches(v, r, depth);
            if (on) covers |= match;
            else if (match) return false;
        }
    }
    if (!covers || q->pi_cons == 0) return covers;
    double total = 0.0, overlap = 0.0;
    const double *outcome = q->data + (size_t)q->conditions * q->cases;
    for (int r = 0; r < q->cases; ++r) {
        double value = candidate_membership(v, r, depth);
        if (v->budget.status != NSAT_OK) return false;
        total += value;
        overlap += value < outcome[r] ? value : outcome[r];
    }
    /* consistency() rejects NaN; never turn unsupported arithmetic into a
     * certificate. The independent comparison matches !altb(score, pi.cons). */
    if (!(total > 0)) return invalid(v);
    double score = overlap / total;
    double scale = fabs(score) < fabs(q->pi_cons) ? fabs(q->pi_cons) : fabs(score);
    return !(q->pi_cons - score > scale * DBL_EPSILON);
}

static bool add_prime(ConsistencyCheck *v) {
    size_t width = (size_t)v->q->conditions;
    /* Checking every retained subterm matters: PI consistency is not
     * monotone, so testing just the immediate parents would be unsound. */
    for (size_t p = 0; p < v->count; ++p) {
        int *prime = v->primes + p * width;
        /* Keep the first coordinate available for duplicate marking in find_prime(). */
        if (!nsat_work(&v->budget, 1)) return false;
        if (prime[0] && prime[0] != v->term[0]) continue;
        int c = 1;
        while (c < v->q->conditions) {
            if (!nsat_work(&v->budget, 1)) return false;
            if (prime[c] <= 0) {
                c = next_fixed(v, prime, c);
            } else {
                if (prime[c] != v->term[c]) break;
                ++c;
            }
        }
        if (c == v->q->conditions) return v->budget.status == NSAT_OK;
    }
    if (v->count == v->capacity) {
        size_t capacity = v->capacity ? v->capacity * 2 : 16, cells;
        if (capacity < v->capacity || capacity > INT_MAX ||
            !nsat_mul(capacity, width, &cells)) return invalid(v);
        int *primes = nsat_alloc(&v->budget, cells, sizeof(int));
        if (!primes) return false;
        if (!nsat_copy(&v->budget, primes, v->primes, v->count * width * sizeof(int))) {
            nsat_free(&v->budget, primes, cells * sizeof(int));
            return false;
        }
        nsat_free(&v->budget, v->primes, v->capacity * width * sizeof(int));
        v->primes = primes;
        v->capacity = capacity;
    }
    if (!nsat_copy(&v->budget, v->primes + v->count * width, v->term, width * sizeof(int)))
        return false;
    ++v->count;
    return true;
}

static bool rebuild_primes(ConsistencyCheck *v, bool index_singletons) {
    if (!validate(v)) return false;
    index_singletons = index_singletons && v->q->depth == 1;
    const int n = v->q->conditions;
    for (int depth = 1; depth <= v->q->depth; ++depth) {
        if (!nsat_work(&v->budget, (uint64_t)depth)) return false;
        for (int c = 0; c < depth; ++c) v->coordinates[c] = c;
        bool more = true;
        while (more) {
            if (!nsat_work(&v->budget, (uint64_t)depth)) return false;
            for (int c = 0; c < depth; ++c) v->term[v->coordinates[c]] = 1;
            for (;;) {
                if (!nsat_work(&v->budget, 1)) return false;
                if (eligible(v, depth) && !add_prime(v)) return false;
                if (v->budget.status != NSAT_OK) return false;
                int i = depth - 1;
                while (i >= 0) {
                    if (!nsat_work(&v->budget, 1)) return false;
                    int c = v->coordinates[i];
                    if (v->term[c] < v->levels[c]) { ++v->term[c]; break; }
                    v->term[c] = 1;
                    --i;
                }
                if (i < 0) break;
            }
            for (int c = 0; c < depth; ++c) {
                if (!nsat_work(&v->budget, 1)) return false;
                v->term[v->coordinates[c]] = 0;
            }
            if (index_singletons) {
                int c = v->coordinates[0];
                size_t start = c ? (size_t)v->levels[c - 1] : 0;
                size_t size = v->count - start;
                if (size > v->largest_singleton_group) v->largest_singleton_group = size;
                /* At depth one this condition's level count has had its last use. */
                v->levels[c] = (int)v->count;
            }
            int i = depth - 1;
            while (i >= 0 && v->coordinates[i] == n - depth + i) --i;
            more = i >= 0;
            if (more) {
                ++v->coordinates[i];
                while (++i < depth) v->coordinates[i] = v->coordinates[i - 1] + 1;
            }
        }
    }
    return true;
}

/* Reconstructed PIs are ordered by depth, then by the positions of their
 * fixed conditions, then by their values. Dropping subsumed PIs preserves
 * that order. Search it directly, using only the private PI array. */
static int find_chart_prime(ConsistencyCheck *v, const int *term) {
    int width = v->q->conditions;
    size_t lo = 0, hi = v->count;
    while (lo < hi) {
        size_t mid = lo + (hi - lo) / 2;
        int *prime = v->primes + mid * (size_t)width;
        int depth = 0, positions = 0, values = 0;
        for (int c = 0; c < width; ++c) {
            if (!nsat_work(&v->budget, 1)) return -1;
            int a = term[c], b = prime[c];
            if (a < 0) return -1;
            if (b < 0) b = c == 0 ? -1 - b : 0;
            depth += (a != 0) - (b != 0);
            if (!positions && (a == 0) != (b == 0)) positions = a ? -1 : 1;
            if (!values && a != b) values = a < b ? -1 : 1;
        }
        int order = depth ? depth : (positions ? positions : values);
        if (order < 0) hi = mid;
        else if (order > 0) lo = mid + 1;
        else {
            if (prime[0] < 0) return -1;
            prime[0] = -1 - prime[0];
            return (int)mid;
        }
    }
    return -1;
}

static int find_prime(ConsistencyCheck *v, const int *term, bool mark) {
    size_t width = (size_t)v->q->conditions;
    for (size_t p = 0; p < v->count; ++p) {
        int *prime = v->primes + p * width;
        size_t c = 0;
        /* During the chart check the rebuilt PIs are private and used only
         * for lookup. Mark a match in its first coordinate, including zero. */
        if (mark) {
            if (!nsat_work(&v->budget, 1)) return -1;
            int first = prime[0] < 0 ? -1 - prime[0] : prime[0];
            if (term[0] != first) continue;
            c = 1;
        }
        for (; c < width; ++c) {
            if (!nsat_work(&v->budget, 1)) return -1;
            int value = prime[c] < 0 ? 0 : prime[c];
            if (term[c] != value) break;
        }
        if (c == width) {
            if (mark) {
                if (prime[0] < 0) return -1;
                prime[0] = -1 - prime[0];
            }
            return (int)p;
        }
    }
    return -1;
}

/* The initial scan pairs each witness coordinate with at most one PI
 * coordinate. Every hit costs the ordinary scan at least width plus its
 * rank; every miss costs at least count. The caller bounds width plus
 * the largest group's remaining comparisons by count. */
static int find_singleton_prime(ConsistencyCheck *v, const int *term) {
    size_t width = (size_t)v->q->conditions;
    size_t first = 0, end = 0, coordinate = width;
    bool first_matches = false;
    for (size_t c = 0; c < width; ++c) {
        if (!nsat_work(&v->budget, 1)) return -1;
        if (!term[c]) continue;
        if (term[c] < 0 || coordinate != width) return -1;
        coordinate = c;
        first = c ? (size_t)v->levels[c - 1] : 0;
        end = (size_t)v->levels[c];
        first_matches = first < end && v->primes[first * width + c] == term[c];
    }
    if (first == end) return -1;
    if (first_matches) return (int)first;
    int value = term[coordinate];
    for (size_t p = first + 1; p < end; ++p) {
        if (!nsat_work(&v->budget, 1)) return -1;
        int found = v->primes[p * width + coordinate];
        if (found == value) return (int)p;
        if (found > value) break;
    }
    return -1;
}

static void *allocate_case_rows(ConsistencyCheck *v, size_t count, size_t row_bytes) {
    /* This is optional and there are no further mandatory allocations.
     * Use a trial budget so memory failure cannot stop the uncached proof.
     * The uninitialized allocation takes no work and adds no poll point. */
    NsatBudget optional = v->budget;
    NsatResources resources;
    if (optional.limits.shared) {
        resources = *optional.limits.shared;
        optional.limits.shared = &resources;
    }
    optional.limits.poll = NULL;
    void *rows = nsat_alloc(&optional, count, row_bytes);
    if (!rows) return NULL;
    v->budget.stats.memory_bytes = optional.stats.memory_bytes;
    v->budget.stats.peak_memory_bytes = optional.stats.peak_memory_bytes;
    if (v->budget.limits.shared) {
        v->budget.limits.shared->memory = resources.memory;
        v->budget.limits.shared->peak_memory = resources.peak_memory;
    }
    return rows;
}

/* Single-prime models are scored in reconstructed PI order. Store the first
 * memberships as they are computed there. A later subset can use the cache
 * only after every singleton has been rejected, so every cached row is ready. */
static void prepare_memberships(ConsistencyCheck *v, int terms) {
    const QcaConsistencyQuery *q = v->q;
    if (terms <= 2 || q->sol_cons == 0 || q->conditions < 2 || FLT_EVAL_METHOD != 0) return;

    size_t row_bytes;
    if (!nsat_mul((size_t)q->cases, sizeof(double), &row_bytes)) return;
    size_t available = nsat_memory_available(&v->budget);
    const size_t max_cache_bytes = 8u * 1024u * 1024u;
    if (available > max_cache_bytes) available = max_cache_bytes;
    size_t count = available / row_bytes;
    if (count > v->count) count = v->count;
    if (!count) return;

    v->memberships = allocate_case_rows(v, count, row_bytes);
    if (v->memberships) v->cached_primes = count;
}

static CoverageCache prepare_coverage(ConsistencyCheck *v) {
    size_t words = ((size_t)v->q->rows + 63) / 64;
    size_t slots, bytes;
    if (!nsat_mul(v->count, words + 1, &slots) || slots > SIZE_MAX - words ||
        !nsat_mul(slots + words, sizeof(CoverageWord), &bytes) ||
        bytes > 8u * 1024u * 1024u) return (CoverageCache){0};
    return (CoverageCache){.masks = allocate_case_rows(v, 1, bytes),
                           .words = words, .bytes = bytes};
}

/* The witness fills the positive masks. A prime's header records the number
 * of initialized words following it, extended only when those words are used. */
static CoverageWord *prime_coverage(const CoverageCache *cache, size_t prime) {
    return cache->masks + cache->words + prime * (cache->words + 1);
}

static bool feasible_with_coverage(ConsistencyCheck *v, const int *selected, int terms,
                                   const CoverageCache *cache) {
    /* Singletons visit every prime before any larger subset uses its masks. */
    if (terms == 1) prime_coverage(cache, (size_t)selected[0])->known = 0;
    for (size_t w = 0; w < cache->words; ++w) {
        if (!nsat_work(&v->budget, 1)) return false;
        /* A word lookup replaces the row visits in this block of the table. */
        uint64_t pending = cache->masks[w].covered;
        while (pending) {
            int low = __builtin_ctzll(pending);
            int row = (int)(64 * w + (size_t)low);
            uint64_t bit = UINT64_C(1) << low;
            for (int t = 0; t < terms; ++t) {
                size_t prime = (size_t)selected[t];
                CoverageWord *header = prime_coverage(cache, prime);
                CoverageWord *word = w <= header->known ? header + w + 1 : NULL;
                if (w == header->known) {
                    *word = (CoverageWord){0};
                    ++header->known;
                }
                if (word && (word->known & bit)) {
                    if (!nsat_work(&v->budget, 1)) return false;
                } else {
                    bool hit = matches(v, v->primes + prime * v->q->conditions, row);
                    if (v->budget.status != NSAT_OK) return false;
                    if (word) {
                        word->known |= bit;
                        if (hit) word->covered |= bit;
                    } else if (hit) pending &= ~bit;
                }
                if (word) pending &= ~word->covered;
                if (!(pending & bit)) break;
            }
            if (pending & bit) return false;
        }
    }
    return true;
}

static bool feasible(ConsistencyCheck *v, const int *selected, int terms,
                     CoverageWord *positive_masks) {
    const QcaConsistencyQuery *q = v->q;
    if (q->sol_cons == 0) {
        uint64_t positive = 0;
        for (int r = 0; r < q->rows; ++r) {
            if (!nsat_work(&v->budget, 1)) return false;
            bool on = q->tt[(size_t)q->conditions * q->rows + r] != 0;
            if (positive_masks) {
                unsigned bit = (unsigned)r % 64;
                if (on) positive |= UINT64_C(1) << bit;
                if (bit == 63 || r == q->rows - 1) {
                    positive_masks[(size_t)r / 64].covered = positive;
                    positive = 0;
                }
            }
            if (!on) continue;
            bool covered = false;
            for (int t = 0; t < terms && !covered; ++t)
                covered = matches(v, v->primes + (size_t)selected[t] * q->conditions, r);
            if (!covered || v->budget.status != NSAT_OK) return false;
        }
        return true;
    }
    double total = 0.0, overlap = 0.0, outcome_sum = 0.0;
    const double *outcome = q->data + (size_t)q->conditions * q->cases;
    for (int r = 0; r < q->cases; ++r) {
        double value = 0.0;
        for (int t = 0; t < terms; ++t) {
            size_t prime = (size_t)selected[t];
            double term;
            if (prime < v->cached_primes && terms > 1) {
                if (!nsat_work(&v->budget, 1)) return false;
                term = v->memberships[prime * (size_t)q->cases + (size_t)r];
            } else {
                term = membership(v, v->primes + prime * (size_t)q->conditions, r);
                if (v->budget.status != NSAT_OK) return false;
                if (prime < v->cached_primes)
                    v->memberships[prime * (size_t)q->cases + (size_t)r] = term;
            }
            if (term > value) value = term;
        }
        total += value;
        overlap += value < outcome[r] ? value : outcome[r];
        outcome_sum += outcome[r];
    }
    if (!(total > 0) || !(outcome_sum > 0)) return false;
    double consistency = overlap / total, coverage = overlap / outcome_sum;
    return (consistency > q->sol_cons || fabs(consistency - q->sol_cons) <= DBL_EPSILON) &&
        (coverage > q->sol_cov || fabs(coverage - q->sol_cov) <= DBL_EPSILON);
}

/* The previous subset's common prefix is still valid. Keep the original
 * case order for all sums. Each lookup replaces one or more cached lookups
 * from the original scorer, and suffix maxima are stored during that pass. */
static bool feasible_with_maxima(ConsistencyCheck *v, const int *model, int terms,
                                 int unchanged, double *maxima) {
    const QcaConsistencyQuery *q = v->q;
    size_t cases = (size_t)q->cases;
    int prefix = unchanged > 1 ? unchanged : 1;
    const double *previous = prefix == 1
        ? v->memberships + (size_t)model[0] * cases
        : maxima + (size_t)(prefix - 2) * cases;
    const double *last = v->memberships + (size_t)model[terms - 1] * cases;
    const double *outcome = q->data + (size_t)q->conditions * cases;
    double total = 0.0, overlap = 0.0, outcome_sum = 0.0;
    for (size_t r = 0; r < cases; ++r) {
        if (!nsat_work(&v->budget, 1)) return false;
        double value = previous[r] > 0.0 ? previous[r] : 0.0;
        for (int t = prefix; t < terms - 1; ++t) {
            if (!nsat_work(&v->budget, 1)) return false;
            double term = v->memberships[(size_t)model[t] * cases + r];
            if (term > value) value = term;
            maxima[(size_t)(t - 1) * cases + r] = value;
        }
        if (!nsat_work(&v->budget, 1)) return false;
        if (last[r] > value) value = last[r];
        total += value;
        overlap += value < outcome[r] ? value : outcome[r];
        outcome_sum += outcome[r];
    }
    if (!(total > 0) || !(outcome_sum > 0)) return false;
    double consistency = overlap / total, coverage = overlap / outcome_sum;
    return (consistency > q->sol_cons || fabs(consistency - q->sol_cons) <= DBL_EPSILON) &&
        (coverage > q->sol_cov || fabs(coverage - q->sol_cov) <= DBL_EPSILON);
}

static NsatStatus finish(ConsistencyCheck *v, NsatStatus result)
{
    if (v->budget.status != NSAT_OK) result = v->budget.status;
    size_t width = v->q && v->q->conditions > 0 ? (size_t)v->q->conditions : 0;
    nsat_free(&v->budget, v->levels, width * sizeof(int));
    nsat_free(&v->budget, v->term, width * sizeof(int));
    nsat_free(&v->budget, v->coordinates, width * sizeof(int));
    nsat_free(&v->budget, v->primes, v->capacity * width * sizeof(int));
    if (v->memberships)
        nsat_free(&v->budget, v->memberships, v->cached_primes * (size_t)v->q->cases * sizeof(double));
    return result;
}

/* A chart lookup is the verifier's last operation. Use the index only
 * when even its worst-case work for all submitted columns fits. This
 * preserves early duplicate and invalid-column results at tight limits. */
static bool index_fits_chart(const ConsistencyCheck *v, size_t columns) {
    if (v->count < 16 || columns < 2) return false;
    const NsatBudget *b = &v->budget;
    uint64_t remaining = (b->limits.work ? b->limits.work : UINT64_MAX) - b->stats.work;
    const NsatResources *shared = b->limits.shared;
    if (shared) {
        uint64_t available = (shared->work_limit ? shared->work_limit : UINT64_MAX) - shared->work;
        if (available < remaining) remaining = available;
    }
    unsigned comparisons = 0;
    for (size_t n = v->count; n; n /= 2) ++comparisons;
    return columns <= remaining / comparisons / (unsigned)v->q->conditions;
}

NsatStatus qca_verify_consistency_chart(const QcaConsistencyQuery *q,
    const int *chart, size_t columns, NsatLimits limits)
{
    ConsistencyCheck v = {.q = q, .budget = {.limits = limits}};
    if (!rebuild_primes(&v, false)) return finish(&v, NSAT_UNKNOWN);
    if ((columns && !chart) || columns > SIZE_MAX / sizeof(int) / (size_t)q->conditions)
        return finish(&v, NSAT_ERROR);
    /* Count equality plus membership excludes both omissions and duplicates. */
    bool indexed = index_fits_chart(&v, columns);
    for (size_t p = 0; p < columns; ++p) {
        const int *term = chart + p * q->conditions;
        int found = indexed ? find_chart_prime(&v, term) : find_prime(&v, term, true);
        if (found < 0)
            return finish(&v, NSAT_ERROR);
    }
    return finish(&v, columns == v.count ? NSAT_UNSAT : NSAT_SAT);
}

NsatStatus qca_verify_consistency_minimum(const QcaConsistencyQuery *q,
    const int *chart, size_t columns, const int *selected, int terms, NsatLimits limits)
{
    ConsistencyCheck v = {.q = q, .budget = {.limits = limits}};
    if (!rebuild_primes(&v, true)) return finish(&v, NSAT_UNKNOWN);
    if (!chart || !selected || terms < 1 || (size_t)terms > v.count ||
        columns > SIZE_MAX / sizeof(int) / (size_t)q->conditions)
        return finish(&v, NSAT_ERROR);
    int *model = nsat_alloc(&v.budget, (size_t)terms, sizeof(int));
    NsatStatus result = NSAT_ERROR;
    double *maxima = NULL;
    size_t maxima_bytes = 0;
    CoverageCache coverage = {0};
    if (!model) goto done;
    size_t width = (size_t)q->conditions;
    bool indexed = v.largest_singleton_group && v.count >= 16 && v.count >= width &&
        v.largest_singleton_group - 1 <= v.count - width;
    for (int t = 0; t < terms; ++t) {
        if (selected[t] < 0 || (size_t)selected[t] >= columns) goto done;
        const int *term = chart + (size_t)selected[t] * width;
        model[t] = indexed ? find_singleton_prime(&v, term) : find_prime(&v, term, false);
        if (model[t] < 0) goto done;
        for (int i = 0; i < t; ++i) {
            if (!nsat_work(&v.budget, 1)) goto done;
            if (model[i] == model[t]) goto done;
        }
    }
    if (q->sol_cons == 0 && q->rows > 64 && terms > 2) coverage = prepare_coverage(&v);
    if (!feasible(&v, model, terms, coverage.masks)) goto done;
    prepare_memberships(&v, terms);
    for (int k = 1; k < terms; ++k) {
        if (!nsat_work(&v.budget, (uint64_t)k)) goto done;
        if (k == 3 && v.cached_primes == v.count) {
            size_t row_bytes = (size_t)q->cases * sizeof(double);
            size_t rows = (size_t)(terms - 3);
            if (rows <= nsat_memory_available(&v.budget) / row_bytes) {
                maxima = allocate_case_rows(&v, rows, row_bytes);
                if (maxima) maxima_bytes = rows * row_bytes;
            }
        }
        for (int i = 0; i < k; ++i) model[i] = i;
        int unchanged = 0;
        for (;;) {
            bool fits;
            if (coverage.masks) fits = feasible_with_coverage(&v, model, k, &coverage);
            else if (maxima) fits = feasible_with_maxima(&v, model, k, unchanged, maxima);
            else fits = feasible(&v, model, k, NULL);
            if (fits) { result = NSAT_SAT; goto done; }
            if (v.budget.status != NSAT_OK) goto done;
            int i = k - 1;
            while (i >= 0 && (size_t)model[i] == v.count - (size_t)k + (size_t)i) --i;
            if (i < 0) break;
            unchanged = i;
            ++model[i];
            while (++i < k) model[i] = model[i - 1] + 1;
        }
    }
    result = NSAT_UNSAT;
done:
    nsat_free(&v.budget, coverage.masks, coverage.bytes);
    nsat_free(&v.budget, maxima, maxima_bytes);
    nsat_free(&v.budget, model, (size_t)terms * sizeof(int));
    return finish(&v, result);
}
