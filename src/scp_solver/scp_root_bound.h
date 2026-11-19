#ifndef QCA_SCP_ROOT_BOUND_H
#define QCA_SCP_ROOT_BOUND_H

#include <math.h>
#include <stdint.h>
#include <string.h>

/* One PI costs one unit. At this scale even INT_MAX row weights of one
 * fit in int64_t. Only the sum of negative reduced costs can underflow. */
#define QCA_ROOT_SCALE (INT64_C(1) << 32)

typedef struct {
    int64_t value;
    int64_t *cost;
} QcaRootBound;

static inline int qca_verify_root_bound(
    QcaRootBound *proof, int rows, int cols, int **covering,
    const int counts[], const double proposed[]
) {
    proof->value = 0;
    if (!proof->cost) return 0;

    for (int c = 0; c < cols; ++c) proof->cost[c] = QCA_ROOT_SCALE;
    int64_t bound = 0;
    for (int r = 0; r < rows; ++r) {
        double y = proposed[r];
        int64_t weight = 0;
        if (isfinite(y) && y > 0) {
            weight = y >= 1 ? QCA_ROOT_SCALE : (int64_t)(y * QCA_ROOT_SCALE);
        }
        bound += weight;
        for (int j = 0; j < counts[r]; ++j) proof->cost[covering[r][j]] -= weight;
    }

    /* L = sum(row weights) + sum(min(0, reduced cost)).
     * Every feasible cover has integer cardinality at least ceil(L). */
    for (int c = 0; c < cols; ++c) {
        int64_t cost = proof->cost[c];
        if (cost >= 0) continue;
        if (bound < INT64_MIN - cost) {
            memset(proof->cost, 0, (size_t)cols * sizeof(*proof->cost));
            return 0;
        }
        bound += cost;
    }
    proof->value = bound;
    return bound > 0 ? (int)((bound - 1) / QCA_ROOT_SCALE + 1) : 0;
}

/* Selecting a PI adds max(0, cost) to L. Omitting it adds max(0, -cost).
 * Equality with the target must retain that choice. */
static inline int qca_root_excludes(const QcaRootBound *proof, int col, int target) {
    int64_t cost = proof->cost ? proof->cost[col] : 0;
    return proof->value > (int64_t)target * QCA_ROOT_SCALE - (cost > 0 ? cost : 0);
}

static inline int qca_root_requires(const QcaRootBound *proof, int col, int target) {
    int64_t cost = proof->cost ? proof->cost[col] : 0;
    return cost < 0 && proof->value > (int64_t)target * QCA_ROOT_SCALE + cost;
}

#endif
