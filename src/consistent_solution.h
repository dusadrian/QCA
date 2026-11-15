#ifndef QCA_CONSISTENT_SOLUTION_H
#define QCA_CONSISTENT_SOLUTION_H

#include <stddef.h>
#include "qca_r.h"

typedef struct {
    const double *data;
    const int *implicants, *depths, *indices, *fuzzy;
    int conditions, cases;
    double consistency, coverage, outcome_sum;
    double *memberships;
} QCAConsistentFit;

typedef struct {
    double *suffix_memberships, *suffix_overlap;
    size_t overlap_stride;
    double slack;
} QCAConsistentBounds;

void qca_cache_pi_memberships(QCAConsistentFit *fit, unsigned int count);
Rboolean consistent_solution(const QCAConsistentFit *fit, int terms, const int selected[]);
void qca_prepare_consistent_bounds(const QCAConsistentFit *fit, unsigned int count,
                                  int max_terms, QCAConsistentBounds *bounds);
Rboolean qca_consistent_prefix_possible(const QCAConsistentFit *fit,
                                        const QCAConsistentBounds *bounds,
                                        const int selected[], int depth, int slots);

#endif
