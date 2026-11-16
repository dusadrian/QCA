#ifndef QCA_CONSISTENCY_VERIFY_H
#define QCA_CONSISTENCY_VERIFY_H

#include "native_sat/native_sat_verify.h"

/* Independent finite-domain verifier. Matrices are column-major: tt has
 * categorical 0-based conditions followed by a 0/1 outcome; data contains
 * calibrated cases followed by their outcome memberships. Terms are rows
 * of 1-based equalities, with zero denoting a wildcard. */
typedef struct {
    const int *tt, *fuzzy;
    const double *data;
    int conditions, rows, cases, depth;
    double pi_cons, sol_cons, sol_cov;
} QcaConsistencyQuery;

/* UNSAT certifies that the supplied chart is exactly the eligible PI family
 * within depth. SAT means an omitted PI; ERROR includes invalid extra terms.
 * This is safe for all-model enumeration as well as minimum-size searches. */
NsatStatus qca_verify_consistency_chart(const QcaConsistencyQuery *query,
    const int *chart, size_t columns, NsatLimits limits);

/* The certificate consists of the proposed model's distinct zero-based
 * indices into chart. The verifier rebuilds the PI domain independently,
 * checks that witness and exhausts every smaller PI model. UNSAT certifies
 * minimum size; SAT supplies a refutation (a smaller model exists); ERROR
 * rejects malformed/infeasible witnesses. UNKNOWN never certifies anything.
 * No production PI-generation, fit, pruning or solver routines are used. */
NsatStatus qca_verify_consistency_minimum(const QcaConsistencyQuery *query,
    const int *chart, size_t columns, const int *selected, int terms, NsatLimits limits);

#endif
