#ifndef QCA_CERTIFY_H
#define QCA_CERTIFY_H

#include "native_sat/native_sat_verify.h"

typedef struct {
    const int *on, *off; /* QCA's column-major truth-table matrices. */
    int conditions, positives, negatives;
    NsatResources resources;
    int verified_minimum; /* Checked bound with a witness for this table. */
} QcaCertification;

/* For binary or categorical equality terms, NSAT_UNSAT rules out every
 * smaller cover. The caller supplies a feasible incumbent of this size.
 * UNKNOWN and ERROR provide no minimum certificate. */
NsatStatus qca_check_smaller_model(QcaCertification *table, int incumbent);

/* Prove that every legal term covering an ON row contains a chart PI. */
bool qca_check_complete_chart(QcaCertification *table, const int *chart,
                              const int *implicants, int columns);

/* Prove that every binary minimum model uses only chart PIs.
 * minimum must already be a global minimum for this table. */
bool qca_check_minimum_models(QcaCertification *table, const int *implicants,
                              int columns, int minimum);

#endif
