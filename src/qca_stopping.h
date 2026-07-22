#ifndef QCA_STOPPING_H
#define QCA_STOPPING_H

#include <stdbool.h>

/*
 * Stopping state for one QCA minimization.  Cardinality certification and
 * model-enumeration completion are deliberately separate: proving that the
 * incumbent uses the fewest terms does not prove that every minimum model has
 * already been generated.
 */
typedef struct {
    bool supported;
    bool diagnostic_checked;
    bool certification_required;
    bool cardinality_certified;
    bool enumeration_complete;
    int agreement_horizon;
    int coverage_horizon;
    int cover_lower_bound;
} QCAStoppingState;

void qca_stopping_state_init(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    const int noflevels[],
    int nconds,
    int posrows,
    int negrows,
    double picons,
    double solcons
);

void qca_stopping_observe_coverage(
    QCAStoppingState *state,
    int level,
    bool covered
);

/*
 * Inspect the first minimum-cardinality plateau.  A delayed private-witness
 * pair is a deterministic warning that a deeper PI may merge incumbent terms.
 * A zero warning leaves the plateau heuristic; it is not an exactness proof.
 */
bool qca_stopping_observe_plateau(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    int nconds,
    int posrows,
    int negrows,
    const int pichart[],
    int foundPI,
    const int selected_indices[],
    int selected_terms
);

/* Prepare the global horizon and incompatibility lower-bound certificate. */
bool qca_stopping_prepare_certificate(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    int nconds,
    int posrows,
    int negrows
);

/* Update and return the cardinality-certificate state. */
bool qca_stopping_cardinality_certified(
    QCAStoppingState *state,
    int level,
    int cover_size,
    bool boundary_exact
);

#endif
