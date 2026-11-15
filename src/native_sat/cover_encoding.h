#ifndef CCUBES_COVER_ENCODING_H
#define CCUBES_COVER_ENCODING_H

#include "native_sat.h"

#define NSAT_COVER_ENCODING_VERSION 1
#define NSAT_WORK_VERSION 2

typedef struct {
    const int *on, *off;
    size_t inputs, positives, negatives, target;
    const size_t *row_order;
    size_t order_count;
} NsatCoverInput;

typedef enum {
    NSAT_COVER_BASELINE, NSAT_COVER_GROUPS, NSAT_COVER_ORDERING, NSAT_COVER_COMBINED
} NsatCoverStrengthening;

typedef struct {
    size_t pairs, triples, staircase, anchors;
    uint64_t prep_work;
    bool discovery_capped;
} NsatCoverDomainStats;

/* Build the independently checked query and solve it. The caller destroys
 * both handles, including partially initialized handles on failure. */
NsatStatus nsat_cover_solve(const NsatCoverInput *input, NsatCoverStrengthening mode,
                            NsatLimits limits, NsatSolver **solver, NsatVerifier **checker,
                            NsatStopReason *stop, NsatCoverDomainStats *stats);

/* Load solver originals from the care rows. Compare each emitted clause
 * against the independent checker reconstruction in emission order. */
NsatStatus nsat_cover_load_solver(NsatSolver **out, const NsatCoverInput *input,
                                  const NsatVerifier *checker, NsatLimits limits,
                                  NsatStopReason *stop);
int nsat_cover_selector(int inputs, int term, int input, int value);

#endif
