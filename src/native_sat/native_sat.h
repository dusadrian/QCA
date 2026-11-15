#ifndef CCUBES_NATIVE_SAT_H
#define CCUBES_NATIVE_SAT_H

#include "native_sat_verify.h"

typedef struct NsatSolver NsatSolver;

typedef struct {
    uint64_t restart_base; /* Zero selects the default Luby scale. */
    bool disable_minimisation; /* For differential tests and comparisons. */
    uint64_t reduction_interval; /* Zero selects the default interval. */
    bool disable_reduction;
    /* Optional consecutive variable IDs to branch on before other IDs.
     * They start with a positive decision phase. */
    int priority_first, priority_count;
} NsatSearchOptions;

typedef struct {
    uint64_t decisions, conflicts, binary_propagations, learned;
    uint64_t minimised_literals, restarts, nonchronological_backjumps;
    uint64_t reductions, deleted, compactions, words_reclaimed, moved;
} NsatSearchStats;

/* Initial clause loading and root propagation. A root fixed point carries
 * no SAT conclusion. The solver borrows a sealed checker for propagation.
 * Budget the solver and checker separately at this stage. */
NsatStatus nsat_init(NsatSolver **out, int variables, NsatLimits limits);
NsatStatus nsat_add_original(NsatSolver *s, const int *lits, size_t count);
/* Submit a checked root-level strengthening before propagation or search. */
NsatStatus nsat_add_root_rup(NsatSolver *s, NsatVerifier *v,
                             const int *lits, size_t count,
                             const NsatClauseId *hints, size_t hint_count,
                             NsatClauseId *id);
NsatStatus nsat_propagate_root(NsatSolver *s, NsatVerifier *v);
/* Start from loaded originals. Both terminal results are checked against
 * the verifier's own formula. Model values are available after SAT. */
NsatStatus nsat_solve(NsatSolver *s, NsatVerifier *v, NsatSearchOptions options);
/* Reopen a checked SAT cover query at a strictly smaller positive target.
 * The previous model expires. This is the only supported resumption. */
NsatStatus nsat_solve_tighter_cover(NsatSolver *s, NsatVerifier *v, size_t target);
int nsat_model_value(const NsatSolver *s, int variable);
NsatSearchStats nsat_search_stats(const NsatSolver *s);
int nsat_root_value(const NsatSolver *s, int variable);
NsatClauseId nsat_root_unit_id(const NsatSolver *s, int variable);
NsatStats nsat_stats(const NsatSolver *s);
void nsat_destroy(NsatSolver *s);

#endif
