#ifndef CCUBES_NATIVE_SAT_VERIFY_H
#define CCUBES_NATIVE_SAT_VERIFY_H

#include <stdbool.h>
#include <stddef.h>
#include <stdint.h>

/* Literals are signed, one-based variable indices. Zero is never a literal. */
typedef uint64_t NsatClauseId;
typedef enum {
    NSAT_OK, NSAT_UNKNOWN, NSAT_ERROR, NSAT_UNSAT, NSAT_ROOT_FIXED_POINT,
    NSAT_SAT
} NsatStatus;
typedef enum {
    NSAT_STOP_NONE, NSAT_STOP_WORK, NSAT_STOP_MEMORY, NSAT_STOP_CANCELLED,
    NSAT_STOP_INPUT, NSAT_STOP_STATE, NSAT_STOP_PROOF, NSAT_STOP_ALLOCATION,
    NSAT_STOP_CONFLICTS
} NsatStopReason;

/* Optional accounting shared by the encoder, solver and checker of one
 * serial query or series of queries. Initialize the counters to zero. */
typedef struct {
    uint64_t work_limit, work;
    size_t memory_limit, memory, peak_memory;
    NsatStopReason stop;
} NsatResources;

typedef struct {
    /* Zero disables a limit. poll may check the caller's monotonic deadline. */
    uint64_t work;
    size_t memory_bytes;
    bool (*poll)(void *context);
    void *context;
    uint64_t conflicts; /* Solver only. Zero is unlimited. */
    NsatResources *shared;
} NsatLimits;

typedef struct {
    uint64_t work;
    size_t memory_bytes;
    size_t peak_memory_bytes;
    NsatStopReason stop;
    uint64_t lemmas_checked, clauses_deleted, index_rebuilds;
} NsatStats;

typedef struct NsatVerifier NsatVerifier;

/* Flat point tables contain 1/2 for Boolean 0/1. A NULL permutation with
 * length zero means file order. Dimensions and point buffers are borrowed
 * during initialization. The cover constructor returns a sealed verifier.
 * In v1 it requires n, on_count, off_count and target all positive.
 *
 * The standalone CNF constructor is for trusted harnesses. Supply normalized
 * originals (no duplicate or complementary literals), then seal. Each
 * original receives the next ID starting at 1. After either constructor's
 * failure, or any UNKNOWN or ERROR, destroy the instance. */
NsatStatus nsat_verify_init_cnf(NsatVerifier **out, int variables, NsatLimits limits);
NsatStatus nsat_verify_init_cover(NsatVerifier **out, const int *on, const int *off,
                                  size_t n, size_t on_count, size_t off_count,
                                  size_t target, const size_t *row_order,
                                  size_t order_count, NsatLimits limits);
NsatStatus nsat_verify_add_original(NsatVerifier *v, const int *lits, size_t count);
NsatStatus nsat_verify_seal_originals(NsatVerifier *v);
/* Bind a freshly sealed checker at the start of a solve session. */
NsatStatus nsat_verify_begin(NsatVerifier *v);

/* Sealed checker events. All buffers are borrowed for this call only.
 * Original IDs cannot be deleted. Conclusion requires a live empty clause. */
NsatStatus nsat_verify_derive(NsatVerifier *v, NsatClauseId id,
                              const int *lits, size_t count,
                              const NsatClauseId *hints, size_t hint_count);
NsatStatus nsat_verify_delete(NsatVerifier *v, NsatClauseId id);
NsatStatus nsat_verify_conclude_unsat(NsatVerifier *v, NsatClauseId id);
/* One signed literal per variable, in arbitrary order. Check the original
 * CNF and all registered target units. */
NsatStatus nsat_verify_conclude_sat(NsatVerifier *v, const int *model, size_t count);
/* After a checked cover model, permanently lower its active target. The
 * checker computes the single new premise from its own cover metadata. */
NsatStatus nsat_verify_tighten_cover(NsatVerifier *v, size_t target,
                                      NsatClauseId *id, int *literal);
size_t nsat_verify_cover_target(const NsatVerifier *v);
size_t nsat_verify_cover_slots(const NsatVerifier *v);

/* Read-only diagnostics. Clause views last until the next modifying call. */
NsatStats nsat_verify_stats(const NsatVerifier *v);
int nsat_verify_variables(const NsatVerifier *v);
NsatClauseId nsat_verify_original_count(const NsatVerifier *v);
bool nsat_verify_clause(const NsatVerifier *v, NsatClauseId id,
                        const int **lits, size_t *count);
void nsat_verify_destroy(NsatVerifier *v);

#endif
