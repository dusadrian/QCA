/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#ifndef CCUBES_IPC_WEIGHTED_VERIFY_H
#define CCUBES_IPC_WEIGHTED_VERIFY_H

#include "cover_encoding.h"

/* UNSAT independently confirms that lower - 1 cubes cannot cover the table.
 * The weights and proposed per-cube maximum are untrusted proof data. */
NsatStatus ipc_verify_weighted_bound(const NsatCoverInput *input,
                                     const uint64_t *weights, size_t count,
                                     uint64_t maximum, int lower,
                                     NsatLimits limits, NsatStats *stats);


/* Optional hints for a chart of at most 64 ON rows. The table has
 * IPC_WEIGHT_HINT_BUCKETS entries. Nonzero indices are one-based into
 * count consecutive vectors of rows weights each. Hints never certify
 * a bound without checking it against independently derived columns. */
#define IPC_WEIGHT_HINT_BUCKETS 1024u
typedef struct {
    uint64_t residual;
    int slots;
    size_t index;
} IpcWeightHint;
typedef struct {
    const IpcWeightHint *entries;
    const uint64_t *weights;
    size_t count;
    size_t rows;
} IpcWeightHints;
static inline size_t ipc_weight_hint_bucket(uint64_t residual, int slots)
{
    uint64_t key = residual ^ ((uint64_t)slots * UINT64_C(0x9e3779b97f4a7c15));
    key ^= key >> 33;
    key *= UINT64_C(0xff51afd7ed558ccd);
    key ^= key >> 33;
    return key & (IPC_WEIGHT_HINT_BUCKETS - 1);
}

NsatStatus ipc_verify_weighted_hint_bound(const NsatCoverInput *input, const uint64_t *weights,
    size_t count, uint64_t maximum, int lower, NsatLimits limits, NsatStats *stats,
    const IpcWeightHints *hints);

#endif
