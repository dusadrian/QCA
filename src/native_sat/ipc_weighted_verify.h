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

#endif
