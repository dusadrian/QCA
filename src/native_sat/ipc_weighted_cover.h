/* Copyright (c) 2016–2026, Adrian Dusa. See LICENSE. */
#ifndef CCUBES_IPC_WEIGHTED_COVER_H
#define CCUBES_IPC_WEIGHTED_COVER_H

#include "cover_encoding.h"

typedef struct {
    int lower, terms;
    int *cover; /* Transposed. The caller frees it. */
    bool supported;
    NsatStatus status;
    NsatStopReason stop;
    uint64_t compatible, nodes;
    size_t maximal;
    unsigned pivots;
    uint64_t repriced, reprice_pivots, reprice_pruned;
    int first_terms; /* Zero until a smaller cover has been materialized. */
    uint64_t first_cover_work, first_cover_nodes;
    NsatStatus check_status; /* OK when no higher bound was proposed. */
    NsatStats check;
} IpcWeightedCoverResult;

/* input->target is the size of a validated incumbent. A limited search can
 * still return earlier proved bounds and completed smaller covers. A bound
 * above the supplied lower bound has passed independent verification. */
IpcWeightedCoverResult ipc_weighted_cover(const NsatCoverInput *input, int lower,
                                          NsatLimits limits);

/* Experimental per-node pricing. The ordinary entry point keeps root weights. */
IpcWeightedCoverResult ipc_weighted_cover_reprice(const NsatCoverInput *input,
                                                  int lower, NsatLimits limits,
                                                  size_t min_branches, unsigned pivots);

/* Search with the root weights until a smaller cover has been materialized.
 * Re-price subsequent targets within the same chart and resource limits. */
IpcWeightedCoverResult ipc_weighted_cover_seed_then_reprice(const NsatCoverInput *input,
                                                            int lower, NsatLimits limits,
                                                            size_t min_branches, unsigned pivots);

typedef struct {
    NsatStatus status; /* UNSAT proves that the returned cover is minimum. */
    int terms;
    uint64_t nodes, repriced, extensions, tasks;
    unsigned threads; /* Actual proof workers; zero means the serial path. */
    NsatStats stats;
} IpcChartCoverResult;

/* Fixed PI chart: column-major coverage masks, ceil(rows / 64) words per PI.
 * solution is a feasible 0/1 incumbent, updated only with validated covers.
 * Unlike the ON/OFF interface, the supplied columns are the entire domain;
 * integer bounds and exhaustive branching prove optimality within this chart.
 * UNKNOWN/ERROR never certifies a minimum. All allocations obey limits. */
IpcChartCoverResult ipc_chart_cover(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots);

/* Extend each probe only after closing a root branch or validating a smaller
 * cover. limits.work is still a hard cap; a stalled probe returns UNKNOWN. */
IpcChartCoverResult ipc_chart_cover_adaptive(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots, uint64_t probe_work);

/* Disjoint warm proof subtrees, with private bases and a single total resource
 * allowance. Closed subtrees permit adaptive extensions. Caller polling stays
 * on the invoking thread. Without pthreads this uses the serial implementation. */
IpcChartCoverResult ipc_chart_cover_parallel(int rows, int columns,
    const unsigned long long *masks, int *solution, NsatLimits limits,
    size_t min_branches, unsigned pivots, uint64_t probe_work, unsigned threads);

#endif
