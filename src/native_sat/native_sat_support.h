#ifndef CCUBES_NATIVE_SAT_SUPPORT_H
#define CCUBES_NATIVE_SAT_SUPPORT_H

#include "native_sat_verify.h"
#include <limits.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    NsatLimits limits;
    NsatStats stats;
    NsatStatus status;
    unsigned until_poll;
} NsatBudget;

static inline bool nsat_stop(NsatBudget *b, NsatStatus status, NsatStopReason why)
{
    if (b->status == NSAT_OK) {
        b->status = status;
        b->stats.stop = why;
        if (b->limits.shared && !b->limits.shared->stop)
            b->limits.shared->stop = why;
    }
    return false;
}

static inline bool nsat_poll(NsatBudget *b)
{
    if (b->status != NSAT_OK) return false;
    b->until_poll = 4096;
    if (b->limits.poll && b->limits.poll(b->limits.context))
        return nsat_stop(b, NSAT_UNKNOWN, NSAT_STOP_CANCELLED);
    return true;
}

static inline bool nsat_work(NsatBudget *b, uint64_t units)
{
    if (b->status != NSAT_OK) return false;
    uint64_t max = b->limits.work ? b->limits.work : UINT64_MAX;
    NsatResources *shared = b->limits.shared;
    uint64_t overall = shared && shared->work_limit ? shared->work_limit : UINT64_MAX;
    if (units > max - b->stats.work || (shared && units > overall - shared->work))
        return nsat_stop(b, NSAT_UNKNOWN, NSAT_STOP_WORK);
    b->stats.work += units;
    if (shared) shared->work += units;
    if (units >= b->until_poll && !nsat_poll(b)) return false;
    b->until_poll = units >= b->until_poll ? 0 : b->until_poll - (unsigned)units;
    return true;
}

static inline bool nsat_mul(size_t a, size_t b, size_t *out)
{
    if (b && a > SIZE_MAX / b) return false;
    *out = a * b;
    return true;
}

/* Optional scratch space must leave the mandatory allocation path usable. */
static inline size_t nsat_memory_available(const NsatBudget *b)
{
    size_t available = b->limits.memory_bytes
        ? b->limits.memory_bytes - b->stats.memory_bytes : SIZE_MAX;
    const NsatResources *shared = b->limits.shared;
    if (shared && shared->memory_limit) {
        size_t remaining = shared->memory_limit - shared->memory;
        if (remaining < available) available = remaining;
    }
    return available;
}

static inline void *nsat_alloc(NsatBudget *b, size_t count, size_t width)
{
    size_t bytes;
    if (!nsat_mul(count, width, &bytes)) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_INPUT);
        return NULL;
    }
    size_t max = b->limits.memory_bytes ? b->limits.memory_bytes : SIZE_MAX;
    NsatResources *shared = b->limits.shared;
    size_t overall = shared && shared->memory_limit ? shared->memory_limit : SIZE_MAX;
    if (bytes > max - b->stats.memory_bytes || (shared && bytes > overall - shared->memory)) {
        nsat_stop(b, NSAT_UNKNOWN, NSAT_STOP_MEMORY);
        return NULL;
    }
    if (!nsat_poll(b)) return NULL;
#ifdef NSAT_TESTING
    extern void *nsat_test_malloc(size_t);
    void *ptr = nsat_test_malloc(bytes);
#else
    void *ptr = malloc(bytes);
#endif
    if (!ptr) {
        nsat_stop(b, NSAT_ERROR, NSAT_STOP_ALLOCATION);
        return NULL;
    }
    b->stats.memory_bytes += bytes;
    if (b->stats.peak_memory_bytes < b->stats.memory_bytes)
        b->stats.peak_memory_bytes = b->stats.memory_bytes;
    if (shared) {
        shared->memory += bytes;
        if (shared->peak_memory < shared->memory) shared->peak_memory = shared->memory;
    }
    return ptr;
}

static inline void nsat_free(NsatBudget *b, void *ptr, size_t bytes)
{
    if (ptr) {
        free(ptr);
        b->stats.memory_bytes -= bytes;
        if (b->limits.shared) b->limits.shared->memory -= bytes;
    }
}

static inline void nsat_release_budget(NsatBudget *b)
{
    if (b->limits.shared) b->limits.shared->memory -= b->stats.memory_bytes;
}

/* Charge and poll during bulk initialization and copying as well as scans. */
static inline bool nsat_copy(NsatBudget *b, void *dest, const void *src, size_t bytes)
{
    unsigned char *d = dest;
    const unsigned char *s = src;
    while (bytes) {
        size_t chunk = bytes < 4096 ? bytes : 4096;
        if (!nsat_work(b, chunk)) return false;
        if (s) {
            memcpy(d, s, chunk);
            s += chunk;
        } else {
            memset(d, 0, chunk);
        }
        d += chunk;
        bytes -= chunk;
    }
    return true;
}

static inline int nsat_var(int lit, int variables)
{
    if (!lit || lit > variables || lit < -variables) return 0;
    return lit < 0 ? -lit : lit;
}

#endif
