#include <assert.h>
#include <stdatomic.h>
#include <stdlib.h>

#include "qca_threads.h"

typedef struct {
    atomic_int *visits;
    atomic_int *worker_ids;
} test_context;

static void visit_range(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
) {
    test_context *context = (test_context *)data;
    atomic_store_explicit(
        &context->worker_ids[worker_id], 1, memory_order_relaxed
    );
    for (unsigned long long i = start; i < end; i++) {
        atomic_fetch_add_explicit(
            &context->visits[i], 1, memory_order_relaxed
        );
    }
}

static void verify_parallel_for(int threads, int count) {
    atomic_int *visits = calloc((size_t)count, sizeof(atomic_int));
    atomic_int *worker_ids = calloc((size_t)threads, sizeof(atomic_int));
    assert(visits != NULL);
    assert(worker_ids != NULL);

    test_context context = {
        .visits = visits,
        .worker_ids = worker_ids
    };
    assert(qca_parallel_for(
        (unsigned long long)count, threads, visit_range, &context
    ));

    for (int i = 0; i < count; i++) {
        assert(atomic_load_explicit(&visits[i], memory_order_relaxed) == 1);
    }
    int expected_worker_ids = 1;
#if defined(HAVE_PTHREAD)
    expected_worker_ids = threads;
#endif
    for (int i = 0; i < expected_worker_ids; i++) {
        assert(atomic_load_explicit(
            &worker_ids[i], memory_order_relaxed
        ) == 1);
    }

    free(visits);
    free(worker_ids);
}

int main(void) {
    verify_parallel_for(1, 17);
    verify_parallel_for(2, 101);
    verify_parallel_for(4, 1003);
    return 0;
}
