#ifndef QCA_THREADS_H
#define QCA_THREADS_H

#if defined(HAVE_PTHREAD)
    #include <pthread.h>
#endif

typedef void (*qca_range_worker)(
    unsigned long long start,
    unsigned long long end,
    int worker_id,
    void *data
);

int qca_default_thread_count(void);

int qca_parallel_for(
    unsigned long long count,
    int requested_threads,
    qca_range_worker worker,
    void *data
);

typedef struct {
#if defined(HAVE_PTHREAD)
    pthread_mutex_t mutex;
#else
    int unused;
#endif
} qca_mutex;

int qca_mutex_init(qca_mutex *mutex);
void qca_mutex_lock(qca_mutex *mutex);
void qca_mutex_unlock(qca_mutex *mutex);
void qca_mutex_destroy(qca_mutex *mutex);

#endif
