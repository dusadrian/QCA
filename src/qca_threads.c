#include "qca_threads.h"

#include <stdlib.h>

#if defined(HAVE_PTHREAD)
    #include <pthread.h>
    #if defined(_WIN32)
        #include <windows.h>
    #else
    #include <unistd.h>
    #endif
#endif

#define QCA_THREAD_LIMIT 64

typedef struct {
    qca_range_worker worker;
    void *data;
    unsigned long long start;
    unsigned long long end;
    int worker_id;
} QCARangeTask;

#if defined(HAVE_PTHREAD)
static void *qca_range_thread_main(void *arg) {
    QCARangeTask *task = (QCARangeTask *) arg;
    task->worker(task->start, task->end, task->worker_id, task->data);
    return NULL;
}
#endif

int qca_default_thread_count(void) {
#if defined(HAVE_PTHREAD)
    #if defined(_WIN32)
    SYSTEM_INFO info;
    GetSystemInfo(&info);
    if (info.dwNumberOfProcessors > 1) {
        return info.dwNumberOfProcessors > QCA_THREAD_LIMIT ?
            QCA_THREAD_LIMIT : (int) info.dwNumberOfProcessors;
    }
    #else
    long nprocs = sysconf(_SC_NPROCESSORS_ONLN);
    if (nprocs > 1) {
        return nprocs > QCA_THREAD_LIMIT ? QCA_THREAD_LIMIT : (int) nprocs;
    }
    #endif
#endif
    return 1;
}

int qca_parallel_for(
    unsigned long long count,
    int requested_threads,
    qca_range_worker worker,
    void *data
) {
    int nthreads = requested_threads;

    if (count == 0) {
        return 1;
    }

    if (nthreads <= 0) {
        nthreads = qca_default_thread_count();
    }
    if (nthreads < 1) {
        nthreads = 1;
    }
    if (nthreads > QCA_THREAD_LIMIT) {
        nthreads = QCA_THREAD_LIMIT;
    }
    if ((unsigned long long) nthreads > count) {
        nthreads = (int) count;
    }

#if defined(HAVE_PTHREAD)
    if (nthreads > 1) {
        pthread_t *threads = (pthread_t *) calloc((size_t) nthreads, sizeof(pthread_t));
        QCARangeTask *tasks = (QCARangeTask *) calloc((size_t) nthreads, sizeof(QCARangeTask));
        unsigned long long base = count / (unsigned long long) nthreads;
        unsigned long long rem = count % (unsigned long long) nthreads;
        unsigned long long start = 0;
        int started = 0;
        int ok = 1;

        if (threads == NULL || tasks == NULL) {
            free(threads);
            free(tasks);
            return 0;
        }

        for (int i = 0; i < nthreads; i++) {
            unsigned long long width = base + ((unsigned long long) i < rem ? 1 : 0);
            tasks[i].worker = worker;
            tasks[i].data = data;
            tasks[i].start = start;
            tasks[i].end = start + width;
            tasks[i].worker_id = i;
            start += width;

            if (pthread_create(&threads[i], NULL, qca_range_thread_main, &tasks[i]) != 0) {
                ok = 0;
                break;
            }
            started++;
        }

        for (int i = 0; i < started; i++) {
            pthread_join(threads[i], NULL);
        }

        free(threads);
        free(tasks);
        return ok;
    }
#endif

    worker(0, count, 0, data);
    return 1;
}

int qca_mutex_init(qca_mutex *mutex) {
#if defined(HAVE_PTHREAD)
    return pthread_mutex_init(&mutex->mutex, NULL) == 0;
#else
    (void) mutex;
    return 1;
#endif
}

void qca_mutex_lock(qca_mutex *mutex) {
#if defined(HAVE_PTHREAD)
    pthread_mutex_lock(&mutex->mutex);
#else
    (void) mutex;
#endif
}

void qca_mutex_unlock(qca_mutex *mutex) {
#if defined(HAVE_PTHREAD)
    pthread_mutex_unlock(&mutex->mutex);
#else
    (void) mutex;
#endif
}

void qca_mutex_destroy(qca_mutex *mutex) {
#if defined(HAVE_PTHREAD)
    pthread_mutex_destroy(&mutex->mutex);
#else
    (void) mutex;
#endif
}
