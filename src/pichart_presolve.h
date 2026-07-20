#ifndef QCA_PICHART_PRESOLVE_H
#define QCA_PICHART_PRESOLVE_H

int qca_reduce_active_columns(
    const int *chart,
    int nrows,
    int ncols,
    unsigned char *active
);

#endif
