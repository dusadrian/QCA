#include "cover_validation.h"

#include <stddef.h>

Rboolean qca_cover_is_feasible(
    const int *chart,
    int nrows,
    int ncols,
    const int *indices,
    int solmin
) {
    if (!chart || nrows <= 0 || ncols <= 0 || !indices ||
        solmin <= 0 || solmin > ncols) {
        return FALSE;
    }

    for (int i = 0; i < solmin; ++i) {
        int col = indices[i];
        if (col < 0 || col >= ncols) return FALSE;
        for (int j = 0; j < i; ++j) {
            if (indices[j] == col) return FALSE;
        }
    }

    for (int row = 0; row < nrows; ++row) {
        Rboolean covered = FALSE;
        for (int i = 0; i < solmin && !covered; ++i) {
            covered = chart[(size_t)indices[i] * (size_t)nrows + (size_t)row] != 0;
        }
        if (!covered) return FALSE;
    }
    return TRUE;
}
