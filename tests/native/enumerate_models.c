#include "enumerate_models.h"

SEXP test_enumerate_fixed(SEXP chart, SEXP size) {
    int *solutions = R_Calloc(1, int), nr = 0, nc = 0;
    qca_enumerate_fixed_models(LOGICAL(chart), nrows(chart), ncols(chart),
                               asInteger(size), &solutions, &nr, &nc);
    SEXP out = PROTECT(allocMatrix(INTSXP, nr, nc));
    if (nr && nc) memcpy(INTEGER(out), solutions, (size_t)nr * nc * sizeof(int));
    R_Free(solutions);
    UNPROTECT(1);
    return out;
}
