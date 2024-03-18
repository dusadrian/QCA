#include <R_ext/Boolean.h>
#include "super_rows.h"
#include "utils.h"
#include "consistent_solution.h"

void find_consistent_models(
    const int p_implicants[],
    const int indx[],
    const int ck[],
    const double p_data[],
    const int p_fuzzy[],
    const int nconds,
    const int nrdata,
    const int posrows,
    const double solcons,
    const double solcov,
    const Rboolean allsol,
    const int soldepth,
    const unsigned int foundPI,
    const double maxcomb,
    
    int **solutions,
    int *nr,
    int *nc
);
