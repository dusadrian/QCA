#include <R_ext/Boolean.h>

Rboolean consistent_solution(
    const double p_data[],      // the raw, calibrated (fuzzy) data
    const int ncdata,           // number of columns in the data
    const int nrdata,           // number of rows in the data
    const int k,                // level of complexity
    const int tempk[],          // combination of PIs at the level of complexity k
    const unsigned int foundPI, // the number of found PIs
    const int p_implicants[],   // the PIs implicants matrix
    const int ck[],             // the vector of complexity levels where the PIs were found
    const int indx[],           // the indices of the PIs in the implicants matrix, at each position in ck
    const int p_fsconds[],      // boolean, if conditions fuzzy
    const double solcons,       // the solution consistency threshold
    const double solcov         // the solution coverage threshold
);
