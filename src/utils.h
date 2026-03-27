#include <R_ext/Boolean.h>
#include "row_dominance.h"
#include "sort_cols.h"

// void printfarray3(int* arr, int nr, int nc);
// void printfarray4(double* arr, int nr, int nc);

int compute_value_bit_width(
    int max_value
);

void resize(
    void **array,
    int type, // 1 = int, 2 = unsigned int
    int increase,
    int size,
    int nrows
);

Rboolean too_complex(
    const unsigned int foundPI,
    const int solmin,
    const double maxcomb
);

void over_transpose(
    int matrix[],
    const int nr,
    const int nc,
    const int type // 0 boolean, 1 int, 2 double
);

Rboolean altb(
    double a,
    double b
);

Rboolean agteb(
    double a,
    double b
);

Rboolean redundant(
    unsigned int p_implicants_pos[],
    unsigned int p_implicants_val[],
    int implicant_words,
    unsigned int fixed_bits[],
    unsigned int value_bits[],
    unsigned int prevfoundPI,
    const int covered[],
    unsigned int covered_limit,
    bool debug
);

void increment(
    int k,
    int *e,
    int *h,
    int nconds,
    int *tempk,
    int minval
);

void populate_posneg(
    int *rowpos,
    int *rowneg,
    int nconds,
    int ttrows,
    int posrows,
    const int p_tt[],
    int posmat[],
    int negmat[],
    int *max_value
);

void get_noflevels(
    int noflevels[],
    const int p_tt[],
    int nconds,
    int ttrows
);

void get_decimals(
    const int posrows,
    const int negrows,
    const int k,
    int decpos[],
    int decneg[],
    const int posmat[],
    const int negmat[],
    const int tempk[],
    const int mbase[]
);

void get_uniques(
    const int posrows,
    int *found,
    int decpos[],
    Rboolean possiblePI[],
    int possiblePIrows[]
);

void verify_possible_PI(
    const int compare,
    const int negrows,
    int *found,
    Rboolean possiblePI[],
    const int possiblePIrows[],
    const int decpos[],
    const int decneg[]
);

void get_frows(
    int frows[],
    const Rboolean possiblePI[],
    const int possiblePIrows[],
    const int compare
);

void fill_matrix(
    int nrows,
    int ncols,
    int nofl[],
    int *matrix,
    int startrow,
    int cols[],
    int plus1
);

void calculate_rows(
    int *nrows,
    int ncols,
    int nofl[],
    int arrange,
    int maxprod
);

Rboolean all_covered(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols
);

long long unsigned int nchoosek(
    int n,
    int k
);
