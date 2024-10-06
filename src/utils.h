#include <R_ext/Boolean.h>
#include "row_dominance.h"
#include "sort_cols.h"

// void printfarray3(int* arr, int nr, int nc);
// void printfarray4(double* arr, int nr, int nc);

void resize(
    int **array,
    const int rows,
    const unsigned int newcols,
    const unsigned int oldcols
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

Rboolean nonredundant(
    const int p_implicants[],
    const int p_indx[],
    const int p_ck[],
    const int tempk[],
    const int tempc[],
    const int nconds,
    const int k,
    const unsigned int prevfoundPI
);

void push_PI(
    int p_implicants[],
    int p_indx[],
    int p_ck[],
    int p_pichart[],
    const int tempk[],
    const int tempc[],
    const int nconds,
    const int k,
    const int f,
    const int decpos[],
    const int frows[],
    const int posrows,
    const unsigned int foundPI
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
    int negmat[]
);

void get_noflevels(
    int noflevels[],
    const int p_tt[],
    int nconds,
    int ttrows
);

void fill_mbase(
    int mbase[],
    const int tempk[],
    const int noflevels[],
    const int k
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

int fill_combination_tasks(
    const int n,
    int k,
    int (*comb)[3],
    const int arrayMaxSize,
    int* outNumItemsFilled
);

int get_next_combination(
    int *arr,
    int k,
    int n,
    int indexLimit
);

Rboolean all_covered(
    const int p_pichart[],
    const int pirows,
    const unsigned int picols
);
