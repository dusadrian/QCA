#include <stdbool.h>
#include "qca_rinternals.h"

bool gurobi_runtime_available(void);
void gurobi_release_env(void);
SEXP C_gurobiRuntimeAvailable(void);
SEXP C_findminExact(SEXP chart);

bool solvePIchart_gurobi(
    const int pichart[],
    int foundPI,
    int on_minterms,
    int indices[],
    int *solmin
);
bool solvePIchart_gurobi_active(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const unsigned char active[],
    int indices[],
    int *solmin
);
bool solvePIchart_gurobi_active_with_incumbent(
    const int pichart[],
    int foundPI,
    int on_minterms,
    const unsigned char active[],
    const int initial_indices[],
    int initial_solmin,
    int indices[],
    int *solmin
);
