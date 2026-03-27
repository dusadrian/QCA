#include <stdbool.h>
#include <Rinternals.h>

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
