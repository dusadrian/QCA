#include "qca_stopping.h"

#include <stdlib.h>

static bool binary_point_set_supported(
    const int posmat[],
    const int negmat[],
    const int noflevels[],
    int nconds,
    int posrows,
    int negrows,
    double picons,
    double solcons
) {
    if (
        !posmat || !negmat || !noflevels ||
        nconds <= 0 || posrows <= 0 || negrows <= 0 ||
        picons != 0.0 || solcons != 0.0
    ) {
        return false;
    }

    for (int c = 0; c < nconds; ++c) {
        if (noflevels[c] != 2) return false;
        for (int r = 0; r < posrows; ++r) {
            const int value = posmat[c * posrows + r];
            if (value < 0 || value > 1) return false;
        }
        for (int r = 0; r < negrows; ++r) {
            const int value = negmat[c * negrows + r];
            if (value < 0 || value > 1) return false;
        }
    }
    return true;
}

void qca_stopping_state_init(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    const int noflevels[],
    int nconds,
    int posrows,
    int negrows,
    double picons,
    double solcons
) {
    if (!state) return;
    *state = (QCAStoppingState){0};
    state->supported = binary_point_set_supported(
        posmat,
        negmat,
        noflevels,
        nconds,
        posrows,
        negrows,
        picons,
        solcons
    );
}

void qca_stopping_observe_coverage(
    QCAStoppingState *state,
    int level,
    bool covered
) {
    if (
        state && state->supported && covered &&
        state->coverage_horizon == 0
    ) {
        state->coverage_horizon = level;
    }
}

/*
 * Two positive rows are compatible when the cube formed by their common
 * coordinates covers no negative row.  QCA stores both matrices column-major.
 */
static bool on_pair_compatible(
    const int posmat[],
    const int negmat[],
    int nconds,
    int posrows,
    int negrows,
    int p,
    int q,
    int *agreements_out
) {
    int agreements = 0;
    for (int c = 0; c < nconds; ++c) {
        if (posmat[c * posrows + p] == posmat[c * posrows + q]) {
            ++agreements;
        }
    }
    if (agreements_out) *agreements_out = agreements;
    if (agreements == 0) return false;

    for (int z = 0; z < negrows; ++z) {
        bool matches = true;
        for (int c = 0; c < nconds; ++c) {
            const int pv = posmat[c * posrows + p];
            if (
                pv == posmat[c * posrows + q] &&
                negmat[c * negrows + z] != pv
            ) {
                matches = false;
                break;
            }
        }
        if (matches) return false;
    }
    return true;
}

static bool pair_joint_in_chart(
    const int pichart[],
    int posrows,
    int foundPI,
    int p,
    int q
) {
    for (int c = 0; c < foundPI; ++c) {
        const int offset = c * posrows;
        if (pichart[offset + p] && pichart[offset + q]) return true;
    }
    return false;
}

bool qca_stopping_observe_plateau(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    int nconds,
    int posrows,
    int negrows,
    const int pichart[],
    int foundPI,
    const int selected_indices[],
    int selected_terms
) {
    if (!state || !state->supported || state->diagnostic_checked) return true;
    if (
        !posmat || !negmat || !pichart || !selected_indices ||
        nconds <= 0 || posrows <= 0 || negrows <= 0 || foundPI <= 0 ||
        selected_terms <= 0
    ) {
        return false;
    }

    int *witness = (int *)calloc((size_t)selected_terms, sizeof(int));
    if (!witness) return false;
    for (int j = 0; j < selected_terms; ++j) witness[j] = -1;

    for (int r = 0; r < posrows; ++r) {
        int owner = -1;
        int covering_terms = 0;
        for (int j = 0; j < selected_terms; ++j) {
            const int c = selected_indices[j];
            if (c < 0 || c >= foundPI) {
                free(witness);
                return false;
            }
            if (pichart[c * posrows + r]) {
                owner = j;
                ++covering_terms;
            }
        }
        if (covering_terms == 1 && witness[owner] < 0) witness[owner] = r;
    }

    bool delayed_pair = false;
    for (int a = 0; a < selected_terms && !delayed_pair; ++a) {
        if (witness[a] < 0) continue;
        for (int b = a + 1; b < selected_terms; ++b) {
            if (witness[b] < 0) continue;
            if (
                on_pair_compatible(
                    posmat,
                    negmat,
                    nconds,
                    posrows,
                    negrows,
                    witness[a],
                    witness[b],
                    NULL
                ) &&
                !pair_joint_in_chart(
                    pichart,
                    posrows,
                    foundPI,
                    witness[a],
                    witness[b]
                )
            ) {
                delayed_pair = true;
                break;
            }
        }
    }

    free(witness);
    state->diagnostic_checked = true;
    state->certification_required = delayed_pair;
    return true;
}

bool qca_stopping_prepare_certificate(
    QCAStoppingState *state,
    const int posmat[],
    const int negmat[],
    int nconds,
    int posrows,
    int negrows
) {
    if (!state || !state->supported || !posmat || !negmat) return false;
    if (state->cover_lower_bound > 0) return true;

    int horizon = 0;
    for (int p = 0; p < posrows; ++p) {
        for (int q = p + 1; q < posrows; ++q) {
            int agreements = 0;
            if (
                on_pair_compatible(
                    posmat,
                    negmat,
                    nconds,
                    posrows,
                    negrows,
                    p,
                    q,
                    &agreements
                ) && agreements > horizon
            ) {
                horizon = agreements;
            }
        }
    }

    int *selected = (int *)calloc((size_t)posrows, sizeof(int));
    if (!selected) return false;

    int selected_count = 0;
    for (int p = 0; p < posrows; ++p) {
        bool incompatible_with_all = true;
        for (int j = 0; j < selected_count; ++j) {
            if (on_pair_compatible(
                posmat,
                negmat,
                nconds,
                posrows,
                negrows,
                p,
                selected[j],
                NULL
            )) {
                incompatible_with_all = false;
                break;
            }
        }
        if (incompatible_with_all) selected[selected_count++] = p;
    }
    free(selected);

    state->agreement_horizon = horizon;
    state->cover_lower_bound = selected_count;
    return true;
}

bool qca_stopping_cardinality_certified(
    QCAStoppingState *state,
    int level,
    int cover_size,
    bool boundary_exact
) {
    if (!state || !state->supported || cover_size <= 0) return false;
    if (cover_size == 1) state->cardinality_certified = true;
    if (
        state->cover_lower_bound > 0 &&
        cover_size == state->cover_lower_bound
    ) {
        state->cardinality_certified = true;
    }

    int horizon = state->agreement_horizon;
    if (state->coverage_horizon > horizon) horizon = state->coverage_horizon;
    if (
        boundary_exact && state->cover_lower_bound > 0 &&
        state->coverage_horizon > 0 && level >= horizon
    ) {
        state->cardinality_certified = true;
    }
    return state->cardinality_certified;
}
