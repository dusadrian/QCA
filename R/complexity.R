complexity <- function(
    n, layers = NULL, noflevels = NULL, ...
) {
    if (
        !is.numeric(n) || length(n) != 1L || !is.finite(n) ||
        n < 0 || n > .Machine$integer.max || n != floor(n)
    ) {
        admisc::stopError("Argument <n> must be a nonnegative whole number.", ... = ...)
    }

    if (is.null(noflevels)) {
        noflevels <- rep(2L, n)
    }
    if (is.null(layers)) {
        layers <- seq_len(n)
    }
    if (
        !is.numeric(noflevels) || length(noflevels) != n ||
        any(!is.finite(noflevels)) || any(noflevels != floor(noflevels)) ||
        any(noflevels < 1 | noflevels > .Machine$integer.max)
    ) {
        admisc::stopError("Specify a positive number of levels for each condition.", ... = ...)
    }
    if (
        !is.numeric(layers) || any(!is.finite(layers)) ||
        any(layers != floor(layers)) || any(layers < 1 | layers > n)
    ) {
        admisc::stopError("Specify whole complexity layers between 1 and n.", ... = ...)
    }

    return(.Call(
        "C_omplexity", list(as.integer(n), as.integer(layers), as.integer(noflevels)),
        PACKAGE = "QCA"
    ))
}
