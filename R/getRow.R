`getRow` <- function(
    row.no = NULL, noflevels = NULL, zerobased = FALSE, ...
) {
    if (
        anyNA(noflevels) || any(noflevels < 1 | noflevels > .Machine$integer.max) ||
        any(noflevels != floor(noflevels))
    ) {
        admisc::stopError("Condition levels must be positive whole numbers.", ... = ...)
    }

    if (is.numeric(row.no) && any(row.no > prod(noflevels) - zerobased, na.rm = TRUE)) {
        admisc::stopError(
            paste0("There cannot be more than ", prod(noflevels), " rows."), ... = ...
        )
    }

    return(.Call("C_getRow", list(row.no, noflevels, as.logical(zerobased)), PACKAGE = "QCA"))
}

configurationRows <- function(data, noflevels) {
    # R's ordinary numeric row labels retain up to 15 significant digits.
    if (prod(noflevels) >= 1e15) {
        return(.Call("C_configurationRows", as.matrix(data), noflevels, PACKAGE = "QCA"))
    }

    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    return(as.vector(as.matrix(data) %*% mbase) + 1)
}
