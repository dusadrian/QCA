complexity <- function(
    n, layers = NULL, noflevels = NULL, ...
) {
    if (!is.numeric(n)) {
        admisc::stopError("Argument <n> should be numeric.", ... = ...)
    }
    
    if (length(n) != 1L) {
        admisc::stopError("Argument <n> should be a scalar of length 1.", ... = ...)
    }
    
    if (n < 0) {
        admisc::stopError("Argument <n> should be positive.", ... = ...)
    }

    if (is.null(noflevels)) noflevels <- rep(2, n)
    if (is.null(layers)) layers <- seq(n)
    
    if (any(layers > n)) {
        admisc::stopError("Argument <layers> cannot be greater than <n>.", ... = ...)
    }
    
    sumk <- .Call("C_omplexity", list(as.integer(n), as.integer(layers), as.integer(noflevels)), PACKAGE = "QCA")
    sumk[sumk < 0] <- Inf
    
    return(sumk)
    
    # the same code in pure R
    sumk <- rep(0, length(layers))
    for (i in seq(length(layers))) {
        sumk[i] <- sum(apply(admisc::combnk(n, layers[i]), 2, function(x) {
            prod(noflevels[x])
        }))
    }
}
