`findSupersets` <- function (
    input, noflevels = NULL, ...
) {
    
    dots <- list(...)
    
    ### 
    ### ### backwards compatibility 
    ### 
        if (is.element("input.combs", names(dots))) {
            input <- dots$input.combs
        }
    ### 
    ### ### backwards compatibility 
    ### 
    
    if (!is.matrix(input)) {
        if (!is.vector(drop(input))) {
            admisc::stopError(
                "input must be either an solution-space matrix or a vector of row numbers.",
                ... = ...
            )
        }
        else {
            if (any(input > prod(noflevels))) {
                admisc::stopError(
                    paste(
                        "Some line numbers do not belong in the solution-space for",
                        length(noflevels),
                        "causal conditions."
                    ),
                    ... = ...
                )
            }
            input <- getRow(input, noflevels)
        }
    }
    
    mbase <- rev(c(1, cumprod(rev(noflevels))))[-1]
    allcombn <- t(createMatrix(rep(2, length(noflevels)))[-1, ])
    
    primes <- sort.int(
        unique(
            as.vector(
                apply(input, 1, function(x) (x*mbase) %*% allcombn + 1)
            )
        )
    )

    if (primes[1] == 1) {
        return(primes[-1])
    }
    else {
        return(primes)
    }
}
