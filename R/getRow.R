`getRow` <-  function(
    row.no = NULL, noflevels = NULL, zerobased = FALSE, ...
) {
    
    dots <- list(...)
    enter <- ifelse (is.element("enter", names(dots)), dots$enter, TRUE)
    
    max.combs <- prod(noflevels)
    if (any(row.no > (max.combs - zerobased))) {
        admisc::stopError(
            paste0(
                "There cannot be more than ",
                max.combs,
                " rows."
            ),
            ... = ...
        )
    }
    
    if (!zerobased) {row.no <- row.no - 1}
    
    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    
    return(
        .Call(
            "C_getRow",
            list(
                row.no,
                noflevels,
                mbase
            ),
            PACKAGE = "QCA"
        )
    )
    
    # t(sapply(row.no, function(x) x %/% mbase) %% noflevels)
}
