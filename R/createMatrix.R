`createMatrix` <- function(
    noflevels = NULL, ...
) {
    
    dots <- list(...)
    
    RAM <- 2
    if (is.element("RAM", names(dots))) {
        if (length(dots$RAM) == 1) {
            if (is.numeric(dots$RAM) & dots$RAM > 0) {
                RAM <- dots$RAM
            }
        }
    }
    
    arrange <- FALSE
    if (is.element("arrange", names(dots))) {
        arrange <- dots$arrange
    }
    
    depth <- length(noflevels)
    if (is.element("depth", names(dots))) {
        if (!is.null(dots$depth)) {
            if (is.numeric(dots$depth)) {
                depth <- dots$depth
                if (depth < length(noflevels)) {
                    arrange <- TRUE
                }
            }
        }
    }
    
    if (any(abs(noflevels) %% 1 > .Machine$double.eps ^ 0.5)) {
        admisc::stopError(
            "The number of levels must be integers.", ... = ...
        )
    }
    
    if (!is.logical(arrange)) {
        admisc::stopError(
            "The argument <arrange> should be logical.", ... = ...
        )
    }
    
    if (abs(depth) %% 1 > .Machine$double.eps ^ 0.5) {
        admisc::stopError(
            "The argument depth has to be an integer number.", ... = ...
        )
    }
    
    if ((mem <- prod(noflevels) * length(levels) * 8 / 1024^3) > RAM) {
        admisc::stopError(
            paste0(
                "Too much memory needed (",
                round(mem, 1),
                " Gb) to create the matrix."
            ),
            ... = ...
        )
    }
    
    noflevels <- as.integer(abs(noflevels))
    arrange <- as.integer(arrange * 1)
    depth <- as.integer(abs(depth))
    nofconds <- as.integer(length(noflevels))
    
    if (arrange) {
        if (depth < 1 | depth > nofconds) {
            depth <- nofconds
        }
    }
    
    tosend <- list(noflevels, arrange, depth)
    
    if (is.element("colnames", names(dots))) {
        colnms <- dots$colnames
        if (is.character(colnms)) {
            if (length(colnms) == length(noflevels)) {
                tosend <- c(tosend, list(colnms))
            }
        }
    }
    
    return(.Call("C_createMatrix", tosend, PACKAGE = "QCA"))
    
    pwr <- unique(noflevels)
    
    # rep.int comes from rep internal, not rep integers
    
    if (length(pwr) == 1) {
        create <- function(idx) {
            rep.int(c(sapply(seq_len(pwr) - 1, function(x) rep.int(x, pwr^(idx - 1)))),
                    pwr^nofconds/pwr^idx)
        }
        retmat <- sapply(rev(seq_len(nofconds)), create)
    }
    else {
        mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
        orep  <- cumprod(rev(c(rev(noflevels)[-1], 1)))
        retmat <- sapply(seq_len(nofconds), function(x) {
           rep.int(rep.int(seq_len(noflevels[x]) - 1, rep.int(mbase[x], noflevels[x])), orep[x])
        })
    }
    
    if (is.vector(retmat)) {
        retmat <- matrix(retmat, nrow = 1)
    }
    
    return(retmat)
}
