`makeChart` <- function(
    primes = "", configs = "", snames = "", mv = FALSE, collapse = "*", ...
) {


    primes <- admisc::recreate(substitute(primes))
    configs <- admisc::recreate(substitute(configs))
    snames <- admisc::recreate(substitute(snames))

    prmat <- is.matrix(primes)
    comat <- is.matrix(configs)
    
    dots <- list(...)

    curly <- ifelse(is.null(dots$curly), FALSE, dots$curly)

    if (prmat & comat) {

        if (!(is.numeric(primes) & is.numeric(configs))) {
            admisc::stopError(
                "Matrices have to be numeric.", ... = ...
            )
        }

        if (any(primes < 0) | any(configs < 0)) {
            admisc::stopError(
                "Matrix values have to be non-negative.", ... = ...
            )
        }

        ### bug fix 22.09.2017, when function minimize returns a complete Boolean elimination
        # due to a complete truth table being fed into the QMS for c.sol
        # example from Michael Baumgartner
        # if (!is.element("getSolution", unlist(lapply(lapply(sys.calls(), as.character), "[[", 1)))) {
        if (
            !is.element("getSolution", names(dots)) &&
            (
                any(
                    apply(primes, 1, sum) == 0
                ) |
                any(
                    apply(configs, 1, sum) == 0
                )
            )
        ) {
            admisc::stopError(
                "Matrices have to be specified at implicants level.", ... = ...
            )
            
        }

        if (nrow(primes == 1) & sum(primes) == 0) {
            mtrx = matrix(nrow = 0, ncol = nrow(configs))
        }
        else {
            primes2 <- matrix(logical(length(primes)), dim(primes))
            primes2[primes > 0] <- TRUE

            mtrx <- sapply(seq(nrow(primes)), function(x) {
                apply(configs, 1, function(y) {
                    all(primes[x, primes2[x, ]] == y[primes2[x, ]])
                })
            })

            if (nrow(configs) == 1) {
                mtrx <- matrix(mtrx)
            }
            else {
                mtrx <- t(mtrx)
            }

            rownames(mtrx) <- admisc::writePrimeimp(
                primes,
                mv = mv,
                collapse = collapse,
                curly = curly
            )
        }

        colnames(mtrx) <- admisc::writePrimeimp(
            configs,
            mv = mv,
            collapse = collapse,
            curly = curly
        )

        class(mtrx) <- c("matrix", "pic")
        return(mtrx)
        
    }
    else if (!prmat & !comat) {
        
        if (!identical(snames, "")) {
            if (length(snames) == 1 & is.character(snames)) {
                snames <- admisc::splitstr(snames)
            }
        }

        noflevels <- dots$noflevels
        if (!identical(snames, "") && is.null(noflevels)) {
            noflevels <- rep(2, length(snames))
        }

        # noflevels is redundant, but solves a bug with the GUI, via findRows/minimize
        tconfigs <- attr(
            admisc::translate(
                expression = configs,
                snames = snames,
                noflevels = noflevels,
                retlist = TRUE
            ),
            "retlist"
        )

        if (identical(snames, "")) {
            snames <- names(tconfigs[[1]])
        }

        tprimes <- attr(
            admisc::translate(
                expression = primes,
                snames = snames,
                noflevels = noflevels,
                retlist = TRUE
            ),
            "retlist"
        )

        # both primes and configs cannot contain multiple levels

        mtrx <- matrix(
            FALSE,
            nrow = length(tprimes),
            ncol = length(tconfigs)
        )

        for (i in seq(nrow(mtrx))) {
            # so we can unlist
            # tp <- unlist(tprimes[[i]])

            for (j in seq(ncol(mtrx))) {
                # tc <- unlist(tconfigs[[j]])
                subset <- TRUE
                s <- 1
                while (subset & s <= length(tprimes[[i]])) {
                    # all(tp[tp >= 0] == tc[tp >= 0])
                    if (tprimes[[i]][[s]] >= 0) {
                        subset <- is.element(
                            tprimes[[i]][[s]],
                            tconfigs[[j]][[s]]
                        )
                    }
                    s <- s + 1
                }
                mtrx[i, j] <- subset
            }
        }

        colnames(mtrx) <- names(tconfigs)
        rownames(mtrx) <- names(tprimes)

        class(mtrx) <- c("matrix", "QCA_pic")
        return(mtrx)
    }
    else {
        admisc::stopError(
            "Both arguments have to be matrices.", ... = ...
        )
    }
}
