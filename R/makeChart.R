# Copyright (c) 2016 - 2024, Adrian Dusa
# All rights reserved.
# 
# Redistribution and use in source and binary forms, with or without
# modification, in whole or in part, are permitted provided that the
# following conditions are met:
#     * Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * The names of its contributors may NOT be used to endorse or promote
#       products derived from this software without specific prior written
#       permission.
# 
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ADRIAN DUSA BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
# (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
# LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

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
        mtrx <- matrix(
            FALSE,
            nrow = length(tprimes),
            ncol = length(tconfigs)
        )
        for (i in seq(nrow(mtrx))) {
            for (j in seq(ncol(mtrx))) {
                subset <- TRUE
                s <- 1
                while (subset & s <= length(tprimes[[i]])) {
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
