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

`findRows` <- function(
    obj = NULL, expression = "", observed = FALSE, type = 1, ...
) {
    expression <- admisc::recreate(substitute(expression))
    dots <- list(...)
    if (is.element("remainders", names(dots))) {
        if (is.logical(dots$remainders)) {
            observed <- !dots$remainders[1]
        }
    }
    if (any(type == 0)) {
        type <- 0
    }
    if (any(is.element(type, 0:1)) & identical(expression, "")) {
        admisc::stopError(
            "The expression is missing, for type 0 or 1.", ... = ...
        )
    }
    noflevels <- NULL
    if (is.element("noflevels", names(dots))) {
        noflevels <- dots$noflevels
    }
    conditions <- NULL
    if (is.element("conditions", names(dots))) {
        conditions <- dots$conditions
    }
    if (missing(obj) && !is.null(noflevels)) {
        admisc::stopError(
            "The truth table object is missing.", ... = ...
        )
    }
    if (methods::is(obj, "QCA_tt")) {
        noflevels <- obj$noflevels
        conditions <- obj$options$conditions
        if (any(is.element(type, c(0, 2, 3)))) {
            call <- as.list(obj$call)[-1]
            call$data <- obj$initial.data
            call$outcome <- obj$options$outcome
            if (admisc::tilde1st(call$outcome)) {
                call$outcome <- admisc::notilde(call$outcome)
            }
            else {
                call$outcome <- paste("~", call$outcome, sep = "")
            }
            call$incl.cut <- rev(obj$options$incl.cut)
            if (length(dots) > 0) {
                if (length(setdiff(names(dots), c("incl.cut", "n.cut", "pri.cut"))) > 0) {
                    admisc::stopError(
                        "Only cutoff arguments can be specified for the negation of the outcome.",
                        ... = ...
                    )
                }
                nms <- names(dots)
                for (i in seq(length(nms))) {
                    call[[nms[i]]] <- dots[[nms[i]]]
                }
            }
            for (i in seq(length(call))) {
                call[[i]] <- eval.parent(call[[i]])
            }
            nobj <- suppressWarnings(do.call("truthTable", call))
        }
    }
    else {
        if (is.null(obj)) {
            if (is.null(noflevels) | is.null(conditions)) {
                admisc::stopError(
                    "The truth table argument <obj> is missing.", ... = ...
                )
            }
        }
        else {
            if (is.matrix(obj) && is.numeric(obj)) {
                conditions <- colnames(obj)
                if (is.null(conditions)) {
                    admisc::stopError(
                        "The <obj> matrix does not have column names.", ... = ...
                    )
                }
            }
            else {
                admisc::stopError(
                    "Argument <obj> is not a truth table object or a numerical matrix.",
                    ... = ...
                )
            }
        }
    }
    SBS <- NULL
    CSA <- NULL
    SSR <- NULL
    if (any(is.element(type, 0:1))) {
        trexp <- attr(
            admisc::translate(
                paste(expression, collapse = "+"),
                snames = conditions,
                retlist = TRUE
            ),
            "retlist"
        )
        result <- matrix(ncol = length(conditions), nrow = 0)
        if (is.matrix(obj)) {
            noflevels <- admisc::getInfo(obj)$noflevels
        }
        for (i in seq(length(trexp))) {
            rowi <- trexp[[i]]
            detected <- !unlist(lapply(rowi, function(x) identical(x, -1)))
            rowi <- rowi[detected]
            rowi <- expand.grid(rowi)
            if (sum(!detected) > 0) {
                restm <- createMatrix(noflevels[!detected])
                colnames(restm) <- conditions[!detected]
                rowi <- apply(rowi, 1, function(x) rep(x, each = nrow(restm)))
                for (r in seq(ncol(rowi))) {
                    detm <- matrix(rowi[, r], nrow = nrow(restm))
                    colnames(detm) <- conditions[detected]
                    temp <- cbind(restm, detm)
                    result <- rbind(result, temp[, conditions])
                }
            }
            else {
                result <- rbind(result, rowi)
            }
        }
        if (methods::is(obj, "QCA_tt")) {
            mbase <- rev(c(1, cumprod(rev(noflevels))))[-1]
            diffwith <- NULL
            if (!observed) {
                diffwith <- obj$indexes
            }
            SBS <- setdiff(drop(result %*% mbase) + 1, diffwith)
        }
        else {
            SBS <- as.vector(which(apply(obj, 1, function(x) {
                x <- as.numeric(x)
                any(apply(result, 1, function(y) {
                    identical(x, as.numeric(y))
                }))
            })))
        }
    }
    if (any(is.element(type, c(0, 2)))) {
        pSA <- minimize(obj, include = "?")$SA
        pSAn <- minimize(nobj, include = "?")$SA
        SA1 <- pSA[[1]]
        if (length(pSA) > 1) {
            for (i in seq(2, length(pSA))) {
                SA1 <- rbind(SA1, pSA[[i]])
            }
        }
        SA2 <- pSAn[[1]]
        if (length(pSAn) > 1) {
            for (i in seq(2, length(pSAn))) {
                SA2 <- rbind(SA2, pSAn[[i]])
            }
        }
        CSA <- as.numeric(
            intersect(
                rownames(unique(SA1)),
                rownames(unique(SA2))
            )
        )
    }
    if (any(is.element(type, c(0, 3)))) {
        SSR <- as.numeric(
            intersect(
                rownames(obj$tt)[obj$tt$OUT == 1],
                rownames(nobj$tt)[nobj$tt$OUT == 1]
            )
        )
    }
    return(sort(unique(c(SBS, CSA, SSR))))
}
