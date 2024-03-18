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

`modelFit` <- function(
    model, theory = "", select = NULL, ...
) {
    if (!(methods::is(model, "QCA_min") | methods::is(model, "admisc_deMorgan"))) {
        admisc::stopError(
            "The model should be a minimization object or its negation.",
            ... = ...
        )
    }
    theory <- admisc::recreate(substitute(theory))
    if (is.character(theory)) {
        if (length(theory) != 1) {
            admisc::stopError(
                "Theory should be a single SOP expression.",
                ... = ...
            )
        }
    }
    else {
        admisc::stopError(
            "Theory should be a SOP expression.",
            ... = ...
        )
    }
    noflevels <- model$tt$noflevels
    snames <- model$tt$options$conditions
    if (model$tt$options$use.letters) {
        snames <- LETTERS[seq(length(snames))]
    }
    pims <- model$pims
    if (is.element("i.sol", names(model))) {
        pims <- lapply(model$i.sol, function(x) x$pims)
        names(pims) <- NULL
        pims <- do.call("cbind", pims)
        solutions <- lapply(model$i.sol, function(x) x$solution)
    }
    else {
        solutions <- list(model$solution)
        if (identical(model$options$include, "")) {
            if (length(solutions[[1]][[1]]) > 2) {
                message("Warning: the negation of the conservative solution is potentially computer intensive.")
            }
        }
    }
    if (!is.null(select)) {
        if (!is.atomic(select) || !(is.numeric(select) | is.character(select))) {
            admisc::stopError(
                "Argument 'select' should be a numerical or character vector.",
                ... = ...
            )
        }
        if (is.character(select)) {
            if (!all(is.element(select, names(solutions)))) {
                admisc::stopError(
                    "Component specified with 'select' not found.",
                    ... = ...
                )
            }
        } else if (is.numeric(select)) {
            if (any(select > length(solutions))) {
                admisc::stopError(
                    "Numbers in 'select' greater than number of existing solutions.",
                    ... = ...
                )
            }
        }
        solutions <- solutions[select]
    }
    if (max(unlist(lapply(solutions, function(x) sapply(x, length)))) > 4) {
        message("Warning: the negation of the such model(s) is potentially computer intensive.")
    }
    models <- unlist(lapply(solutions, function(x) unlist(lapply(x, paste, collapse = " + "))))
    slengths <- unlist(lapply(solutions, length))
    if (is.null(names(solutions))) {
        names(models) <- "M"
        if (slengths > 1) {
            names(models) <- paste("M", seq(slengths), sep = "")
        }
    }
    else {
        mnum <- unlist(lapply(slengths, function(x) {
            mnum <- ""
            if (x > 1) {
                mnum <- seq(x)
            }
            paste("M", mnum, sep = "")
        }))
        names(models) <- paste(mnum, rep(names(solutions), slengths), sep = "-")
    }
    result <- intersections <- vector(mode = "list", length = length(models))
    arglist <- list(snames = snames, noflevels = noflevels)
    for (i in seq(length(models))) {
        expression <- models[i]
        cpims <- pims[, unlist(strsplit(expression, split = " \\+ ")), drop = FALSE]
        cpims$model <- admisc::compute(expression, data = model$tt$initial.data)
        testheory <- admisc::tryCatchWEM(
            cpims$theory <- admisc::compute(
                theory,
                data = model$tt$initial.data,
                enter = ""
            )
        )
        if (!is.null(testheory$error)) {
            if (grepl("multi-value", testheory$error)) {
                admisc::stopError(
                    "Theory expression should be specified in multi-value notation.",
                    ... = ...
                )
            }
        }
        intersections <- rep("", 4)
        negtheory <- negate(theory, snames = snames, noflevels = noflevels)[[1]][1]
        negexp <- negate(expression, snames = snames, noflevels = noflevels)[[1]][1]
        intersections[1] <- do.call(
            admisc::intersection,
            c(
                list(theory, expression),
                arglist
            )
        )
        intersections[2] <- do.call(
            admisc::intersection,
            c(
                list(
                    negtheory,
                    expression
                ),
                arglist
            )
        )
        intersections[3] <- do.call(
            admisc::intersection,
            c(
                list(
                    theory,
                    negexp
                ),
                arglist
            )
        )
        intersections[4] <- do.call(
            admisc::intersection,
            c(
                list(
                    negtheory,
                    negexp
                ),
                arglist
            )
        )
        intnms <- c("model*theory", "model*~theory", "~model*theory", "~model*~theory")
        for (nm in seq(4)) {
            int <- intersections[nm]
            if (int == "") {
                cpims[[intnms[nm]]] <- rep(0, nrow(model$tt$initial.data))
            }
            else {
                cpims[[intnms[nm]]] <- admisc::compute(int, data = model$tt$initial.data)
            }
        }
        intersections[intersections == ""] <- "-"
        names(intersections) <- intnms
        neg.out <- admisc::hastilde(model$tt$options$outcome)
        pofobj <- pof(
            cpims,
            as.numeric(model$tt$initial.data[, admisc::notilde(model$tt$options$outcome)]),
            relation = "sufficiency",
            neg.out = neg.out
        )
        pofobj$incl.cov <- pofobj$incl.cov[, 1:3]
        pofobj$incl.cov[is.na(pofobj$incl.cov[, 1]), 3] <- NA
        pofobj$modelfit <- list(
            model = expression,
            theory = theory,
            intersections = intersections
        )
        result[[i]] <- pofobj
    }
    if (length(result) == 1) {
        return(result[[1]])
    }
    return(structure(result, class = "QCA_modelFit"))
}
