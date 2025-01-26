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

`getSolution` <- function(
    expressions, mv, collapse, inputt, row.dom, initial, all.sol, indata, curly, use.labels, ...
) {
    mtrx <- NULL
    sol.matrix <- NULL
    dots <- list(...)
    pi.cons <- if (is.element("pi.cons", names(dots))) dots$pi.cons else 0
    outcome <- if (is.element("outcome", names(dots))) dots$outcome else ""
    complex <- FALSE
    if (is.list(expressions)) {
        mtrx <- expressions[[2]]
        sol.matrix <- expressions[[3]]
        if (length(expressions) > 3) {
            complex <- expressions[[4]]
        }
        if (is.null(sol.matrix)) {
            admisc::stopError(
                "There are no solutions, given these constraints.",
                ... = ...
            )
        }
        expressions <- expressions[[1]]
        if (nrow(unique(expressions)) != nrow(expressions)) {
            expressions <- unique(expressions)
            mtrx <- NULL
            sol.matrix <- NULL
        }
    }
    if (
        nrow(expressions) == 1 &&
        identical(
            unique(
                as.vector(expressions)
            ),
            0L
        )
    ) {
        admisc::stopError(
            paste(
                "All truth table configurations are used, all conditions are minimized.",
                "       Please check the truth table.",
                sep = "\n"
            ),
            ... = ...
        )
    }
    if (FALSE) {
        if (!missing(indata)) {
            hastime <- logical(ncol(expressions))
            for (i in seq(ncol(expressions))) {
                if (any(is.element(indata[, i], c("-", "dc", "?")))) {
                    hastime[i] <- TRUE
                }
            }
            indata <- indata[, !hastime, drop = FALSE]
            expressions <- expressions[, !hastime, drop = FALSE]
            inputt <- inputt[, !hastime, drop = FALSE]
            relevant <- apply(expressions, 1, sum) > 0
            if (any(!relevant)) {
                sol.matrix <- NULL
                mtrx <- mtrx[relevant, , drop = FALSE]
                expressions <- expressions[relevant, , drop = FALSE]
            }
        }
    }
    PI <- admisc::writePrimeimp(
        impmat = expressions,
        mv = mv,
        collapse = collapse,
        curly = curly
    )
    rownames(expressions) <- PI
    if (pi.cons > 0 & outcome != "") {
        pofPI <- pof(
            paste(
                PI,
                collapse = " + "
            ),
            outcome = indata[, dots$outcome],
            data = indata,
            relation = "sufficiency",
            use.labels = use.labels
        )
        inclS <- pofPI$incl.cov[seq(length(PI)), 1]
        filterPI <- admisc::agteb(inclS, pi.cons)
        expressions <- expressions[filterPI, , drop = FALSE]
        PI <- PI[filterPI]
        mtrx <- makeChart(
            expressions,
            inputt,
            mv = mv,
            collapse = collapse,
            getSolution = TRUE,
            curly = curly
        )
        if (any(colSums(mtrx) == 0)) {
            admisc::stopError(
                "There are no solutions, given these constraints.",
                ... = ...
            )
        }
    }
    if (is.null(mtrx)) {
        mtrx <- makeChart(
            expressions,
            inputt,
            mv = mv,
            collapse = collapse,
            getSolution = TRUE,
            curly = curly
        )
    }
    else {
        rownames(mtrx) <- PI
    }
    notempty <- apply(mtrx, 1, any)
        expressions <- expressions[notempty, , drop = FALSE]
        mtrx <- mtrx[notempty, , drop = FALSE]
    admisc::setColnames(mtrx, initial)
    reduced <- list(expressions = expressions, mtrx = mtrx)
    if (nrow(mtrx) > 0) {
        if (row.dom & is.null(sol.matrix)) {
            reduced.rows <- rowDominance(mtrx)
            if (length(reduced.rows) > 0) {
                reduced$mtrx <- mtrx[reduced.rows, , drop = FALSE]
                reduced$expressions <- expressions[reduced.rows, , drop = FALSE]
            }
        }
        mtrx <- reduced$mtrx
        admisc::setColnames(mtrx, initial)
        if (is.null(sol.matrix)) {
            sol.matrix <- solveChart(mtrx, all.sol = all.sol, ... = ...)
        }
        tokeep <- sort(unique(as.vector(unique(sol.matrix))))
        all.PIs <- rownames(mtrx)[tokeep]
        solm <- matrix(as.integer(sol.matrix), nrow = nrow(sol.matrix))
        sol.matrix[sol.matrix == 0] <- NA
        sol.matrix <- matrix(rownames(mtrx)[sol.matrix], nrow = nrow(sol.matrix))
        reduced$expressions <- reduced$expressions[tokeep, , drop = FALSE]
        solution.list <- writeSolution(sol.matrix, mtrx)
    }
    else {
        all.PIs <- NA
        solution.list <- NA
        solm <- NA
    }
    return(
        list(
            expressions = expressions,
            mtrx = mtrx,
            reduced = reduced,
            all.PIs = all.PIs,
            solution.list = solution.list,
            sol.matrix = solm,
            complex = complex
        )
    )
}
