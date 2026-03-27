# Copyright (c) 2016 - 2026, Adrian Dusa
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

`findmin` <- function(chart, type = c("exact", "lagrangian"), ...) {
    dots <- list(...)
    verbose <- isTRUE(dots$verbose)
    if (missing(type) || is.null(type)) {
        type <- attr(chart, "type")
    }
    if (is.null(type)) {
        type <- "exact"
    }
    type <- match.arg(type, c("exact", "lagrangian"))
    if (!methods::is(chart, "QCA_pic")) {
        if (!is.matrix(chart) | (!is.logical(chart) & length(setdiff(chart, 0:1)) > 0)) {
            admisc::stopError(
                "Use a logical, TRUE/FALSE matrix. See makeChart()'s output.",
                ... = ...
            )
        }
    }
    cpi <- !is.null(attr(chart, "C_PI"))
    if (!cpi) {
        chart <- t(chart)
    }
    if (all(colSums(chart) > 0)) {
        just_minima <- !isTRUE(attr(chart, "solind")) && !isTRUE(dots$solind)
        gurobi <- !isFALSE(attr(chart, "gurobi")) && !isFALSE(dots$gurobi)
        solution <- NULL
        if (identical(type, "lagrangian")) {
            solution <- .Call(
                "C_findminLagrangian",
                matrix(as.logical(chart), nrow = nrow(chart)),
                PACKAGE = "QCA"
            )
        } else if (gurobi) {
            native_gurobi_available <- getOption("native.gurobi.available", NULL)
            if (is.null(native_gurobi_available)) {
                native_gurobi_available <- isTRUE(.Call("C_gurobiRuntimeAvailable", PACKAGE = "QCA"))
                options(native.gurobi.available = native_gurobi_available)
            }
            if (isTRUE(native_gurobi_available)) {
                solution <- .Call(
                    "C_findminExact",
                    matrix(as.logical(chart), nrow = nrow(chart)),
                    PACKAGE = "QCA"
                )
            } else if (verbose) {
                message("Gurobi not available, falling back to lpSolve.")
            }
        }
        if (is.null(solution)) {
            solution <- lpSolve::lp(
                direction = "min",
                objective.in = rep(1, ncol(chart)),
                const.mat = chart,
                const.dir = rep(">=", nrow(chart)),
                const.rhs = rep(1, nrow(chart)),
                all.bin = TRUE
            )$solution
        }
        if (just_minima | any(solution > 0 & solution < 1)) {
            solution <- as.integer(ceiling(sum(solution)))
        }
    }
    else {
        solution <- 0L
    }
    class(solution) <- c(class(solution), "QCA_findmin")
    return(solution)
}
