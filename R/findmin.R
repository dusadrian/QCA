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

`findmin` <- function(chart, ...) {
    dots <- list(...)
    verbose <- isTRUE(dots$verbose)
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
        gurobi <- !isFALSE(attr(chart, "gurobi")) &&
                !isFALSE(dots$gurobi) &&
                eval(parse(
                    text = "requireNamespace('gurobi', quietly = TRUE)"
                ))
        just_minima <- !isTRUE(attr(chart, "solind")) && !isTRUE(dots$solind)
        if (gurobi) {
            chart <- matrix(as.numeric(chart), nrow = nrow(chart))
            model <- list(
                A = chart,
                obj = rep(1, ncol(chart)),
                modelsense = "min",
                rhs = rep(1, nrow(chart)),
                sense = rep(">=", nrow(chart)),
                vtype = "B"
            )
            params <- list(
                OutputFlag = verbose * 1,
                LogToConsole = verbose * 1
            )
            tc <- admisc::tryCatchWEM(
                solution <- eval(parse(text = "gurobi::gurobi(model, params)"))$x
            )
            if (!is.null(tc$error)) {
                gurobi <- FALSE
                if (cpi) {
                    message(
                        sprintf(
                            "%s, using lpSolve instead.",
                            ifelse(
                                grepl("license", tc$error),
                                "No valid Gurobi license found",
                                "Gurobi threw an error"
                            )
                        )
                    )
                }
            }
            gurobi <- is.null(tc)
        }
        if (!gurobi) {
            solution <- lpSolve::lp(
                "min",
                rep(1, ncol(chart)),
                chart,
                ">=",
                1,
                int.vec = seq(nrow(chart)),
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
