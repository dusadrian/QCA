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

`solveChart` <- function(
    chart, row.dom = FALSE, all.sol = FALSE, depth = NULL, max.comb = 0,
    first.min = FALSE, ...
) {
    if (!is.logical(chart) && length(setdiff(chart, 0:1)) > 0) {
        admisc::stopError(
            "Use a logical, T/F matrix. See makeChart()'s output.", ... = ...
        )
    }
    dots <- list(...)
    if (is.element("min.dis", names(dots))) {
        if (is.logical(dots$min.dis)) {
            all.sol <- !dots$min.dis
        }
    }
    if (all.sol) {
        row.dom <- FALSE
    }
    row.numbers <- seq(nrow(chart))
    if (row.dom) {
        row.numbers <- rowDominance(chart)
        chart <- chart[row.numbers, ]
    }
    foundm <- findmin(chart, ... = ...) 
    if (foundm == 0) {
        admisc::stopError(
            "The PI chart cannot be solved.", ... = ...
        )
    }
    if (is.null(depth)) depth <- 0L
    output <- .Call(
        "C_solveChart",
        matrix(as.logical(chart),
        nrow = nrow(chart)),
        all.sol,
        as.integer(depth),
        as.integer(foundm),
        max.comb,
        first.min,
        PACKAGE = "QCA"
    )
    if (output[[2]]) {
        warning(
            simpleWarning(
                "The PI chart is exceedingly complex, solution(s) not guaranteed to be exhaustive.\n\n"
            )
        )
    }
    output <- output[[1]]
    output[output == 0] <- NA
    output <- matrix(as.integer(row.numbers[output]), nrow = nrow(output))
    output[is.na(output)] <- 0L
    return(output)
}
