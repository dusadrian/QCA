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
