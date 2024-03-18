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

`Xplot` <- function(
    x, jitter = FALSE, at = pretty(x), ...
) {
    dots <- list(...)
    funargs <- unlist(lapply(match.call(), deparse)[-1])
    xname <- admisc::getName(funargs[["x"]])
    linex <- 1.75
    jitfactor <- 0.5
    jitamount <- 0.5
    cexpoints <- 1
    cexaxis <- 0.8
    pch <- 21
    bgpoints <- NA
    if (length(testarg <- which(names(dots) == "line")) > 0) {
        linex <- dots$line
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "factor")) > 0) {
        jitfactor <- dots$factor
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "amount")) > 0) {
        jitamount <- dots$amount
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "cex")) > 0) {
        cexpoints <- dots$cex
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "cex.axis")) > 0) {
        cexaxis <- dots$cex.axis
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "pch")) > 0) {
        pch <- dots$pch
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "bg")) > 0) {
        bgpoints <- dots$bg
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "xlab")) > 0) {
        xname <- dots$xlab
        dots <- dots[-testarg]
    }
    y <- rep(1, length(x))
    if (jitter) {
        y <- jitter(y, jitfactor, jitamount)
    }
    toplot <- list(as.name("plot"), x, y)
    toplot$type <- "n"
    if (!is.null(at)) {
        toplot$xlim <- range(at)
    }
    toplot$ylim <- c(0, 2)
    toplot$xlab <- ""
    toplot$ylab <- ""
    toplot$axes <- FALSE
    if (length(dots) > 0) {
        toplot <- c(toplot, dots)
    }
    par(mar = c(ifelse(xname == "", 2, 3), 0.3, 0, 0))
    suppressWarnings(eval(as.call(toplot)))
    axis(1, at = at, cex.axis = cexaxis)
    title(xlab = xname, cex.lab = cexaxis + 0.1, font.lab = 2, line = linex)
    plotpoints <- list(as.name("points"), x, y, pch = pch, cex = cexpoints, bg = bgpoints)
    suppressWarnings(eval(as.call(c(plotpoints, dots))))
}
