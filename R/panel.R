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

`as.panel` <- function(
    x, row.names, ...
) {
    if (!missing(row.names)) {
        if (is.character(row.names)) {
            if (length(row.names) == 1L) {
                rowvar <- (1L:ncol(x))[match(colnames(x), row.names, 0L) == 1L]
                row.names <- x[[rowvar]]
                x <- x[-rowvar]
            }
        }
        else if (is.numeric(row.names) && length(row.names) == 1L) {
            rowvar <- row.names
            row.names <- x[[rowvar]]
            x <- x[-rowvar]
        }
        else {
            admisc::stopError(
                "invalid 'row.names' specification.", ... = ...
            )
        }
        if (is.object(row.names) || !(is.integer(row.names))) {
            row.names <- as.character(row.names)
        }
        if (anyNA(row.names)) {
            admisc::stopError(
                "missing values in 'row.names' are not allowed.", ... = ...
            )
        }
        attr(x, "row.names") <- row.names
    }
    if (!is.data.frame(x)) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    structure(x, class = c("panel", "data.frame"))
}
`[.panel` <- function(x, i, j, ...) {
    funargs <- unlist(lapply(match.call(), deparse)[-1])
    drop <- TRUE
    if (any(names(funargs) == "drop")) {
        drop <- as.logical(funargs["drop"])
    }
    rownms <- row.names(x)
    class(x) <- "data.frame"
    classes <- lapply(x, class)
    clevels <- lapply(x, levels)
    cordered <- lapply(x, is.ordered)
    x <- eval(
        parse(
            text = sprintf(
                "x[%s, %s, drop = %s]",
                funargs["i"],
                funargs["j"],
                drop
            )
        )
    )
    if (!is.null(dim(x))) {
        x <- as.matrix(x)
        rownms <- eval(
            parse(
                text = sprintf(
                    "rownms[%s]",
                    funargs["i"]
                )
            )
        )
        row.names(x) <- rownms
        x <- rebuild(
            as.data.frame(x, stringsAsFactors = FALSE),
            classes[colnames(x)],
            clevels[colnames(x)],
            cordered[colnames(x)]
        )
        class(x) <- c("QCA_panel", "data.frame")
    }
    return(x)
}
`row.names<-.panel` <- function(x, value) {
    classes <- lapply(x, class)
    clevels <- lapply(x, levels)
    cordered <- lapply(x, is.ordered)
    x <- as.matrix(x)
    admisc::setRownames(x, value)
    x <- rebuild(
        as.data.frame(x, stringsAsFactors = FALSE),
        classes,
        clevels,
        cordered
    )
    class(x) <- c("QCA_panel", "data.frame")
    return(x)
}
`rebuild` <- function(x, classes, clevels, cordered) {
    for (i in seq(ncol(x))) {
        x[[i]] <- if (is.element("factor", classes[[i]])){
            factor(x[[i]], levels = clevels[[i]], ordered = cordered[[i]])
        }
        else if (is.element("Date", classes[[i]])) {
            as.Date(x[[i]])
        }
        else if (is.element("POSIXct", classes[[i]])) {
            as.POSIXct(x[[i]])
        }
        else {
            methods::as(x[[i]], classes[[i]])
        }
    }
    return(x)
}
