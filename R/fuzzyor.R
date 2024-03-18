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

`fuzzyor` <- function(
    ..., na.rm = FALSE
) {
    funargs <- unlist(lapply(
        lapply(match.call(), deparse)[-1],
        function(x) gsub(paste0("'|\"|[[:space:]|", "\u00a0", "]"), "", x)
    ))
    if (!is.na(rem <- match("na.rm", names(funargs)))) {
        funargs <- funargs[-rem]
    }
    dots <- vector(mode = "list", length = length(funargs))
    funargs <- gsub(paste(admisc::dashes(), collapse = "|"), "-", funargs)
    negated <- grepl("1-", funargs)
    funargs <- gsub("1-", "", funargs)
    tildenegated <- badnames <- logical(length(funargs))
    cols <- vector(mode = "list", length = length(funargs))
    for (i in seq(length(funargs))) {
        badnames[i] <- grepl("\\(|:", funargs[i])
        cols[[i]] <- admisc::getName(admisc::notilde(funargs[i]))
        tildenegated[i] <- admisc::tilde1st(funargs[i])
        funargs[i] <- admisc::notilde(funargs[i])
    }
    if (sum(badnames) > 0) {
        if (
            sum(badnames) > length(LETTERS) |
            any(is.element(unlist(cols), LETTERS))
        ) {
            cols[[badnames]] <- paste("X", seq(sum(badnames)), sep = "")
        }
        else {
            cols[[badnames]] <- LETTERS[seq(sum(badnames))]
        }
    }
    for (i in seq(length(funargs))) {
        tc <- tryCatchWEM(
            funi <- eval.parent(
                parse(text = funargs[i]), n = 1
            )
        )
        if (!is.null(tc$error)) {
            admisc::stopError(sprintf(
                "Object '%s' not found.", funargs[i]), ... = ...
            )
        }
        else {
            dots[[i]] <- drop(funi)
        }
    }
    if (!is.null(dots[[1]]) && is.atomic(dots[[1]])) {
        if (
            any(
                !unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Input vectors should be numeric or logical.", ... = ...
            )
        }
        dots <- as.data.frame(dots)
    }
    else if (is.matrix(dots[[1]])) {
        dots <- dots[[1]]
        if (is.null(colnames(dots))) {
            if (ncol(dots) > length(LETTERS)) {
                cols <- list(paste("X", seq(ncol(dots)), sep = ""))
            }
            else {
                cols <- list(LETTERS[seq(ncol(dots))])
            }
        }
        dots <- as.data.frame(dots)
        negated <- logical(ncol(dots))
        tildenegated <- logical(ncol(dots))
        if (
            !all(
                unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Input should be numeric or logical.", ... = ...
            )
        }
    }
    else if (is.data.frame(dots[[1]])) {
        dots <- dots[[1]]
        negated <- logical(ncol(dots))
        tildenegated <- logical(ncol(dots))
        cols <- list(colnames(dots))
        if (
            !all(
                unlist(
                    lapply(
                        dots,
                        function(x) is.numeric(x) | is.logical(x)
                    )
                )
            )
        ) {
            admisc::stopError(
                "Some columns are not numeric or logical.", ... = ...
            )
        }
    }
    else {
        admisc::stopError(
            "The input should be vectors, or a matrix or a dataframe.",
            ... = ...
        )
    }
    for (i in seq(length(cols))) {
        if (tildenegated[i]) {
            dots[[i]] <- 1 - dots[[i]]
        }
        if (negated[i]) {
            dots[[i]] <- 1 - dots[[i]]
        }
        if (negated[i] + tildenegated[i] == 1) {
            cols[[i]] <- paste("~", cols[[i]], sep = "")
        }
    }
    result <- apply(dots, 1, max, na.rm = na.rm)
    attr(result, "names") <- NULL
    attr(result, "name") <- paste(unlist(cols), collapse = " + ")
    class(result) <- c("numeric", "QCA_fuzzy")
    return(result)
}
