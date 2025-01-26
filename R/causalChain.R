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

`causalChain` <- function(
    data, ordering = NULL, strict = FALSE, pi.cons = 0, pi.depth = 0,
    sol.cons = 0, sol.cov = 1, sol.depth = 0, ...
) {
    metacall <- match.call(expand.dots = TRUE)
    allargs <- as.list(metacall)[-1]
    allargs <- allargs[
        -which(
            is.element(
                c("data", "ordering", "strict"),
                names(allargs)
            )
        )
    ]
    sol.cons <- 1
    if (is.element("sol.cons", names(allargs))) {
        sol.cons <- allargs$sol.cons
    }
    sol.cons <- ifelse(sol.cons > 0 & sol.cons < 1, sol.cons, 1)
    pi.cons <- 1
    if (is.element("pi.cons", names(allargs))) {
        pi.cons <- allargs$pi.cons
    }
    if (
        any(c(pi.cons, sol.cons) < 1) & 
        !is.element("incl.cut", names(allargs))
    ) {
        allargs$incl.cut <- 0.5
    }
    verify.qca(data)
    noflevels  <- admisc::getInfo(data)$noflevels
    mv <- noflevels > 2
    names(noflevels) <- names(mv) <- colnames(data)
    if (methods::is(ordering, "character")) {
        ordering <- gsub("[[:space:]]", "", ordering)
        if (length(ordering) == 1) {
            ordering <- unlist(strsplit(ordering, split = "<"))
        }
        else {
            if (any(grepl("<", ordering))) {
                admisc::stopError(
                    "Causal ordering character \"<\" requires a single string.",
                    ... = ...
                )
            }
        }
        ordering <- lapply(ordering, admisc::splitstr)
    }
    if (length(allout <- unlist(ordering)) > 0) {
        if (length(setdiff(allout, colnames(data))) > 0) {
            admisc::stopError(
                "Some elements in the argument <ordering> not found in the data.",
                ... = ...
            )
        }
    }
    allargs <- c(list(input = data), allargs)
    allargs$causalChain <- TRUE
    checkpos <- function(x, arg) {
        pos <- pmatch(names(allargs), arg)
        return(pos[!is.na(pos)])
    }
    pos <- checkpos(allargs, "include")
    if (length(pos) == 0) {
        allargs$include <- "?"
    }
    pos <- checkpos(allargs, "all.sol")
    if (length(pos) == 0) {
        allargs$all.sol <- TRUE
    }
    pos <- checkpos(allargs, "SA")
    if (length(pos) == 0) {
        allargs$SA <- FALSE
    }
    minimizeit <- function(allargs) {
        return(
            tryCatch(
                do.call("minimize", allargs),
                error = function(e) NA
            )
        )
    }
    allargs$enter <- ""
    minimize.list <- list()
    if (length(ordering) > 0) {
        if (any(table(unlist(ordering)) > 1)) {
            admisc::stopError(
                "Same condition(s) in multiple ordering levels.",
                ... = ...
            )
        }
        allcols <- colnames(data)
        if (length(restcols <- setdiff(allcols, unlist(ordering))) > 0) {
            ordering <- c(list(restcols), ordering)
        }
        for (i in seq(length(ordering))) {
            nextcols <- ordering[[i]]
            if (i == 1) {
                if (!strict & length(nextcols) > 1) {
                    for (j in seq(length(nextcols))) {
                        allargs$input <- subset(data, select = nextcols)
                        if (mv[nextcols[j]]) {
                            uniqv <- sort(unique(data[, nextcols[j]]))
                            for (v in seq(noflevels[nextcols[j]] - 1)) {
                                if (is.element(v, uniqv)) {
                                    allargs$outcome <- sprintf("%s[%s]", nextcols[j], v)
                                    minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                                }
                            }
                        }
                        else {
                            allargs$outcome <- nextcols[j]
                            minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                        }
                    }
                }
            }
            else {
                restcols <- unlist(ordering[seq(i - 1)])
                for (j in seq(length(nextcols))) {
                    if (strict) {
                        allcols <- c(restcols, nextcols[j])
                    }
                    else {
                        allcols <- c(restcols, nextcols)
                    }
                    allcols <- allcols[order(match(allcols, colnames(data)))]
                    allargs$input <- subset(data, select = allcols)
                    if (mv[nextcols[j]]) {
                        uniqv <- sort(unique(data[, nextcols[j]]))
                        for (v in seq(noflevels[nextcols[j]] - 1)) {
                            if (is.element(v, uniqv)) {
                                allargs$outcome <- sprintf("%s[%s]", nextcols[j], v)
                                minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                            }
                        }
                    }
                    else {
                        allargs$outcome <- nextcols[j]
                        minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                    }
                }
            }
        }
    }
    else {
        for (x in colnames(data)) {
            if (mv[x]) {
                uniqv <- sort(unique(data[, x]))
                for (v in seq(noflevels[x] - 1)) {
                    if (is.element(v, uniqv)) {
                        allargs$outcome <- sprintf("%s[%s]", x, v)
                        minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                    }
                }
            }
            else {
                allargs$outcome <- x
                minimize.list[[x]] <- minimizeit(allargs)
            }
        }
    }
    attr(minimize.list, "call") <- metacall
    return(structure(minimize.list, class = "QCA_chain"))
}
