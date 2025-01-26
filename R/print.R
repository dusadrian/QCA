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

`print.QCA_aE` <- function(x, ...) {
    rownames(x) <- format(seq.nrow <- seq(nrow(x)))
    if (attr(x, "raw")) {
        x[x >= 0] <- paste("", x[x >= 0])
    }
    else {
        x[x < 0] <- " "
    }
    cat("\n")
    for (i in seq.nrow) {
        cat(paste(c(rownames(x)[i], x[i, ]), collapse = ifelse(attr(x, "raw"), "   ", "    ")), "\n")
    }
    cat("\n")
}
`print.QCA_chain` <- function(x, ...) {
    dots <- list(...)
    line.length <- getOption("width")
    if (any(names(x) == "via.web")) {
        line.length <- 10000
    }
    cat("\n")
    x <- lapply(x, function(x) {
        toreturn <- FALSE
        if (!identical(x, NA)) {
            sol.cons <- x$options$sol.cons
            sol.cov <- x$options$sol.cov
            if (!identical(x, NULL)) {
                outcome <- x$tt$options$outcome
                if (grepl(mvregexp, outcome)) {
                    if (x$options$neg.out) {
                        outcome <- paste("~", admisc::notilde(outcome), sep = "")
                    }
                }
                else {
                    if (x$options$neg.out) {
                        outcome <- paste("~", admisc::notilde(outcome), sep = "")
                    }
                }
                if (is.element("i.sol", names(x))) {
                    sufnec <- logical(length(x$i.sol))
                    for (i in seq(length(x$i.sol))) {
                        if (is.element("overall", names(x$i.sol[[i]]$IC))) {
                            sufnec[i] <- all(admisc::agteb(x$i.sol[[i]]$IC$overall$sol.incl.cov[c(1, 3)], c(sol.cons, sol.cov)))
                        }
                        else {
                            sufnec[i] <- all(admisc::agteb(x$i.sol[[i]]$IC$sol.incl.cov[c(1, 3)], c(sol.cons, sol.cov)))
                        }
                    }
                    if (any(sufnec)) {
                        for (i in seq(length(x$i.sol))) {
                            if (sufnec[i]) {
                                prettyNums <- formatC(seq(length(x$solution)), digits = nchar(length(x$solution)) - 1, flag = 0)
                                for (sol in seq(length(x$i.sol[[i]]$solution))) {
                                    preamble <- sprintf("M%s-%s:", i, sol)
                                    preamble <- paste(preamble, paste(rep(" ", 7 - nchar(preamble)), collapse = ""), sep = "")
                                    cat(preamble)
                                    cat(admisc::prettyString(x$i.sol[[i]]$solution[[sol]], line.length - 7, 7, "+", "<->", outcome), "\n")
                                }
                            }
                        }
                        toreturn <- TRUE
                    }
                }
                else {
                    sufnec <- logical(length(x$solution))
                    if (length(x$solution) == 1) {
                        sufnec <- all(admisc::agteb(x$IC$sol.incl.cov[c(1, 3)], c(sol.cons, sol.cov)))
                        if (is.na(sufnec)) sufnec <- FALSE
                        if (sufnec) {
                            cat(paste("M1: ", admisc::prettyString(x$solution[[1]], line.length - 4, 4, "+", "<->", outcome), "\n", sep = ""))
                            toreturn <- TRUE
                        }
                    }
                    else {
                        for (i in seq(length(x$solution))) {
                            sufnec[i] <- all(admisc::agteb(x$IC$individual[[i]]$sol.incl.cov[c(1, 3)], c(sol.cons, sol.cov)))
                            if (is.na(sufnec[i])) sufnec[i] <- FALSE
                        }
                        if (any(sufnec)) {
                            prettyNums <- formatC(seq(length(x$solution)), digits = nchar(length(x$solution)) - 1, flag = 0)
                            for (i in seq(length(x$solution))) {
                                if (sufnec[i]) {
                                    cat(paste("M", prettyNums[i], ": ", sep = ""))
                                    cat(sprintf("%s <-> %s\n", paste(x$solution[[i]], collapse = " + "), outcome))
                                }
                            }
                            toreturn <- TRUE
                        }
                    }
                }
            }
            if (toreturn) {
                cat("\n")
            }
        }
        return(toreturn)
    })
    if (all(!unlist(x))) {
        cat("There are no causal chains in this data.\n\n")
    }
}
`print.QCA_findmin` <- function(x, ...) {
    attr(x, "solution") <- NULL
    print(unclass(x))
}
`print.QCA_fuzzy` <- function(x, ...) {
    attr(x, "name") <- NULL
    print(unclass(x))
    cat("\n")
}
`print.QCA_modelFit` <- function(x, ...) {
    for (i in seq(length(x))) {
        print(x[[i]])
    }
}
`print.QCA_loopmin` <- function(x, ...) {
    cat("\n")
    for (i in seq(length(x))) {
        print.QCA_min(x[[i]], details = FALSE, mqca = TRUE)
    }
}
`print.QCA_panel` <- function(x, ...) {
    dots <- list(...)
    quote <- FALSE
    if (is.element("quote", names(dots))) {
        if (is.logical(dots$quote)) {
            quote <- dots$quote[1]
        }
    }
    right <- TRUE
    if (is.element("right", names(dots))) {
        if (is.logical(dots$right)) {
            right <- dots$right[1]
        }
    }
    n <- length(row.names(x))
    if (length(x) == 0L) {
        cat(sprintf(ngettext(n, "cluster data frame with 0 columns and %d row",
            "cluster data frame with 0 columns and %d rows"), n), "\n",
            sep = "")
    }
    else if (n == 0L) {
        print.default(names(x), quote = FALSE)
        cat(gettext("<0 rows> (or 0-length row.names)\n"))
    }
    else {
        x <- as.matrix(x)
        if (is.null(row.names(x))) {
            row.names(x) <- seq(nrow(x))
        }
        print(x, quote = quote, right = right)
    }
}
`print.QCA_pic` <- function(x, ...) {
    if (all(dim(x) > 0)) {
        rownms  <- rownames(x)
        colnms  <- colnames(x)
        x2 <- matrix(as.vector(x), nrow = nrow(x))
        rownames(x2) <- paste(rownms, "")
        colnames(x2) <- format(colnms, width = 2)
        x2[as.vector(x)]  <- "x"
        x2[!as.vector(x)] <- "-"
        x <- x2
    }
    cat("\n")
    print(admisc::prettyTable(x))
    cat("\n")
}
`print.QCA_pof` <- function(x, ...) {
    if (is.element("modelfit", names(x))) {
        cat("\n")
        if (names(x$modelfit$model) != "M") {
            cat(names(x$modelfit$model), "\n", sep = "")
        }
        cat("model:          ", x$modelfit$model, "\n", sep = "")
        cat("theory:         ", x$modelfit$theory, "\n", sep = "")
        if (any(grepl("~", names(x$modelfit$intersections)))) {
            cat(names(x$modelfit$intersections)[1], ":   ", x$modelfit$intersections[1], "\n", sep = "")
            cat(names(x$modelfit$intersections)[2], ":  ", x$modelfit$intersections[2], "\n", sep = "")
            cat(names(x$modelfit$intersections)[3], ":  ", x$modelfit$intersections[3], "\n", sep = "")
            cat(names(x$modelfit$intersections)[4], ": ", x$modelfit$intersections[4], "\n", sep = "")
        }
        else {
            for (int in seq(4)) {
                cat(names(x$modelfit$intersections)[int], ": ", x$modelfit$intersections[int], "\n", sep = "")
            }
        }
    }
    essential.PIs <- NULL
    if (is.element("essential", names(x))) {
        essential.PIs <- x$essential
    }
    essentials <- length(essential.PIs) > 0
    overall <- FALSE
    if (is.element("overall", names(x))) {
        overall <- TRUE
    }
    cases.column <- sol.exists <- FALSE
    valid.covU <- TRUE
    dots <- list(...)
    max.nchar.cases <- 0
    line.length <- getOption("width")
    if (is.element("line.length", names(dots))) {
        line.length <- dots$line.length
    }
    if (!is.element("show.cases", names(x$options))) {
        x$options$show.cases <- FALSE
    }
    if (is.element("show.cases", names(dots))) {
        if (is.logical(dots$show.cases)) {
            x$options$show.cases <- dots$show.cases
        }
    }
    if (overall) {
        incl.cov <- x$overall$incl.cov
        nrow.incl.cov <- nrow(incl.cov)
        nchar.nrow <- nchar(nrow.incl.cov)
        ind.len <- length(x$individual)
        if (essentials) {
            essential.PIs.rows <- is.element(rownames(incl.cov), essential.PIs)
        }
        if (x$options$show.cases) {
            max.nchar.cases <- max(nchar(encodeString(incl.cov$cases)))
            cases.column <- TRUE
            incl.cov.cases <- incl.cov$cases
            incl.cov$cases <- NULL
            if (essentials) {
                incl.cov.e.cases <- incl.cov.cases[essential.PIs.rows]
                incl.cov.cases <- incl.cov.cases[!essential.PIs.rows]
            }
        }
        else {
            incl.cov$cases <- NULL
        }
        prettyNums <- format(seq(nrow.incl.cov))
        for (i in seq(ncol(incl.cov))) {
            NAs <- is.na(incl.cov[, i])
            incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits = 3, format = "f")
            incl.cov[NAs, i] <- "  -  "
        }
        colnames(incl.cov) <- format(colnames(incl.cov), justify = "centre")
        if (essentials) {
            which.essential <- seq(length(which(essential.PIs.rows)))
            prettyNums.e <- prettyNums[which.essential]
            prettyNums <- prettyNums[-which.essential]
            incl.cov.e <- incl.cov[essential.PIs.rows, , drop = FALSE]
            incl.cov <- incl.cov[!essential.PIs.rows, , drop = FALSE]
            for (i in seq(ind.len)) {
                essentialPIs <- is.element(
                    rownames(x$individual[[i]]$incl.cov),
                    essential.PIs
                )
                unique.coverages <- formatC(
                    x$individual[[i]]$incl.cov$covU[essentialPIs],
                    digits = 3,
                    format = "f"
                )
                incl.cov.e <- cbind(incl.cov.e, S = unique.coverages, stringsAsFactors = FALSE)
                x$individual[[i]]$incl.cov <- x$individual[[i]]$incl.cov[!essentialPIs, ]
                if (x$options$use.labels && length(x$categories) > 0) {
                    rownames(x$individual[[i]]$incl.cov) <- replaceCategories(
                        rownames(x$individual[[i]]$incl.cov),
                        categories = x$categories
                    )
                }
            }
        }
        if (nrow(incl.cov) > 0) {
            for (i in seq(ind.len)) {
                incl.cov <- cbind(incl.cov, "     ", stringsAsFactors = FALSE)
                colnames(incl.cov)[ncol(incl.cov)] <- format(
                    ifelse(
                        ind.len < line.length,
                        paste("(M", i, ")", sep = ""),
                        paste("M", i, sep = "")
                    ),
                    width = 5
                )
                if (length(x$individual[[i]]$incl.cov$covU) > 0) {
                    incl.cov[rownames(x$individual[[i]]$incl.cov), ncol(incl.cov)] <- formatC(
                        x$individual[[i]]$incl.cov$covU,
                        digits = 3,
                        format = "f"
                    )
                }
            }
        }
        sol.incl.cov <- matrix(
            unlist(lapply(x$individual, "[", "sol.incl.cov")),
            nrow = length(x$individual),
            ncol = 3,
            byrow = TRUE
        )
        rownames(sol.incl.cov) <- paste("M", seq(length(x$individual)), sep = "")
        sol.exists <- TRUE
    }
    else {
        incl.cov <- x$incl.cov
        nrow.incl.cov <- nrow(incl.cov)
        nchar.nrow <- nchar(nrow.incl.cov)
        prettyNums <- format(seq(nrow.incl.cov))
        incl.cov[incl.cov == "  NA"] <- "     "
        colnames(incl.cov) <- format(colnames(incl.cov), justify = "centre")
        if (x$options$show.cases) {
            max.nchar.cases <- max(5, max(nchar(encodeString(incl.cov$cases))))
            cases.column <- TRUE
            incl.cov.cases <- incl.cov$cases
            incl.cov$cases <- NULL
        }
        incl.cov$cases <- NULL
        for (i in seq(ncol(incl.cov))) {
            NAs <- is.na(incl.cov[, i])
            incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits = 3, format = "f")
            incl.cov[NAs, i] <- "  -  "
        }
        for (i in seq(ncol(x$optionals))) {
            NAs <- is.na(x$optionals)
            x$optionals[!NAs, i] <- formatC(x$optionals[!NAs, i], digits = 3, format = "f")
            x$optionals[NAs, i] <- "  -  "
        }
        if (is.element("sol.incl.cov", names(x))) {
            sol.incl.cov <- as.matrix(x$sol.incl.cov)
            rownames(sol.incl.cov) <- "M1"
            sol.exists <- TRUE
        }
    }
    if (is.null(rownames(incl.cov))) {
        rownames(incl.cov) <- rep("  ", nrow(incl.cov))
    }
    if (x$options$use.labels && length(x$categories) > 0) {
        rownames(incl.cov) <- replaceCategories(
            rownames(incl.cov),
            categories = x$categories
        )
    }
    nchar.rownames <- max(nchar(rownames(incl.cov)))
    if (nchar.rownames == 1) {
        nchar.rownames <- 2
    }
    if (essentials) {
        nchar.rownames <- max(
            nchar.rownames,
            max(nchar(rownames(incl.cov.e)))
        )
        rownames(incl.cov.e) <- sprintf(
            paste("% ", nchar.rownames, "s", sep = ""),
            rownames(incl.cov.e)
        )
    }
    rownames(incl.cov) <- sprintf(
        paste("% ", nchar.rownames, "s", sep = ""),
        rownames(incl.cov)
    )
    if (sol.exists) {
        rownames(sol.incl.cov) <- sprintf(
            paste("% ", nchar.rownames, "s", sep = ""),
            rownames(sol.incl.cov)
        )
        NAs <- is.na(sol.incl.cov)
        sol.incl.cov[!NAs] <- formatC(sol.incl.cov[!NAs], digits = 3, format = "f")
        sol.incl.cov[NAs] <- ""
    }
    incl.cov[incl.cov == "  NA"] <- "  -  "
    max.chars <- 1
    if (is.element(x$options$relation, c("sufficiency", "suf"))) {
        if (
            ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov)))) &
            is.null(x$options$add)
        ) {
            first.printed.row <- paste(
                c(
                    rep(" ", nchar.rownames + nchar.nrow + 25),
                    rep("-", 7 * (ncol(incl.cov) - (2 + valid.covU)) - 2)
                ),
                collapse = ""
            )
            max.chars <- nchar(first.printed.row)
        }
    }
    if (max.chars < line.length) {
        if (cases.column) {
            max.chars <- max.nchar.cases
        }
        sep.row <- paste(
            rep(
                "-",
                nchar.rownames +
                7 * ncol(incl.cov) +
                ifelse(cases.column, max.nchar.cases, 0) +
                nchar.nrow
            ),
            collapse = ""
        )
        if (nchar(sep.row) < line.length) {
            if (
                ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov)))) &
                length(intersect(colnames(incl.cov), "pval1")) == 0 &
                is.null(x$options$add)
            ) {
                cat(first.printed.row, "\n")
            }
            else {
                cat("\n")
            }
            colstoprint <- colnames(incl.cov)
            colnames.row <- cat(
                paste(
                    c(
                        paste(
                            rep(" ", nchar.rownames + nchar.nrow + 2),
                            collapse = ""
                        ),
                        format(colstoprint)
                    ),
                    collapse = "  "
                )
            )
            cat(
                paste(
                    colnames.row,
                    ifelse(cases.column, "  cases", ""),
                    sep = ""
                ),
                "\n"
            )
            sep.row <- paste(
                rep(
                    "-",
                    nchar.rownames +
                    7 * ncol(incl.cov) +
                    ifelse(cases.column, max.nchar.cases + 2, 0) +
                    nchar.nrow +
                    2
                ),
                collapse = ""
            )
            cat(sep.row, "\n")
            if (essentials) {
                for (i in seq(nrow(incl.cov.e))) {
                    i.row <- paste(
                        prettyNums.e[i],
                        paste(
                            c(rownames(incl.cov.e)[i], incl.cov.e[i, ]),
                            collapse = "  "
                        ),
                        sep = "  "
                    )
                    if (cases.column) {
                        i.row <- paste(i.row, incl.cov.e.cases[i], sep = "  ")
                    }
                    cat(i.row, "\n")
                }
                cat(sep.row, "\n")
            }
            for (i in seq(nrow(incl.cov))) {
                rowtoprint <- c(rownames(incl.cov)[i], incl.cov[i, ])
                i.row <- paste(
                    prettyNums[i],
                    paste(rowtoprint, collapse = "  "),
                    sep = "  "
                )
                if (cases.column) {
                    i.row <- paste(i.row, incl.cov.cases[i], sep = "  ")
                }
                cat(i.row, "\n")
            }
            cat(sep.row, "\n")
            if (sol.exists) {
                for (i in seq(nrow(sol.incl.cov))) {
                    cat(
                        paste(
                            paste(
                                rep(" ", nchar.nrow),
                                collapse = ""
                            ),
                            paste(
                                c(
                                    rownames(sol.incl.cov)[i],
                                    sol.incl.cov[i, ]
                                ),
                                collapse = "  "
                            ),
                            sep = "  "
                        ),
                        "\n"
                    )
                }
            }
            cat("\n")
        }
        else {
            if (
                ncol(incl.cov) > (3 + any(grepl("PRI|RoN", colnames(incl.cov)))) &
                is.null(x$options$add)
            ) {
                cat(first.printed.row, "\n")
            } else {
                cat("\n")
            }
            cat(
                paste(
                    c(
                        paste(
                            rep(" ", nchar.rownames + nchar.nrow + 2),
                            collapse = ""
                        ),
                        colnames(incl.cov)
                    ),
                    collapse = "  "
                ),
                "\n"
            )
            sep.row <- paste(
                rep("-", nchar.rownames + 7 * ncol(incl.cov) + nchar.nrow + 2),
                collapse = ""
            )
            cat(sep.row, "\n")
            if (essentials) {
                for (i in seq(nrow(incl.cov.e))) {
                    cat(
                        paste(
                            prettyNums.e[i],
                            paste(
                                c(
                                    rownames(incl.cov.e)[i],
                                    incl.cov.e[i, ]
                                ),
                                collapse = "  "
                            ),
                            "\n"
                        ),
                        sep = "  "
                    )
                }
                cat(sep.row, "\n")
            }
            for (i in seq(nrow(incl.cov))) {
                cat(
                    paste(
                        prettyNums[i],
                        paste(
                            c(
                                rownames(incl.cov)[i],
                                incl.cov[i, ]
                            ),
                            collapse = "  "
                        ),
                        sep = "  "
                    ),
                    "\n"
                )
            }
            cat(sep.row, "\n")
            if (sol.exists) {
                for (i in seq(nrow(sol.incl.cov))) {
                    cat(
                        paste(
                            paste(
                                rep(" ", nchar.nrow),
                                collapse = ""
                            ),
                            paste(
                                c(
                                    rownames(sol.incl.cov)[i],
                                    sol.incl.cov[i, ]
                                ),
                                collapse = "  "
                            ),
                            sep = "  "
                        ),
                        "\n"
                    )
                }
            }
            if (cases.column) {
                cat(
                    "\n",
                    paste(
                        paste(
                            rep(" ", nchar.rownames + nchar.nrow + 2),
                            collapse = ""
                        ),
                        "cases"
                    ),
                    "\n"
                )
                cat(
                    paste(
                        rep("-", nchar.rownames + 7 + nchar.nrow + 2),
                        collapse = ""
                    ),
                    "\n"
                )
                if (essentials) {
                    for (i in seq(nrow(incl.cov.e))) {
                        cat(
                            paste(
                                prettyNums.e[i],
                                paste(
                                    rownames(incl.cov.e)[i],
                                    " "
                                ),
                                sep = "  "
                            )
                        )
                        cases <- unlist(
                            strsplit(
                                incl.cov.e.cases[i],
                                split = "; ",
                                useBytes = TRUE
                            )
                        )
                        cat(
                            admisc::prettyString(
                                cases,
                                getOption("width") - nchar.rownames - nchar.nrow - 4,
                                nchar.rownames + nchar.nrow + 4,
                                ";",
                                cases = TRUE
                            )
                        )
                        cat("\n")
                    }
                    cat(
                        paste(
                            rep("-", nchar.rownames + nchar.nrow + 9),
                            collapse = ""
                        ),
                        "\n"
                    )
                }
                for (i in seq(nrow(incl.cov))) {
                    cat(
                        paste(
                            prettyNums[i],
                            paste(
                                rownames(incl.cov)[i],
                                " "
                            ),
                            sep = "  "
                        )
                    )
                    cases <- unlist(
                        strsplit(
                            incl.cov.cases[i],
                            split = "; ",
                            useBytes = TRUE
                        )
                    )
                    cat(
                        admisc::prettyString(
                            cases,
                            getOption("width") - nchar.rownames - nchar.nrow - 4,
                            nchar.rownames + nchar.nrow + 4,
                            ";",
                            cases = TRUE
                        )
                    )
                    cat("\n")
                }
                cat(
                    paste(
                        rep("-", nchar.rownames + nchar.nrow + 9),
                        collapse = ""
                    ),
                    "\n\n"
                )
            }
        }
    }
    else {
        ncols <- floor((line.length - nchar.rownames)/7)
        if (ncols < 0) {
            admisc::stopError(
                "Too complex to print. Try using single letters for all conditions."
            )
        }
        chunks <- ceiling(ncol(incl.cov)/ncols)
        colsplits <- seq(1, ncol(incl.cov), by=ncols)
        for (chunk in seq(chunks)) {
            sel.cols <- seq(
                colsplits[chunk],
                ifelse(
                    chunk == chunks,
                    ncol(incl.cov),
                    colsplits[chunk + 1] - 1
                )
            )
            incl.cov.temp <- incl.cov[, sel.cols, drop = FALSE]
            if (essentials) {
                incl.cov.e.temp <- incl.cov.e[, sel.cols, drop = FALSE]
            }
            if (chunk < chunks) {
                if (ncols > 3) { 
                    cat(
                        paste(
                            c(
                                "\n",
                                rep(
                                    ifelse(chunk == 1, " ", "-"),
                                    nchar.rownames + nchar.nrow + 18
                                ),
                                rep("-", 7 * (ncols - 2) - 2),
                                "\n"
                            ),
                            collapse = ""
                        )
                    )
                }
                cat(
                    paste(
                        c(
                            paste(
                                rep(" ", nchar.rownames + nchar.nrow + 2),
                                collapse = ""
                            ),
                            colnames(incl.cov.temp)
                        ),
                        collapse = "  "
                    ),
                    "\n"
                )
                sep.row <- paste(
                    rep(
                        "-",
                        nchar.rownames + 7 * ncol(incl.cov.temp) + nchar.nrow + 2
                    ),
                    collapse = ""
                )
                cat(sep.row, "\n")
                if (essentials) {
                    for (i in seq(nrow(incl.cov.e.temp))) {
                        cat(
                            paste(
                                prettyNums.e[i],
                                paste(
                                    c(
                                        rownames(incl.cov.e.temp)[i],
                                        incl.cov.e.temp[i, ]
                                    ),
                                    collapse = "  "
                                ),
                                sep = "  "
                            ),
                            "\n"
                        )
                    }
                    cat(sep.row, "\n")
                }
                for (i in seq(nrow(incl.cov.temp))) {
                    cat(
                        paste(
                            prettyNums[i],
                            paste(
                                c(
                                    rownames(incl.cov.temp)[i],
                                    incl.cov.temp[i, ]
                                ),
                                collapse = "  "
                            ),
                            sep = "  "
                        ),
                        "\n"
                    )
                }
                cat(sep.row, "\n")
                if (chunk == 1 & sol.exists) {
                    for (i in seq(nrow(sol.incl.cov))) {
                        cat(
                            paste(
                                paste(
                                    rep(" ", nchar.nrow),
                                    collapse = ""
                                ),
                                paste(
                                    c(
                                        rownames(sol.incl.cov)[i],
                                        sol.incl.cov[i, ]
                                    ),
                                    collapse = "  "
                                ),
                                sep = "  "
                            ),
                            "\n"
                        )
                    }
                }
                cat("\n")
            }
            else {
                max.chars <- nchar.rownames + 7 * ncol(incl.cov.temp) + nchar.nrow + 2
                sep.row <- paste(c(rep("-", max.chars)), collapse = "")
                if (cases.column) {
                    max.chars <- max.chars + max.nchar.cases
                }
                if (max.chars < line.length) {
                    cat(sep.row, "\n")
                    sep.row <- paste(
                        sep.row,
                        ifelse(
                            cases.column,
                            paste(
                                rep("-", max.nchar.cases + 2),
                                collapse = ""
                            ),
                            ""
                        ),
                        sep = ""
                    )
                    colnames.row <- paste(
                        c(
                            paste(
                                rep(" ", nchar.rownames + nchar.nrow + 2),
                                collapse = ""
                            ),
                            colnames(incl.cov.temp)
                        ),
                        collapse = "  "
                    )
                    cat(
                        paste(
                            colnames.row,
                            ifelse(cases.column, "  cases", ""),
                            sep = ""
                        ),
                        "\n"
                    )
                    cat(sep.row, "\n")
                    if (essentials) {
                        for (i in seq(nrow(incl.cov.e.temp))) {
                            i.row <- paste(
                                prettyNums.e[i],
                                paste(
                                    c(
                                        rownames(incl.cov.e.temp)[i],
                                        incl.cov.e.temp[i, ]
                                    ),
                                    collapse = "  "
                                ),
                                sep = "  "
                            )
                            if (cases.column) {
                                i.row <- paste(
                                    i.row,
                                    incl.cov.e.cases[i],
                                    sep = "  "
                                )
                            }
                            cat(i.row, "\n")
                        }
                        cat(sep.row, "\n")
                    }
                    for (i in seq(nrow(incl.cov.temp))) {
                        i.row <- paste(
                            prettyNums[i],
                            paste(
                                c(
                                    rownames(incl.cov.temp)[i],
                                    incl.cov.temp[i, ]
                                ),
                                collapse = "  "
                            ),
                            sep = "  "
                        )
                        if (cases.column) {
                            i.row <- paste(
                                i.row,
                                incl.cov.cases[i],
                                sep = "  "
                            )
                        }
                        cat(i.row, "\n")
                    }
                    cat(sep.row, "\n")
                }
                else {
                    cat(sep.row, "\n")
                    cat(
                        paste(
                            c(
                                paste(
                                    rep(" ", nchar.rownames + nchar.nrow + 2),
                                    collapse = ""
                                ),
                                colnames(incl.cov.temp)
                            ),
                            collapse = "  "
                        ),
                        "\n"
                    )
                    cat(sep.row, "\n")
                    if (essentials) {
                        for (i in seq(nrow(incl.cov.e.temp))) {
                            i.row <- paste(
                                prettyNums.e[i],
                                paste(
                                    c(
                                        rownames(incl.cov.e.temp)[i],
                                        incl.cov.e.temp[i, ]
                                    ),
                                    collapse = "  "
                                ),
                                sep = "  "
                            )
                            cat(i.row, "\n")
                        }
                        cat(sep.row, "\n")
                    }
                    for (i in seq(nrow(incl.cov.temp))) {
                        cat(
                            paste(
                                prettyNums[i],
                                paste(
                                    c(
                                        rownames(incl.cov.temp)[i],
                                        incl.cov.temp[i, ]
                                    ),
                                    collapse = "  "
                                ),
                                sep = "  "
                            ),
                            "\n"
                        )
                    }
                    cat(sep.row, "\n")
                    if (cases.column) {
                        cat("\n", paste(paste(rep(" ", nchar.rownames + nchar.nrow + 2), collapse = ""), "cases"), "\n")
                        sep.row <- paste(rep("-", nchar.rownames + nchar.nrow + 9), collapse = "")
                        cat(sep.row, "\n")
                        if (essentials) {
                            for (i in seq(nrow(incl.cov.e.temp))) {
                                cat(paste(prettyNums[i], paste(rownames(incl.cov.e.temp)[i], " "), sep = "  "))
                                cases <- unlist(strsplit(incl.cov.e.cases[i], split = "; ", useBytes = TRUE))
                                cat(
                                    admisc::prettyString(
                                        cases,
                                        getOption("width") - nchar.rownames - 2,
                                        nchar.rownames + 2,
                                        ";",
                                        cases = TRUE
                                    )
                                )
                                cat("\n")
                            }
                            cat(sep.row, "\n")
                        }
                        for (i in seq(nrow(incl.cov.temp))) {
                            cat(
                                paste(
                                    prettyNums[i],
                                    paste(
                                        rownames(incl.cov.temp)[i],
                                        " "
                                    ),
                                    sep = "  "
                                )
                            )
                            cases <- unlist(
                                strsplit(
                                    incl.cov.cases[i],
                                    split = "; ",
                                    useBytes = TRUE
                                )
                            )
                            cat(
                                admisc::prettyString(
                                    cases,
                                    getOption("width") - nchar.rownames - 2,
                                    nchar.rownames + 2,
                                    ";",
                                    cases = TRUE
                                )
                            )
                            cat("\n")
                        }
                        cat(sep.row, "\n")
                    }
                }
                cat("\n")
            }
        }
    }
}
`print.QCA_min` <- function(x, ...) {
    enter <- ifelse (
        is.element("enter", names(as.list(x$call))),
        as.list(x$call)$enter,
        "\n"
    )
    dots <- list(...)
    enter <- identical(enter, "\n")
    essential_replaced <- FALSE
    replace <-  length(x$tt$categories) > 0 &&
                (
                    isTRUE(dots$use.labels) ||
                    (
                        !isFALSE(dots$use.labels) && x$tt$options$use.labels
                    )
                )
    line.length <- getOption("width")
    if (any(names(x) == "via.web")) {
        line.length <- 10000
    }
    details <- x$options$details
    mqca <- FALSE
    if (is.element("mqca", names(dots))) {
        if (is.logical(dots$mqca)) {
            mqca <- dots$mqca
        }
    }
    sol.cons <- x$options$sol.cons
    sol.cov  <- x$options$sol.cov
    outcome <- x$tt$options$outcome
    if (x$options$curly) {
        outcome <- gsub("\\[", "{", gsub("\\]", "}", outcome))
    }
    if (grepl(mvregexp, outcome)) {
        if (x$options$neg.out) {
            outcome <- paste("~", admisc::notilde(outcome), sep = "")
        }
        if (any(x$options$explain != 1)) {
            outcome <- ""
        }
    }
    else {
        if (x$options$neg.out) {
            if (admisc::tilde1st(outcome)) {
                outcome <- admisc::notilde(outcome)
            }
            else {
                outcome <- paste("~", outcome, sep = "")
            }
        }
    }
    if (replace) {
        outcome <- replaceCategories(
            outcome,
            categories = x$tt$categories
        )
    }
    if (is.element("show.cases", names(dots))) {
        if (is.logical(dots$show.cases)) {
            x$options$show.cases <- dots$show.cases
        }
    }
    if (!x$options$show.cases) {
        if (is.element("cases", names(x$IC$incl.cov))) {
            x$IC$incl.cov$cases <- NULL
        }
    }
    if (is.element("details", names(dots))) {
        if (is.logical(dots$details)) {
            details <- dots$details
            x$options$details <- details
        }
    }
    if (!mqca & enter) {
        cat("\n")
    }
    if (is.element("i.sol", names(x))) {
        sufnec <- logical(length(x$i.sol))
        for (i in seq(length(x$i.sol))) {
            if (is.element("overall", names(x$i.sol[[i]]$IC))) {
                sufnec[i] <- all(admisc::agteb(x$i.sol[[i]]$IC$overall$sol.incl.cov[3], sol.cov))
            }
            else {
                sufnec[i] <- all(
                    admisc::agteb(
                        x$i.sol[[i]]$IC$sol.incl.cov[3],
                        sol.cov
                    )
                )
            }
        }
        sufnec.char <- rep("", length(sufnec))
        uniques <- unique(lapply(x$i.sol, function(x) x$solution))
        for (j in seq(length(uniques))) {
            indices <- unlist(lapply(x$i.sol, function(x) identical(uniques[[j]], x$solution)))
            isols <- names(indices)[indices]
            cat(paste(ifelse(j > 1, "\n", ""), "From ", paste(isols, collapse = ", "), ": ", sep = ""))
            if (enter) cat("\n")
            i <- which(names(x$i.sol) == isols[1])
            if (!mqca & enter) {
                cat("\n")
            }
            for (sol in seq(length(x$i.sol[[i]]$solution))) {
                prettyNums <- formatC(
                    seq(length(x$i.sol[[i]]$solution)),
                    digits = nchar(length(x$i.sol[[i]]$solution)) - 1,
                    flag = 0
                )
                preamble <- paste("M", prettyNums[sol], ": ", sep = "")
                preamble <- paste(
                    preamble,
                    paste(
                        rep(" ", 7 - nchar(preamble)),
                        collapse = ""
                    ),
                    sep = ""
                )
                cat(preamble)
                xsol <- x$i.sol[[i]]$solution[[sol]]
                sufnec.char[i] <- paste(
                    ifelse(sufnec[i], "<", ""),
                    "->",
                    sep = ""
                )
                if (length(x$i.sol[[i]]$essential) > 0) {
                    xsol <- xsol[!is.element(xsol, x$i.sol[[i]]$essential)]
                    if (replace) {
                        x$i.sol[[i]]$essential <- replaceCategories(
                            x$i.sol[[i]]$essential,
                            categories = x$tt$categories
                        )
                    }
                    xsol <- paste(
                        paste(
                            x$i.sol[[i]]$essential,
                            collapse = "@"
                        ),
                        ifelse(
                            length(xsol) > 0,
                            paste(
                                "@(",
                                paste(xsol, collapse = "@"),
                                ")",
                                sep = ""
                            ),
                            ""
                        ),
                        sep = ""
                    )
                    cat(
                        admisc::prettyString(
                            unlist(
                                strsplit(xsol, split="@")
                            ),
                            line.length - 7,
                            7,
                            "+",
                            sufnec.char[i],
                            outcome
                        ),
                        "\n"
                    )
                }
                else {
                    if (replace) {
                        x$i.sol[[i]]$solution[[sol]] <- replaceCategories(
                            x$i.sol[[i]]$solution[[sol]],
                            categories = x$tt$categories
                        )
                    }
                    cat(
                        admisc::prettyString(
                            x$i.sol[[i]]$solution[[sol]],
                            line.length - 7,
                            7,
                            "+",
                            sufnec.char[i],
                            outcome
                        ),
                        "\n"
                    )
                }
            }
            if (x$options$details) {
                print.QCA_pof(x$i.sol[[i]]$IC, show.cases = x$options$show.cases)
            }
        }
    }
    else { 
        if (length(x$solution) == 1) {
            sufnec <- all(admisc::agteb(x$IC$sol.incl.cov[3], sol.cov))
            sufnec <- paste(ifelse(sufnec, "<", ""), "->", sep = "")
            if (replace) {
                x$solution[[1]] <- replaceCategories(
                    x$solution[[1]],
                    categories = x$tt$categories
                )
            }
            cat(
                sprintf(
                    "M1: %s\n",
                    admisc::prettyString(
                        x$solution[[1]],
                        line.length - 4,
                        4,
                        "+",
                        sufnec,
                        outcome
                    )
                )
            )
        }
        else {
            prettyNums <- formatC(
                seq(length(x$solution)),
                digits = nchar(length(x$solution)) - 1,
                flag = 0
            )
            sufnec <- logical(length(x$solution))
            for (i in seq(length(x$solution))) {
                sufnec[i] <- all(
                    admisc::agteb(
                        x$IC$individual[[i]]$sol.incl.cov[3],
                        sol.cov
                    )
                )
            }
            sufnec.char <- rep("", length(sufnec))
            for (i in seq(length(x$solution))) {
                cat(paste("M", prettyNums[i], ": ", sep = ""))
                xsol <- x$solution[[i]]
                sufnec.char[i] <- paste0(ifelse(sufnec[i], "<", ""), "->")
                if (length(x$essential) > 0) {
                    xsol <- xsol[!is.element(xsol, x$essential)]
                    if (replace) {
                        if (!essential_replaced) {
                            x$essential <- replaceCategories(
                                x$essential,
                                categories = x$tt$categories
                            )
                            essential_replaced <- TRUE
                        }
                        xsol <- replaceCategories(
                            xsol,
                            categories = x$tt$categories
                        )
                    }
                    xsol <- paste0(
                        paste(x$essential, collapse = "@"),
                        ifelse(
                            length(xsol) > 0,
                            paste(
                                "@(",
                                paste(xsol, collapse = "@"),
                                ")",
                                sep = ""
                            ),
                            ""
                        )
                    )
                    cat(
                        admisc::prettyString(
                            unlist(strsplit(xsol, split="@")),
                            line.length - nchar(prettyNums[i]) - 3,
                            nchar(prettyNums[i]) + 3,
                            "+",
                            sufnec.char[i],
                            outcome
                        ),
                        "\n"
                    )
                }
                else {
                    if (replace) {
                        x$solution[[i]] <- replaceCategories(
                            x$solution[[i]],
                            categories = x$tt$categories
                        )
                    }
                    cat(
                        admisc::prettyString(
                            x$solution[[i]],
                            line.length - nchar(prettyNums[i]) - 3,
                            nchar(prettyNums[i]) + 3,
                            "+",
                            sufnec.char[i],
                            outcome
                        ),
                        "\n"
                    )
                }
            }
            if (!mqca & x$options$details & enter) {
                cat("\n")
            }
        }
        if (x$options$details) {
            print.QCA_pof(
                x$IC,
                show.cases = x$options$show.cases,
                line.length = line.length
            )
        }
    }
    if (x$complex) {
        warning(
            simpleWarning(
                "The PI chart is exceedingly complex, solution(s) not guaranteed to be exhaustive.\n\n"
            )
        )
    }
    if (!x$options$details & enter) {
        cat("\n")
    }
}
`print.QCA_sS` <- function(x, ...) {
    dots <- list(...)
    xletters <- NULL
    if (x$options$use.letters) {
        xletters <- x$letters
    }
    if (x$options$use.labels && length(x$categories) > 0) {
        rownames(x$incl.cov) <- replaceCategories(
            rownames(x$incl.cov),
            categories = x$categories
        )
    }
    if (x$options$use.letters) {
        conditions <- names(xletters)
        xletters <- as.vector(xletters)
        if (!all(is.element(conditions, xletters))) {
            cat("\n")
            for (i in seq(length(xletters))) {
                cat(
                    "    ",
                    paste(xletters[i], ": ", sep = ""),
                    conditions[i],
                    "\n",
                    sep = ""
                )
            }
        }
    }
    incl.cov <- x$incl.cov
    cat("\n")
    prettyNums <- format(seq(nrow(incl.cov)))
    rownames(incl.cov) <- format(rownames(incl.cov))
    colnames(incl.cov) <- format(
        colnames(incl.cov),
        justify = "centre",
        width = 5
    )
    for (i in seq(ncol(incl.cov))) {
        NAs <- is.na(incl.cov[, i])
        incl.cov[!NAs, i] <- formatC(incl.cov[!NAs, i], digits = 3, format = "f")
        incl.cov[NAs, i] <- "  -  "
    }
    nchar.rownames <- nchar(rownames(incl.cov)[1])
    cat(
        paste(
            c(
                paste(
                    rep(" ", nchar.rownames + nchar(nrow(incl.cov)) + 2),
                    collapse = ""
                ),
                format(
                    colnames(incl.cov),
                    justify = "centre"
                )
            ),
            collapse = "  "
        ),
        "\n"
    )
    sep.row <- paste(
        rep("-", nchar.rownames + nchar(nrow(incl.cov)) + 7 * ncol(incl.cov) + 2),
        collapse = ""
    )
    cat(sep.row, "\n")
    for (i in seq(nrow(incl.cov))) {
        cat(
            paste(
                prettyNums[i],
                paste(
                    c(
                        rownames(incl.cov)[i],
                        incl.cov[i, ]
                    ),
                    collapse = "  "
                ),
                sep = "  "
            ),
            "\n"
        )
    }
    cat(sep.row, "\n")
    cat("\n")
}
`print.QCA_tt` <- function(x, ...) {
    dots <- list(...)
    enter <- ifelse (
        is.element("enter", names(as.list(x$call))),
        as.list(x$call)$enter,
        TRUE
    )
    if (!is.null(x$rowsorder)) {
        x$tt <- x$tt[x$rowsorder, ]
    }
    complete <- x$options$complete
    if (is.element("complete", names(dots))) {
        if (is.logical(dots$complete)) {
            complete <- dots$complete[1]
        }
    }
    show.cases <- x$options$show.cases
    if (is.element("show.cases", names(dots))) {
        if (is.logical(dots$show.cases)) {
            show.cases <- dots$show.cases[1]
        }
    }
    if (!complete) {
        if (!is.element("removed", names(x$options))) {
            x$tt <- x$tt[x$tt$OUT != "?", , drop = FALSE]
        }
    }
    if (show.cases) {
        if (x$options$dcc) {
            x$tt$cases <- ""
            x$tt[names(x$DCC), "cases"] <- x$DCC
            colnames(x$tt)[colnames(x$tt) == "cases"] <- "DCC"
        }
    }
    else {
        x$tt$cases <- NULL
    }
    if (nrow(x$tt) > 1024) {
        if (enter) cat("\n")
        cat(
            paste(
                "Warning: The truth table is too large (",
                nrow(x$tt),
                " rows). ",
                "Printing it on the screen is impractical.\n         ",
                "N.B.: You can still use its internal components (see ?str).",
                "\n\n",
                sep = ""
            )
        )
    }
    else {
        rownames(x$tt) <- paste(format(as.numeric(rownames(x$tt))), "")
        nofconditions <- length(x$noflevels)
        names.mydata <- colnames(x$recoded.data)[seq(nofconditions + 1)]
        if (!is.element("removed", names(x$options))) {
            if (enter) cat("\n")
            if (
                !all(
                    is.element(
                        names(x$tt)[seq(nofconditions)],
                        names(x$recoded.data)[seq(nofconditions)]
                    )
                )
            ) {
                for (i in seq(nofconditions)) {
                    cat("    ", paste(names(x$tt)[i], ": ", sep = ""), names.mydata[i], "\n", sep = "")
                }
            }
        }
        x$tt[, "n"] <- paste(" ", x$tt[, "n"], "")
        colnames(x$tt)[colnames(x$tt) == "n"] <- "  n "
        inclusion <- x$tt[, "incl"]
        missincl <- x$tt[, "incl"] == "-"
        x$tt[!missincl, "incl"] <- formatC(as.numeric(inclusion[!missincl]), digits = 3, format = "f")
        whichpri <- tail(which(colnames(x$tt) == "PRI"), 1)
        pri <- x$tt[, whichpri]
        misspri <- x$tt[, whichpri] == "-"
        x$tt[!misspri, whichpri] <- formatC(as.numeric(pri[!misspri]), digits = 3, format = "f")
        if (any(names(x$tt) == "pval1")) {
            x$tt[x$tt[, "pval1"] != "-", "pval1"] <- formatC(
                as.numeric(x$tt[x$tt[, "pval1"] != "-", "pval1"]),
                digits = 3,
                format = "f"
            )
            if (length(x$options$incl.cut) > 1) {
                x$tt[x$tt[, "pval0"] != "-", "pval0"] <- formatC(
                    as.numeric(x$tt[x$tt[, "pval0"] != "-", "pval0"]),
                    digits = 3,
                    format = "f"
                )
            }
        }
        if (any(missincl)) {
            x$tt[missincl, "incl"] <- "  -"
        }
        if (any(misspri)) {
            x$tt[misspri, "PRI"] <- "  -"
        }
        if (!is.element("removed", names(x$options))) {
            cat("  OUT: output value\n")
            cat("    n: number of cases in configuration\n")
            cat(" incl: sufficiency inclusion score\n")
            cat("  PRI: proportional reduction in inconsistency\n")
            if (show.cases & x$options$dcc) {
                cat("  DCC: deviant cases consistency\n")
            }
            if (any(names(x$tt) == "pval1")) {
                cat(
                    paste(
                        "pval1: p-value for alternative hypothesis inclusion > ",
                        x$options$incl.cut[1],
                        "\n",
                        sep = ""
                    )
                )
                if (length(x$options$incl.cut) > 1) {
                    cat(
                        paste(
                            "pval0: p-value for alternative hypothesis inclusion > ",
                            x$options$incl.cut[2],
                            "\n",
                            sep = ""
                        )
                    )
                }
            }
            cat("\n")
        }
        alloutzero <- all(x$tt$OUT == 0)
        wOUT <- which(colnames(x$tt) == "OUT")
        if (x$options$use.labels && length(x$categories) > 0) {
            outcome <- x$options$outcome
            for (fcond in setdiff(names(x$categories), outcome)) {
                x$tt[, fcond] <- x$categories[[fcond]][x$tt[, fcond] + 1]
            }
            if (is.element(outcome, names(x$categories))) {
                obs <- x$tt[, wOUT] != "?"
                x$tt[, wOUT][obs] <- x$categories[[outcome]][
                    as.numeric(x$tt[, wOUT][obs]) + 1
                ]
            }
        }
        x$tt[, wOUT] <- paste(" ", x$tt[, wOUT], "")
        colnames(x$tt)[wOUT] <- "  OUT "
        print(admisc::prettyTable(x$tt))
        if (alloutzero) {
            if (enter) cat("\n")
            cat(
                paste(
                    "It seems that all output values have been coded to zero.",
                    "Suggestion: lower the inclusion score for the presence of the outcome,",
                    sprintf(
                        "the relevant argument is \"incl.cut\" which now has a value of %s.\n", x$options$incl.cut[1]
                    ),
                    sep = "\n"
                )
            )
        }
        if (enter) cat("\n")
    }
}
