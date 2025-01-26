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

`XYplot` <- function(
    x, y, data, relation = "sufficiency", mguides = TRUE,
    jitter = FALSE, clabels = NULL, enhance = FALSE, model = FALSE, ...
) {
    dots <- list(...)
    funargs <- unlist(lapply(match.call(), deparse)[-1])
    if (missing(x)) {
        admisc::stopError(
            "Argument x is mandatory.", ... = ...
        )
    }
    x <- admisc::recreate(substitute(x))
    if (!missing(y)) {
        y <- admisc::recreate(substitute(y))
    }
    via.web <- FALSE
    if (length(testarg <- which(names(dots) == "via.web")) > 0) {
        via.web <- dots$via.web
        dots <- dots[-testarg]
    }
    negated <- logical(2)
    xname <- yname <- ""
    minus <- rawToChar(as.raw(c(226, 128, 147)))
    if (is.vector(drop(x)) & is.character(x) & any(grepl("\\$solution", funargs["x"]))) {
        x <- list(x)
    }
    if (is.list(x)) {
        if (any(grepl("\\$solution", funargs["x"]))) {
            model <- TRUE
            obj <- get(unlist(strsplit(funargs["x"], split = "[$]"))[1])
            data <- obj$tt$initial.data
            y <- obj$tt$options$outcome
            if (obj$tt$options$neg.out) {
                y <- paste("~", y, sep = "")
            }
            x <- paste(unlist(x), collapse = " + ") 
        }
    }
    if (!is.data.frame(x) & !is.matrix(x) & !missing(y)) {
        testit <- capture.output(tryCatch(eval(y), error = function(e) e))
        if (length(testit) == 1 & is.character(testit)) {
            if (grepl("Error", testit)) {
                y <- as.vector(funargs["y"])
            }
        }
    }
    if (is.character(x)) {
        if (length(x) == 1) {
            x <- admisc::splitstr(x)
        }
        if (length(x) == 1) {
            x <- unlist(strsplit(x, split = "->|=>"))
            if (length(x) == 1) {
                x <- unlist(strsplit(x, split = "<-|<="))
                if (length(x) > 1) {
                    relation <- "necessity"
                    y <- admisc::trimstr(x[2])
                    x <- admisc::trimstr(x[1])
                }
            }
            else {
                y <- admisc::trimstr(x[2])
                x <- admisc::trimstr(x[1])
            }
            if (missing(y)) {
                admisc::stopError(
                    "The outcome's name is missing.", ... = ...
                )
            }
            else if (!is.character(y)) {
                admisc::stopError(
                    "Unknown <x> and/or <y> arguments.", ... = ...
                )
            }
        }
        else {
            if (!missing(y)) {
                if (is.data.frame(y)) {
                    data <- y
                }
            }
            y <- x[2]
            x <- x[1]
        }
        if (missing(data)) {
            admisc::stopError(
                "Data is missing.", ... = ...
            )
        }
        else {
            verify.qca(data)
        }
        xname <- as.character(parse(text = x))
        yname <- as.character(parse(text = y))
        x <- gsub(minus, "-", gsub("[[:space:]]", "", x))
        y <- gsub(minus, "-", gsub("[[:space:]]", "", y))
        negated <- logical(2)
        negated[1] <- identical(unname(substring(x, 1, 2)), "1-")
        negated[2] <- identical(unname(substring(y, 1, 2)), "1-")
        if (any(checks <- grepl("1-", c(x, y)) & !negated)) {
            admisc::stopError(
                paste0(
                    "Incorrect expression in \"",
                    paste(
                        c(x, y)[checks],
                        collapse = "\" and \""
                    ),
                    "\"."
                ),
                ... = ...
            )
        }
        x <- admisc::compute(x, data = data)
        y <- admisc::compute(y, data = data)
        negated <- logical(2)
    }
    else if (is.data.frame(x) | is.matrix(x)) {
        verify.qca(as.data.frame(x))
        if (ncol(x) < 2) {
            admisc::stopError(
                "At least two columns are needed.", ... = ...
            )
        }
        xname <- colnames(x)[1]
        yname <- colnames(x)[2]
        y <- x[, 2]
        x <- x[, 1]
    }
    else if (!missing(y)) {
        if (length(x) > 1 & is.numeric(x)) { 
            oneminus <- identical(unname(substring(gsub("[[:space:]]", "", funargs[1]), 1, 2)), "1-")
            if (any((admisc::hastilde(funargs[1])    & !admisc::tilde1st(funargs[1])) | 
                    (grepl("1-", funargs[1]) & !oneminus)
                   )) {
                admisc::stopError(
                    paste0("Incorrect expression in \"", funargs[1], "\"."),
                    ... = ...
                )
            }
            negated[1] <- oneminus | admisc::tilde1st(funargs[1])
            xname <- "X"
            tc <- capture.output(
                tryCatch(
                    admisc::getName(funargs[1]),
                    error = function(e) e,
                    warning = function(w) w
                )
            )
            if (!grepl("simpleError", tc)) {
                xname <- admisc::notilde(admisc::getName(funargs[1]))
            }
        }
        if (length(y) > 1 & is.numeric(y)) { 
            oneminus <- identical(unname(substring(gsub("[[:space:]]", "", funargs[2]), 1, 2)), "1-")
            if (any((admisc::hastilde(funargs[2])    & !admisc::tilde1st(funargs[2])) | 
                    (grepl("1-", funargs[2]) & !oneminus)
                   )) {
                admisc::stopError(
                    paste0("Incorrect expression in \"", funargs[2], "\"."),
                    ... = ...
                )
            }
            negated[2] <- oneminus | admisc::tilde1st(funargs[2])
            yname <- "Y"
            tc <- capture.output(
                tryCatch(
                    admisc::getName(funargs[2]),
                    error = function(e) e,
                    warning = function(w) w
                )
            )
            if (!grepl("simpleError", tc)) {
                yname <- admisc::notilde(admisc::getName(funargs[2]))
            }
        }
        if (length(y) == 1 & is.character(y)) {
            if (missing(data)) {
                admisc::stopError(
                    "Data is missing.", ... = ...
                )
            }
            else {
                verify.qca(data)
            }
            yname <- as.character(parse(text = y))
            y <- gsub(minus, "-", gsub("[[:space:]]", "", y))
            negated[2] <- identical(unname(substring(y, 1, 2)), "1-")
            if (grepl("1-", y) & !negated[2]) {
                admisc::stopError(
                    paste0("Incorrect expression in \"", y, "\"."),
                    ... = ...
                )
            }
            y <- admisc::compute(y, data = data)
            negated[2] <- FALSE
        }
    }
    else {
        admisc::stopError(
            "Either a dataframe with two columns or two vectors are needed.",
            ... = ...
        )
    }
    if (any(x > 1) | any(y > 1)) {
        admisc::stopError(
            "Values should be bound between 0 and 1.", ... = ...
        )
    }
    xcopy <- x
    ycopy <- y
    if (is.element("QCA_fuzzy", class(xcopy))) {
        attributes(xcopy) <- NULL
    }
    if (is.element("QCA_fuzzy", class(ycopy))) {
        attributes(ycopy) <- NULL
    }
    jitfactor <- 0.01
    jitamount <- 0.01
    cexaxis <- 0.8
    hadj <- 1.1
    padj <- 0
    linex <- 1.75
    liney <- 2
    linet <- 1.5
    pch <- rep(21, length(x))
    cexpoints <- rep(0.8, length(x))
    bgpoints <- rep("#707070", length(x)) # "#ababab"
    if (length(testarg <- which(names(dots) == "pch")) > 0) {
        pch <- dots$pch
        if (length(pch) == 1) {
            pch <- rep(pch, length(x))
        }
        else {
            if (length(pch) != length(x)) {
                admisc::stopError(
                    sprintf(
                        "Length of argument \"pch\" different from the %s.",
                        ifelse(
                            missing(data),
                            "length of \"x\"",
                            "number of rows in the data"
                        )
                    ),
                    ... = ...
                )
            }
        }
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "cex")) > 0) {
        cexpoints <- dots$cex
        if (length(cexpoints) == 1) {
            cexpoints <- rep(cexpoints, length(x))
        }
        else {
            if (length(cexpoints) != length(x)) {
                admisc::stopError(
                    sprintf(
                        "Length of argument \"cex\" different from the %s.",
                        ifelse(
                            missing(data),
                            "length of \"x\"",
                            "number of rows in the data"
                        )
                    ),
                    ... = ...
                )
            }
        }
        dots <- dots[-testarg]
    }
    bginput <- is.element("bg", names(dots))
    if (length(testarg <- which(names(dots) == "bg")) > 0) {
        bgpoints <- dots$bg
        if (length(bgpoints) == 1) {
            bgpoints <- rep(bgpoints, length(x))
        }
        else {
            if (length(bgpoints) != length(x)) {
                admisc::stopError(
                    sprintf(
                        "Length of argument \"bg\" different from the %s.",
                        ifelse(
                            missing(data),
                            "length of \"x\"",
                            "number of rows in the data"
                        )
                    ),
                    ... = ...
                )
            }
        }
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
    if (length(testarg <- which(names(dots) == "hadj")) > 0) {
        hadj <- dots$hadj
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "padj")) > 0) {
        padj <- dots$padj
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "line")) > 0) {
        linex <- dots$line[1]
        liney <- ifelse(is.na(dots$line[2]), dots$line[1], dots$line[2])
        linet <- ifelse(is.na(dots$line[3]), dots$line[1], dots$line[3])
        dots <- dots[-testarg]
    }
    if (!is.null(clabels)) {
        if (is.numeric(clabels)) {
            if (length(clabels) < length(x)) {
                if (all(clabels <= length(x))) {
                    rownms <- rep("", length(x))
                    rownms[clabels] <- clabels
                    clabels <- rownms
                }
                else {
                    admisc::stopError(
                        "Values in the argument <clabels> outside the data rows.",
                        ... = ...
                    )
                }
            }
            clabels <- as.character(clabels)
        }
        if (length(clabels) != length(x)) {
            admisc::stopError(
                sprintf(
                    "Length of argument <clabels> larger than %s.",
                    ifelse(
                        missing(data),
                        "length of <x>",
                        "number of rows in the data"
                    )
                ),
                ... = ...
            )
        }
        if (is.logical(clabels)) {
            if (missing(data)) {
                rownms <- seq(length(x))
            }
            else {
                rownms <- rownames(data)
            }
            rownms[!clabels] <- ""
            clabels <- rownms
        }
    }
    cexlabels <- cexpoints
    if (enhance) {
        if (is.null(clabels)) {
            caselabels <- rep("", length(x))
        }
        if (relation == "sufficiency") {
            if (any(selection <- x >= 0.5 & y >= 0.5 & x <= y)) {
                if (is.null(clabels) & !model) {
                    if (missing(data)) {
                        caselabels[selection] <- which(selection)
                    }
                    else {
                        caselabels[selection] <- rownames(data)[selection]
                    }
                }
                xs <- x[selection]
                ys <- y[selection]
                pch[which(selection)][which.min((ys - xs)/xs)] <- 3
            }
            if (any(selection <- x >= 0.5 & y >= 0.5 & x > y)) {
                if (is.null(clabels) & !model) {
                    if (missing(data)) {
                        caselabels[selection] <- which(selection)
                    }
                    else {
                        caselabels[selection] <- rownames(data)[selection]
                    }
                }
                if (!bginput) {
                    bgpoints[selection] <- "#cccccc"
                }
            }
            if (any(selection <- x >= 0.5 & y < 0.5)) {
                xs <- x[selection]
                ys <- y[selection]
                pch[selection] <- 23
                if (!bginput) {
                    bgpoints[which(selection)][
                        which.min(1 - (ys - xs)/xs)
                    ] <- "#cccccc"
                }
            }
            if (any(selection <- x < 0.5 & y < 0.5)) {
                if (is.null(clabels) & model) {
                    caselabels[selection] <- rownames(data)[selection]
                }
                cexpoints[selection] <- 0.875 * cexpoints[selection]
                pch[selection] <- 24
                if (!bginput) {
                    bgpoints[selection] <- "#cccccc"
                }
            }
            if (any(selection <- x < 0.5 & y >= 0.5)) {
                if (is.null(clabels) & model) {
                    caselabels[selection] <- rownames(data)[selection]
                }
                pch[selection] <- 22
                if (!bginput) {
                    bgpoints[selection] <- "#cccccc"
                }
            }
        }
        if (is.null(clabels)) {
            clabels <- caselabels
        }
    }
    if (jitter) {
        x <- jitter(x, jitfactor, jitamount)
        y <- jitter(y, jitfactor, jitamount)
    }
    toplot <- list(x = x, y = y)
    xlabel <- paste0(ifelse(negated[1], "~", ""), xname)
    ylabel <- paste0(ifelse(negated[2], "~", ""), yname)
    if (model) xlabel <- "MODEL"
    if (length(testarg <- which(names(dots) == "xlab")) > 0) {
        xlabel <- dots$xlab
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "ylab")) > 0) {
        ylabel <- dots$ylab
        dots <- dots[-testarg]
    }
    if (length(testarg <- which(names(dots) == "cex.axis")) > 0) {
        cexaxis <- dots$cex.axis
        dots <- dots[-testarg]
    }
    toplot$type <- "n"
    toplot$xlim <- c(0, 1)
    toplot$ylim <- c(0, 1)
    toplot$xlab <- ""
    toplot$ylab <- ""
    toplot$axes <- FALSE
    if (length(dots) > 0) {
        toplot <- c(toplot, dots)
    }
    par(mar = c(3, 3.1, 2.5, 0.5), cex.axis = cexaxis, tck = -.015,
        las = 1, xpd = FALSE, mgp = c(1.5, 0.5, 0))
    suppressWarnings(do.call("plot", toplot))
    box()
    axis(1, xaxp = c(0, 1, 10), padj = padj)
    axis(2, yaxp = c(0, 1, 10), hadj = hadj)
    title(xlab = xlabel, cex.lab = cexaxis + 0.1, font.lab = 2, line = linex)
	
    title(ylab = ylabel, cex.lab = cexaxis + 0.1, font.lab = 2, line = liney)
	title(
        main = paste(
            ifelse(
                nec(relation),
                "Necessity",
                "Sufficiency"
            ),
            "relation"
            ),
        cex.main = cexaxis/0.8,
        font.main = 2,
        line = linet
    )
    if (mguides) {
        abline(v = .5, lty = 2, col = "gray")
        abline(h = .5, lty = 2, col = "gray")
    }
    abline(0, 1, col = "gray")
    plotpoints <- list(
        x,
        y,
        pch = pch,
        cex = cexpoints,
        bg = bgpoints
    ) 
    suppressWarnings(do.call("points", c(plotpoints, dots)))
    inclcov <- round(
        pof(
            setms = xcopy,
            outcome = ycopy,
            relation = relation
        )$incl.cov[1, 1:3],
        3
    )
    inclcov[is.na(inclcov)] <- 0
    inclcov <- sprintf("%.3f", inclcov)
    mtext(
        paste(
            c(
                "Inclusion:",
                "Coverage:",
                ifelse(nec(relation), "Relevance:", "PRI:")
            ),
            inclcov[c(1, 3, 2)],
            collapse = "   "
        ),
        at = 0,
        adj = 0,
        cex = cexaxis
    )
    cexl <- ifelse(any(names(dots) == "cex"), dots$cex, 1)
    srtl <- ifelse(any(names(dots) == "srt"), dots$srt, 0)
    if (!is.null(clabels)) {
        text(
            x,
            y + 0.02,
            labels = clabels,
            srt = srtl,
            cex = cexlabels * cexl
        )
    }
}
