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

library(shiny)
library(QCA)
library(admisc)
library(tools)
library(venn)
setwd(Sys.getenv("userwd"))
options(help_type = "html")
listFiles <- function(dirpath, filetype = "*") {
    result <- list(dirs = NULL, files = NULL, filepath = as.matrix(filepath), ok = oktoset)
    found <- list.files(dirpath)
    temp <- toupper(found)
    found <- found[match(sort(temp), temp)]
    if (length(found) > 0) {
        isdir <- sapply(found, function(x) file_test("-d", file.path(dirpath, x)))
        if (any(isdir)) {
            result$dirs <- found[isdir]
        }
        if (any(!isdir)) {
            if (filetype != "*") {
                extensions <- unlist(lapply(strsplit(found, split="[.]"), function(x) x[length(x)]))
                found <- found[which(toupper(extensions) == toupper(filetype))]
                if (length(found) > 0) {
                    result$files <- found
                }
            }
            else {
                result$files <- found[!isdir]
            }
        }
        result$dirs <- as.list(result$dirs)
        result$files <- as.list(result$files)
    }
    resfilename <- ""
    if (!identical(filepath, "")) {
        if (file_test("-f", filepath)) {
            extension <- file_ext(basename(filepath))
            resfilename <- gsub("[[:space:]]", "_", file_path_sans_ext(basename(filepath)))
            if (possibleNumeric(substr(resfilename, 1, 1))) {
                resfilename <- paste("x", resfilename, sep="")
            }
        }
    }
    result$filename <- resfilename
    result$extension <- extension
    result$wd <- getwd()
    return(result)
}
infobjs <- function(env, objs, scrollvh) {
    misscroll <- missing(scrollvh)
    toreturn <- list(data = NULL, tt = NULL, qmc = NULL)
    if (length(objs) > 0) {
        objs <- unlist(lapply(mget(objs, env), function(x) {
            if (is.data.frame(x)) {
                return(1)
            }
            else if (methods::is(x, "QCA_tt")) {
                return(2)
            }
            else if (methods::is(x, "QCA_min")) {
                return(3)
            }
            else {
                return(0)
            }
        }))
        if (any(objs == 1)) {
            toreturn$data <- lapply(names(objs[objs == 1]), function(n) {
                x <- env[[n]]
                dscrollvh <- c(1, 1)
                if (!misscroll) {
                    if (is.element(n, names(scrollvh))) {
                        dscrollvh <- scrollvh[[n]]
                    }
                }
                nrowd <- nrow(x)
                ncold <- ncol(x)
                srow <- min(dscrollvh[1], nrowd - min(nrowd, visiblerows) + 1)
                scol <- min(dscrollvh[2], ncold - min(ncold, visiblecols) + 1)
                erow <- min(srow + visiblerows - 1, nrowd)
                ecol <- min(scol + visiblecols - 1, ncold)
                list(
                    nrows = nrowd,
                    ncols = ncold,
                    rownames = as.list(rownames(x)),
                    colnames = as.list(colnames(x)),
                    numerics = as.list(as.vector(unlist(lapply(x, admisc::possibleNumeric)))),
                    calibrated = as.list(as.vector(unlist(lapply(x, function(x) {
                        if (is.numeric(x)) {
                            return(all(na.omit(x) >= 0 & na.omit(x) <= 1))
                        }
                        return(FALSE)
                    })))),
                    binary = as.list(as.vector(unlist(lapply(x, function(x) all(is.element(x, 0:1)))))),
                    scrollvh = c(srow, scol) - 1, 
                    theData = unname(as.list(x[seq(srow, erow), seq(scol, ecol), drop = FALSE])),
                    dataCoords = paste(srow, scol, erow, ecol, ncol(x), sep = "_")
                )
            })
            names(toreturn$data) <- names(objs[objs == 1])
        }
        if (any(objs == 2)) {
            toreturn$tt <- lapply(mget(names(objs[objs == 2]), env), function(x) {
                components <- c("indexes", "noflevels", "cases", "options", "colnames", "numerics")
                x$indexes <- x$indexes - 1 
                cnds <- x$options$conditions
                if (x$options$use.letters) {
                    cnds <- LETTERS[seq(length(cnds))]
                }
                x$options$outcome <- list(notilde(x$options$outcome))
                if (length(x$options$incl.cut) == 1) {
                    x$options$incl.cut <- list(x$options$incl.cut)
                }
                if (length(cnds) <= 7) {
                    x$id <- apply(x$tt[, cnds], 1, function(x) {
                        ifelse(any(x == 1), paste(which(x == 1), collapse=""), "0")
                    })
                    components <- c(components, "id", "tt")
                }
                x$colnames <- colnames(x$initial.data)
                x$numerics <- as.vector(unlist(lapply(x$initial.data, admisc::possibleNumeric)))
                return(x[components])
            })
        }
        if (any(objs == 3)) {
            toreturn$qmc <- lapply(mget(names(objs[objs == 3]), env), function(x) {
                components <- c("indexes", "noflevels", "cases", "options")
                x <- x$tt
                cnds <- x$options$conditions
                if (x$options$use.letters) {
                    cnds <- LETTERS[seq(length(cnds))]
                }
                if (length(cnds) <= 7) {
                    x$id <- apply(x$tt[, cnds], 1, function(x) {
                        ifelse(any(x == 1), paste(which(x == 1), collapse=""), "0")
                    })
                    components <- c(components, "id", "tt")
                }
                x$indexes <- x$indexes - 1 
                return(x[components])
            })
        }
    }
    return(toreturn)
}
calibrateit <- function(foo) {
    scrollvh <- lapply(foo$scrollvh, function(x) unlist(x) + 1)
    checks <- rep(TRUE, 9)
    checks[1] <- !is.null(ev[[foo$dataset]])
    foo$thresholds <- unlist(foo$thresholds)
    nms <- unlist(foo$thnames)[foo$thresholds != ""]
    foo$thresholds <- foo$thresholds[foo$thresholds != ""]
    thrs <- suppressWarnings(as.numeric(foo$thresholds))
    if (any(!is.na(thrs))) {
        foo$thresholds <- as.numeric(thrs[!is.na(thrs)])
    }
    if (all(foo$thresholds == "")) {
        foo$thresholds <- NA
    }
    else {
        if (!is.null(nms)) {
            if (length(nms) == length(foo$thresholds)) {
                if (foo$type == "fuzzy") {
                    names(foo$thresholds) <- nms
                }
            }
            else {
                foo$thresholds <- NA
            }
        }
    }
    checks[4] <- foo$x != ""
    if (checks[1] & checks[4]) {
        if (foo$x %in% names(ev[[foo$dataset]])) {
            checks[6] <- is.numeric(ev[[foo$dataset]][, foo$x])
        }
    }
    if (possibleNumeric(foo$idm)) {
        foo$idm <- as.numeric(foo$idm)
    }
    if (possibleNumeric(foo$below)) {
        foo$below <- as.numeric(foo$below)
    }
    if (possibleNumeric(foo$above)) {
        foo$above <- as.numeric(foo$above)
    }
    if (all(checks)) {
        textoutput <- capture.output(tryCatch(
            calibrate(
                ev[[foo$dataset]][, foo$x],
                type = foo$type,
                thresholds = foo$thresholds,
                logistic = foo$logistic,
                idm = foo$idm,
                ecdf = foo$ecdf,
                below = foo$below,
                above = foo$above), error = function(e) e)
        )
        response <- list()
        response$origin <- "calibrate"
        response$error <- FALSE
        response$toprint <- ""
        if (any(error <- grepl("Error", textoutput))) {
            errmessage <- paste0("Error:", unlist(strsplit(textoutput[which(error)], split=":"))[2])
            errmessage <- substr(errmessage, 1, nchar(errmessage) - 1)
            response$error <- TRUE
            response$toprint <- errmessage
        }
        else {
            fuzzyvals <- calibrate(
                ev[[foo$dataset]][, foo$x],
                type = foo$type,
                thresholds = foo$thresholds,
                logistic = foo$logistic,
                idm = foo$idm,
                ecdf = foo$ecdf,
                below = foo$below,
                above = foo$above)
            if (!foo$thsetter) {
                ev[[foo$dataset]][, ifelse(foo$newvar != "", foo$newvar, foo$x)] <- fuzzyvals
            }
        }
        if (foo$thsetter) {
            if (!response$error) {
                return(list(fuzzyvals = fuzzyvals))
            }
        }
        else {
            response$infobjs <- infobjs(ev, foo$dataset, scrollvh)
            response$toprint <- paste(response$toprint, collapse = "<br>") 
            response$dataset <- foo$dataset
            if (foo$newvar == "") { 
                response$poinths <- list(dataset = foo$dataset,
                                     condition = foo$x,
                                     vals = unname(ev[[foo$dataset]][, foo$x]),
                                     fuzzyvals = fuzzyvals)
                if (foo$findth) {
                    response$poinths$thvals <- as.list(findTh(ev[[foo$dataset]][, foo$x], n = foo$nth))
                    response$poinths$message <- "OK"
                }
            }
            return(list(response = response))
        }
    }
}
tryCatchWEM <- function(expr) {
    toreturn <- list()
    output <- withVisible(withCallingHandlers(
        tryCatch(expr, error = function(e) {
            toreturn$error <<- e$message
            NULL
        }), warning = function(w) {
            toreturn$warning <<- c(toreturn$warning, w$message)
            invokeRestart("muffleWarning")
        }, message = function(m) {
            toreturn$message <<- paste(toreturn$message, m$message, sep = "")
            invokeRestart("muffleMessage")
        }
    ))
    if (output$visible) {
        if (!is.null(output$value)) {
            toreturn$output <- capture.output(output$value)
        }
    }
    if (length(toreturn) > 0) {
        return(toreturn)
    }
}
continue <- lapply(c("ls(", "'ls"), function(x) {
    x <- unlist(strsplit(unlist(strsplit(tryCatch(eval(parse(text = x)), error = identity)$message, "\n"))[1], ":"))
    return(admisc::trimstr(x[length(x)]))
})
evalparse <- function(foo) {
    forbidden <- "dev.new\\(|plot.new\\(|plot.window\\(|X11\\(|quartz\\(|dev.set\\(|windows\\("
    evaluateit <- list()
    tocheck <- rep(c("\n", "@$%$@"), each = 3)
    names(tocheck) <- c("output", "view", "message", "warning", "error", "library")
    if (any(ggplot <- grep("ggplot\\(|qplot\\(|quickplot\\(", foo))) {
        foo[ggplot] <- paste("SOME__BIG__NAME <-", foo[ggplot])
        foo <- c(foo, "print(SOME__BIG__NAME)")
    }
    ggplot <- any(ggplot)
    nextpos <- 1
    temptxt <- foo[1]
    parsem <- NULL
    for (i in seq_along(foo)) {
        endcom <- TRUE
        temp <- tryCatch(parse(text = temptxt), error = identity)
        if (inherits(temp, "error")) {
            if (any(unlist(lapply(continue, grepl, temp)))) {
                if (i < length(foo)) {
                    temptxt <- paste(foo[seq(nextpos, i + 1)], collapse = "\n")
                    endcom <- FALSE
                }
            }
        }
        if (endcom) {
            parsem <- c(parsem, temptxt)
            nextpos <- i + 1
            temptxt <- foo[nextpos]
        }
    }
    parsem_length <- length(parsem)
    if (identical(parsem[parsem_length], "print(SOME__BIG__NAME)")) {
        parsem_length <- parsem_length - 1
    }
    evaluateit <- vector(mode = "list", length = parsem_length)
    recall <- vector(mode = "list", length = parsem_length)
    noerror <- TRUE
    i <- 1
    while (i <= length(parsem) & noerror) {
        evaluateit[[i]] <- list()
        if (grepl(forbidden, gsub(" ", "", parsem[i]))) {
            evaluateit[[i]][["error"]] <- " Opening multiple graphics devices is not supported."
        }
        else {
            sf <- srcfile("parsem[i]")
            temp <- tryCatch(parse(text = parsem[i], srcfile = sf), error = identity, warning = identity)
            if (inherits(temp, "error")) {
                gpd <- getParseData(sf)
                parsesplit <- parsem[i]
                if (any(gpd$token == "';'")) {
                    pos <- c(0, gpd$col1[gpd$token == "';'"], nchar(parsem[i]) + 1)
                    parsesplit <- c()
                    for (p in seq(length(pos) - 1)) {
                        parsesplit <- c(parsesplit, substr(parsem[i], pos[p] + 1, pos[p + 1] - 1))
                    }
                }
                errors <- logical(length(parsesplit))
                for (j in seq(length(parsesplit))) {
                    errors[j] <- inherits(tryCatch(parse(text = parsesplit[j]), error = identity), "error")
                }
                parsesplit <- parsesplit[seq(which(errors)[1])]
            }
            else {
                ptfep <- function(x) {
                    x <- unlist(strsplit(x, split = ""))
                    what <- unlist(strsplit("eval(parse(text=", split = ""))
                    counter <- length(what)
                    while (counter > 0) {
                        if (x[1] == " ") {
                            x <- x[-1]
                        }
                        else {
                            if (x[1] == what[1]) {
                                x <- x[-1]
                                if (what[1] == "(") {
                                    x <- x[-length(x)]
                                }
                                what <- what[-1]
                                counter <- counter - 1
                            }
                        }
                    }
                    x <- admisc::trimstr(paste(x, collapse = ""))
                    x <- gsub("^\\\"", "", x)
                    x <- gsub("\\\"$", "", x)
                    x <- gsub("\\\\\"", "\"", x)
                    x <- admisc::trimstr(x)
                    if (grepl("^eval\\(parse\\(text=", gsub(" ", "", x))) {
                        return(Recall(x))
                    }
                    else {
                        x <- gsub("\\\\\"", "\"", x)
                        return(x)
                    }
                }
                parsesplit <- as.character(parse(text = parsem[i]))
            }
            if (length(parsesplit) > 1) {
                recalsplit <- vector(mode = "list", length = length(parsesplit))
                for (j in seq(length(recalsplit))) {
                    recalsplit[[j]] <- Recall(parsesplit[j])
                    temp <- list()
                    for (tc in names(tocheck)) {
                        tempcheck <- unlist(lapply(recalsplit[[j]], function(x) return(x[[tc]])))
                        if (!is.null(tempcheck)) {
                            temp[[tc]] <- paste(tempcheck, collapse = tocheck[[tc]])
                        }
                    }
                    recalsplit[[j]] <- temp
                }
                temp <- list()
                for (tc in names(tocheck)) {
                    tempcheck <- unlist(lapply(recalsplit, function(x) return(x[[tc]])))
                    if (!is.null(tempcheck)) {
                        temp[[tc]] <- paste(tempcheck, collapse = tocheck[[tc]])
                    }
                }
                errors <- unlist(lapply(recalsplit, function(x) return("error" %in% names(x))))
                if (any(errors)) {
                    temp$partial <- paste(paste(parsesplit[seq(which(errors)[1] - 1)], collapse = ";"), ";", sep="")
                }
                recall[[i]] <- temp
            } 
            else {
                evaluateit[[i]][["source"]] <- grepl("^source\\(", gsub(" ", "", parsem[i]))
                if (evaluateit[[i]][["source"]]) {
                    sourcefile <- gsub("\\\"", "", insideBrackets(parsem[i], type = "("))
                    if (file.exists(sourcefile)) {
                        sourcefile <- readLines(sourcefile, warn = FALSE)
                        temp <- tryCatch(parse(text = sourcefile), error = identity, warning = identity)
                        if (inherits(temp, "error")) {
                            errmsg <- temp$message
                            if (grepl("<text>", errmsg)) {
                                errmsg <- unlist(strsplit(errmsg, split = "\n"))
                                errmsg[1] <- paste("in ", parsem[i], ": ", unlist(strsplit(errmsg[1], split = ": "))[2], sep = "")
                                if (!grepl("unexpected symbol", errmsg[1])) {
                                    errmsg[length(errmsg)] <- substring(errmsg[length(errmsg)], 2)
                                }
                                errmsg <- paste(errmsg, collapse = "\n")
                            }
                            recall[[i]] <- list(error = errmsg)
                            noerror <- FALSE
                        }
                        else {
                            recall[[i]] <- Recall(sourcefile)
                            temp <- list()
                            for (tc in names(tocheck)) {
                                tempcheck <- unlist(lapply(recall[[i]], function(x) return(x[[tc]])))
                                if (!is.null(tempcheck)) {
                                    temp[[tc]] <- paste(tempcheck, collapse = tocheck[[tc]])
                                }
                            }
                            recall[[i]] <- temp
                            if (any(names(temp) == "error")) {
                                noerror <- FALSE
                            }
                        }
                    }
                    else {
                        errmsg <- paste(c("Error in file(filename, \"r\", encoding = encoding) : ",
                                          "   cannot open the connection",
                                          "In addition: Warning message:",
                                          "In file(filename, \"r\", encoding = encoding) :",
                                          paste("  cannot open file '", filename, "': No such file or directory", sep = "")),
                                        collapse = "\n")
                    }
                }
                else {
                    comnd <- parsem[i]
                    if (grepl("^print\\(", gsub(" ", "", comnd)) & !ggplot) {
                        comnd <- insideBrackets(comnd, type = "(")
                    }
                    object <- ""
                    objerror <- FALSE
                    toview <- FALSE
                    lib <- FALSE
                    if (grepl("^View\\(", gsub(" ", "", comnd))) {
                        object <- gsub("\\\"", "", admisc::insideBrackets(comnd, type = "("))
                        if (object != "") {
                            if (exists(object, envir = ev)) {
                                if (is.data.frame(ev[[object]])) {
                                    evaluateit[[i]][["view"]] <- object
                                    toview <- TRUE
                                }
                                else {
                                    objerror <- TRUE
                                    evaluateit[[i]][["error"]] <- paste("The object to View (", object, ") is not a dataframe.", sep="")
                                }
                            }
                        }
                    }
                    if (!objerror & !toview) {
                        if (length(temp <- tryCatchWEM(eval(parse(text = comnd), envir = ev))) > 0) {
                            if (any(names(temp) == "warning")) {
                                temp$warning <- paste(paste("In", parsem[i], ":"), temp$warning, sep = "\n  ")
                            }
                            if (any(names(temp) == "error")) {
                                temp$error <- paste("Error in", comnd, ":\n ", temp$error)
                            }
                            evaluateit[[i]] <- temp
                        }
                    }
                    if (grepl("^library\\(", gsub(" ", "", comnd))) {
                        object <- gsub("\\\"", "", admisc::insideBrackets(comnd, type = "("))
                        if (object != "") {
                            evaluateit[[i]][["library"]] <- object
                        }
                    }
                }
            }
            if (!is.null(evaluateit[[i]][["error"]])) {
                if (grepl("plot.new has not been called yet", evaluateit[[i]][["error"]])) {
                    if (!identical(grafic, emptyplot)) {
                        pdf(templotfile)
                        dev.control("enable")
                        replayPlot(grafic)
                        evaluateit[[i]][["output"]] <- capture.output(suppressWarnings(eval(parse(text = parsem[i]), envir = ev)))
                        if (length(evaluateit[[i]][["output"]]) == 0) {
                            evaluateit[[i]][["output"]] <- NULL
                        }
                        if (length(dev.list()) > 0) {
                            testplot <<- recordPlot()
                            sapply(dev.list(), dev.off)
                        }
                        evaluateit[[i]][["error"]] <- NULL
                        svg(filename = svgfile, width = plotsize[1], height = plotsize[2])
                        replayPlot(grafic)
                        dev.off()
                    }
                }
                else {
                    noerror <- FALSE
                }
            }
        }
        if (length(recall[[i]]) > 0) {
            for (tc in names(tocheck)) {
                if (tc %in% names(recall[[i]])) {
                    if (any(tc %in% names(evaluateit[[i]]))) {
                        evaluateit[[i]][[tc]] <- paste(c(evaluateit[[i]][[tc]], recall[[i]][[tc]]), collapse = tocheck[[tc]])
                    }
                    else {
                        evaluateit[[i]][[tc]] <- paste(recall[[i]][[tc]], collapse = tocheck[[tc]])
                    }
                }
            }
            if ("partial" %in% names(recall[[i]])) {
                evaluateit[[i]][["partial"]] <- recall[[i]][["partial"]]
            }
            if ("library" %in% names(recall[[i]])) {
                evaluateit[[i]][["library"]] <- recall[[i]][["library"]]
            }
        }
        if (!ggplot) {
            evaluateit[[i]][["command"]] <- gsub("SOME__BIG__NAME <- ", "", parsem[i])
            evaluateit[[i]][["continue"]] <- FALSE
            errmsg <- evaluateit[[i]][["error"]]
            if (!is.null(errmsg)) {
                if (any(unlist(lapply(continue, grepl, errmsg)))) {
                    sf <- srcfile("parsem[i]")
                    temp <- tryCatch(parse(text = parsem[i], srcfile = sf), error = identity, warning = identity)
                    if (inherits(temp, "error")) {
                        evaluateit[[i]][["continue"]] <- TRUE
                    }
                }
                if (grepl("<text>", errmsg)) {
                    errmsg <- unlist(strsplit(errmsg, split = "\n"))
                    errmsg[1] <- paste("Error:", unlist(strsplit(errmsg[1], split = ": "))[2], "in:")
                    if (!grepl("unexpected symbol", errmsg[1])) {
                        errmsg[length(errmsg)] <- substring(errmsg[length(errmsg)], 2)
                    }
                    errmsg <- paste(errmsg, collapse = "\n")
                }
                evaluateit[[i]][["error"]] <- errmsg
            }
            if (length(evaluateit[[i]][["output"]]) == 0) {
                evaluateit[[i]][["output"]] <- NULL
            }
            else {
                evaluateit[[i]][["output"]] <- paste(evaluateit[[i]][["output"]], collapse = "\n")
            }
            if ("view" %in% names(evaluateit[[i]])) {
                views <- unlist(strsplit(evaluateit[[i]][["view"]], split = "\n"))
                evaluateit[[i]][["view"]] <- views[length(views)]
            }
            if ("source" %in% names(evaluateit[[i]])) {
                if (evaluateit[[i]][["source"]]) {
                    evaluateit[[i]][["continue"]] <- FALSE
                }
            }
        }
        i <- i + 1
    }
    return(evaluateit)
}
getDatasets <- function() {
    attached <- data()$results[, -2]
    packages <- unique(attached[, "Package"])
    attached <- lapply(packages, function(x) {
        x <- attached[attached[, "Package"] == x, 2:3, drop = FALSE]
        x <- x[x[, 2] != "Internal Functions", , drop = FALSE] 
        titles <- as.list(x[, 2])
        names(titles) <- x[, 1]
        return(titles)
    })
    names(attached) <- packages
    return(attached)
}
getXasp <- function(x, type = "default") {
    clr <- range(x, na.rm = TRUE)
    pdf(templotfile)
    if (type == "default") {
        plot(seq(clr[1], clr[2], length.out = 100), seq(100), xlim = clr, type = "n", axes = FALSE)
        axis(1)
    }
    else if (type == "pretty") {
        plot(x, type = "n", axes = FALSE)
        axis(1, at = pretty(x))
    }
    else if (type == "margins") {
        plot(x, type = "n", axes = FALSE)
        axis(1, at = clr)
    }
    xasp <- par("xaxp")
    dev.off()
    file.remove(templotfile)
    return(xasp)
}
getXYplot <- function(foo) {
    X <- ev[[foo$dataset]][, foo$x]
    Y <- ev[[foo$dataset]][, foo$y]
    notX <- 1 - X
    notY <- 1 - Y
    rpofsuf <- list(pof(   X,    Y, rel = "suf"),
                    pof(notX,    Y, rel = "suf"),
                    pof(   X, notY, rel = "suf"),
                    pof(notX, notY, rel = "suf"))
    rpofsuf <- lapply(rpofsuf, function(x) {
        formatC(c(x$incl.cov$inclS, x$incl.cov$covS, x$incl.cov$PRI), format="f", digits = 3)
    })
    rpofnec <- list(pof(   X,    Y), 
                    pof(notX,    Y),
                    pof(   X, notY),
                    pof(notX, notY))
    rpofnec <- lapply(rpofnec, function(x) {
        formatC(c(x$incl.cov$inclN, x$incl.cov$covN, x$incl.cov$RoN), format="f", digits = 3)
    })
    return(list(rownames(ev[[foo$dataset]]), ev[[foo$dataset]][, foo$x], ev[[foo$dataset]][, foo$y], rpofsuf, rpofnec))
}
numhash <- function(x) {
    mean(as.integer(charToRaw(paste(capture.output(.Internal(inspect(x))), collapse = ""))))
}
warningstack <- NULL
ev <- new.env(parent = globalenv())
hashes <- list()
templotfile <- file.path(tempdir(), "plot.pdf")
visiblerows <- 17
visiblecols <- 8
pdf(templotfile)
emptyplot <- testplot <- recordPlot()
plotsize <- rep(5.729167, 2) 
sapply(dev.list(), dev.off)
file.remove(templotfile)
grafic <- emptyplot
testplot <- emptyplot
svgfile <- file.path(path.package("QCA"), "gui", "www", "css", "images", "plot.svg")
tempdata <- NULL
current_path <- getwd()
oktoset <- TRUE
filepath <- ""
extension <- ""
tcisdata <- TRUE
shinyServer(function(input, output, session) {
    session$onSessionEnded(stopApp)
    observe({
        dirfilist <- input$dirfilist
        session$sendCustomMessage(type = "dirfile", listFiles(current_path))
    })
    observe({ 
        read_table <- input$read_table
        filepath <<- ""
        oktoset <<- TRUE
        if (!is.null(input$dirfile_chosen)) {
            dfchosen <- input$dirfile_chosen
            oktoset <<- TRUE
            if (dfchosen[1] == "file") {
                filepath <<- file.path(gsub("[/]$", "", current_path), dfchosen[2])
            }
            else {
                splitpath <- unlist(strsplit(current_path, split=.Platform$file.sep))
                if (dfchosen[2] == ".." | dfchosen[2] == "...") {
                    if (length(splitpath) > 1) {
                        splitpath <- splitpath[-length(splitpath)]
                    }
                    if (identical(splitpath, "")) {
                        splitpath <- "/"
                    }
                    pathtobe <- paste(splitpath, collapse = .Platform$file.sep)
                    if (length(list.files(pathtobe)) > 0) {
                        current_path <<- pathtobe
                    }
                    else {
                        oktoset <<- FALSE
                    }
                }
                else {
                    current_dirs <- listFiles(current_path)$dirs
                    if (is.element(dfchosen[2], current_dirs)) {
                        pathtobe <- file.path(current_path, dfchosen[2])
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                        else {
                            oktoset <<- FALSE
                        }
                    }
                    else if (dfchosen[2] != "__stdir__") {
                        if (dfchosen[2] == "root") {
                            dfchosen[2] <- ""
                        }
                        splitpath <- splitpath[seq(which(splitpath == dfchosen[2]))]
                        pathtobe <- ifelse(length(splitpath) == 1,
                                           ifelse(identical(splitpath, ""), "/", splitpath),
                                           paste(splitpath, collapse = .Platform$file.sep))
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                        else {
                            oktoset <<- FALSE
                        }
                    }
                }
            }
            if (oktoset) {
                current_path <<- gsub("//", "/", current_path)
                if (!grepl("/", current_path)) {
                    current_path <<- paste(current_path, "/", sep = "")
                }
            }
            if (dfchosen[1] == "dir" & dfchosen[3] != "") {
                if (grepl("cannot change working directory", tryCatch(setwd(dfchosen[3]), error = function(e) e))) {
                    listtosend <- listFiles(current_path)
                    listtosend$filename <- "error!"
                    session$sendCustomMessage(type = "dirfile", listtosend)
                }
                else {
                    if (dfchosen[3] %in% listFiles(current_path)$dirs) {
                        pathtobe <- file.path(current_path, dfchosen[3])
                        if (length(list.files(pathtobe)) > 0) {
                            current_path <<- pathtobe
                        }
                    }
                    else {
                        current_path <<- dfchosen[3]
                    }
                    if (!grepl("/", current_path)) {
                        current_path <<- paste(current_path, "/", sep = "")
                    }
                    setwd(current_path)
                    session$sendCustomMessage(type = "dirfile", listFiles(current_path))
                }
            }
            else {
                setwd(current_path)
                session$sendCustomMessage(type = "dirfile", listFiles(current_path))
            }
        }
        else {
            current_path <<- getwd()
        }
        if (!identical(filepath, "")) {
            numevars <- ""
            header <- read_table$header
            colsep <- read_table$sep 
            row_names <- read_table$row_names
            decimal <- read_table$dec
            fileEncoding <- read_table$fileEncoding
            filename <- unlist(strsplit(basename(filepath), split="\\."))
            filename <- filename[-length(filename)]
            if (length(filename) > 1) {
                filename <- paste(filename, collapse=".")
            }
            if (!identical(row_names, "")) { 
                if (possibleNumeric(row_names)) {
                    row_names <- as.numeric(row_names)
                }
                tc <- capture.output(
                    tryCatch(
                        read.table(
                            filepath,
                            header = header,
                            ifelse(colsep == "tab", "\t", colsep),
                            row.names = row_names,
                            as.is = TRUE,
                            dec = decimal,
                            nrows = 2,
                            fileEncoding = fileEncoding
                        ),
                    error = function(e) e, warning = function(w) w)
                )
            }
            else {
                tc <- capture.output(
                    tryCatch(
                        read.table(
                            filepath,
                            header = header,
                            ifelse(colsep == "tab", "\t", colsep),
                            as.is = TRUE,
                            dec = decimal,
                            nrows = 2,
                            fileEncoding = fileEncoding
                        ),
                        error = function(e) e, warning = function(w) w)
                    )
            }
            if (any(grepl("subscript out of bounds", tc))) {
                mesaj <- paste("The data doesn't have ", row_names, " columns.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = mesaj, rownames="error!"))
                return(invisible())
            }
            else if (any(grepl("are not allowed", tc))) {
                mesaj <- paste("The row.names column has duplicated values.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = mesaj, rownames="error!"))
                return(invisible())
            }
            else if (any(grepl("data frame with 0 columns", tc))) {
                mesaj <- paste("The data has only 1 column.", sep = "")
                tc <- tryCatch(
                    read.table(
                        filepath,
                        header = header,
                        ifelse(colsep == "tab", "\t", colsep),
                        as.is = TRUE,
                        dec = decimal,
                        nrows = 2,
                        fileEncoding = fileEncoding
                    ),
                    error = function(e) e
                )
                session$sendCustomMessage(type = "tempdatainfo", list(ncols = 2, nrows = 2, colnames=c(colnames(tc), mesaj), rownames=""))
                return(invisible())
            }
            else if (any(grepl("attempt to select less than one element", tc))) {
                mesaj <- paste("The column \"", row_names, "\" was not found.", sep = "")
                session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = mesaj, rownames = "error!"))
                return(invisible())
            }
            tc <- tryCatch(
                read.table(
                    filepath,
                    header = header,
                    ifelse(colsep == "tab", "\t", colsep),
                    as.is = TRUE,
                    dec = decimal,
                    nrows = 2,
                    fileEncoding = fileEncoding
                ),
                error = function(e) e, warning = function(w) w
            )
            tcisdata <<- TRUE
            if (is.null(dim(tc))) {
                if (is.list(tc)) {
                    if (identical(names(tc), c("message", "call"))) {
                        tcisdata <<- FALSE
                        session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = tc$message, rownames = "error!"))
                    }
                }
            }
            else {
                if (grepl("X.PDF", names(tc)[1])) {
                    tcisdata <<- FALSE
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = "not a dataframe, this is a PDF file", rownames="error!"))
                }
            }
            if (tcisdata) {
                if (row_names != "") {
                    tc <- tryCatch(
                        read.table(
                            filepath,
                            header = header,
                            ifelse(colsep == "tab", "\t", colsep),
                            row.names = row_names,
                            as.is = TRUE,
                            dec = decimal,
                            fileEncoding = fileEncoding
                        ),
                        error = function(e) e, warning = function(w) w
                    )
                }
                else {
                    tc <- tryCatch(
                        read.table(
                            filepath,
                            header = header,
                            ifelse(colsep == "tab", "\t", colsep),
                            as.is = TRUE,
                            dec = decimal,
                            fileEncoding = fileEncoding
                        ),
                        error = function(e) e,
                        warning = function(w) w
                    )
                }
                if (identical(names(tc), c("message", "call"))) {
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols = 1, nrows = 1, colnames = tc$message, rownames = "error!"))
                }
                else {
                    tempdata <<- tc
                    cnames <- colnames(tempdata)
                    if (length(cnames) == 1) {
                        cnames = list(cnames)
                    }
                    rnames <- rownames(tempdata)
                    if (length(rnames) == 1) {
                        rnames = list(rnames)
                    }
                    session$sendCustomMessage(type = "tempdatainfo", list(ncols = ncol(tempdata),
                                                                     nrows = nrow(tempdata),
                                                                     colnames = cnames,
                                                                     rownames = rnames))
                }
            }
        }
    })
    observe({ 
        foo <- input$import
        if (!is.null(foo) & tcisdata) {
            result <- list(infobjs = NULL, console = NULL)
            if (foo$nameit) {
                ev[[foo$objname]] <- tempdata
            }
            else {
                result$console <- c(capture.output(tempdata), "")
            }
            result$infobjs <- infobjs(ev, ls(ev))
            session$sendCustomMessage(type = "fullinfo", result)
        }
    })
    observe({ 
        foo <- input$scrollobj
        if (!is.null(foo)) {
            scrollvh <- lapply(foo$scrollvh, function(x) unlist(x) + 1)
            visiblerows <<- foo$visiblerows + 1
            visiblecols <<- foo$visiblecols + 1
            if (!foo$alldata) {
                scrollvh <- scrollvh[foo$dataset]
            }
            tosend <- vector(mode = "list", length = length(scrollvh))
            names(tosend) <- names(scrollvh)
            for (n in names(scrollvh)) {
                nrowd <- nrow(ev[[n]])
                ncold <- ncol(ev[[n]])
                dscrollvh <- scrollvh[[n]]
                srow <- min(dscrollvh[1], nrowd - min(nrowd, visiblerows) + 1)
                scol <- min(dscrollvh[2], ncold - min(ncold, visiblecols) + 1)
                erow <- min(srow + visiblerows - 1, nrowd)
                ecol <- min(scol + visiblecols - 1, ncold)
                tosend[[n]] <- list(
                    theData = unname(as.list(ev[[n]][seq(srow, erow), seq(scol, ecol), drop = FALSE])),
                    dataCoords = paste(srow, scol, erow, ecol, ncold, sep="_"),
                    scrollvh = c(srow, scol) - 1
                )
            }
            session$sendCustomMessage(type = "scrollData", tosend)
        }
    })
    observe({ 
        foo <- input$thinfo
        if (!is.null(foo)) {
            response <- list()
            response$message <- "OK"
            if (foo$x != "") {
                foo$nth <- as.numeric(foo$nth)
                if (possibleNumeric((ev[[foo$dataset]][, foo$x]))) {
                    xasp <- getXasp(ev[[foo$dataset]][, foo$x])
                    response$vals <- unname(ev[[foo$dataset]][, foo$x])
                    response$thvals <- vector(length = 0)
                    response$fuzzyvals <- list()
                    response$prettyx <- seq(xasp[1], xasp[2], length.out = xasp[3] + 1)
                    if (foo$findth) {
                        response$thvals <- as.list(findTh(ev[[foo$dataset]][, foo$x], n = foo$nth))
                    }
                }
                else {
                    response$message <- "notnumeric"
                }
                session$sendCustomMessage(type = "dataPoints", response)
            }
        }
    })
    observe({ 
        foo <- input$calibrate
        if (!is.null(foo)) {
            checkit <- calibrateit(foo)
            if ("fuzzyvals" %in% names(checkit)) {
                session$sendCustomMessage(type = "fuzzyvals", checkit$fuzzyvals)
            }
            else { 
                session$sendCustomMessage(type = "calibrate", checkit$response)
            }
        }
    })
    observe({ 
        foo <- input$xyplot
        if (!is.null(foo)) {
            if (all(is.element(c(foo$x, foo$y), names(ev[[foo$dataset]])))) {
                session$sendCustomMessage(type = "xyplot", getXYplot(foo))
            }
        }
    })
    observe({ 
        foo <- input$Rcommand
        if (!is.null(foo)) {
            if (length(foo$plotsize) > 0) {
                 plotsize <<- unlist(foo$plotsize)
            }
            calib_before <- calib_after <- ""
            xyplot_before <- xyplot_after <- ""
            recalibrate <- FALSE
            if (length(foo$calibrate) > 0) {
                calib_before <- numhash(ev[[foo$calibrate$dataset]][, foo$calibrate$x, drop = FALSE])
                if (foo$calibrate$thsetter) {
                    recalibrate <- TRUE
                    tocalibrate <- foo$calibrate
                }
            }
            if (length(foo$xyplot) > 0) {
                xyplot_before <- numhash(ev[[foo$xyplot$dataset]][, c(foo$xyplot$x, foo$xyplot$y), drop = FALSE])
            }
            scrollvh <- lapply(foo$scrollvh, function(x) unlist(x) + 1)
            thinfo <- foo$thinfo
            hashes_before <- lapply(ev, function(x) {
                numhash(x)
            })
            if (length(dev.list()) > 0) {
                sapply(dev.list(), dev.off)
            }
            command <- unlist(strsplit(gsub("&nbsp;", " ", foo$command), split = "\n"))
            tosend <- list(evaluate = NULL, plot = FALSE, added = NULL, modified = NULL, deleted = NULL)
            rm(list = ls(envir = globalenv()), envir = globalenv())
            testplot <<- emptyplot 
            pdf(templotfile)
            dev.control("enable")
            if (identical(command, "warnings()")) {
                if (is.null(warningstack)) {
                    warnmsg <- "NULL"
                }
                else {
                    if (length(warningstack) == 1) {
                        warnmsg <- paste("Warning message: ", warningstack, sep = "\n")
                    }
                    else {
                        warnmsg <- paste(c("Warning messages: ", paste(seq(length(warningstack)), warningstack, sep = ": ")), collapse = "\n")
                    }
                }
                tosend$evaluate <- list(list(command = command, output = warnmsg, continue = FALSE))
            }
            else {
                tosend$evaluate <- evalparse(command)
            }
            lib <- unlist(lapply(tosend$evaluate, function(x) {
                if (any(names(x) == "library")) {
                    return(unlist(strsplit(x$library, split = "\\@\\$\\%\\$\\@")))
                }
            }))
            libcar <- "car" %in% lib
            if (!identical(lib, NULL)) {
                session$sendCustomMessage(type = "packages", getDatasets())
            }
            tosend$evaluate <- lapply(tosend$evaluate, function(x) {
                if (any(names(x) == "warning")) {
                    x$warning <- unlist(strsplit(x$warning, split = "\\@\\$\\%\\$\\@"))
                    if (length(x$warning) == 1) {
                        warningstack <<- unlist(x$warning)
                        x$warning <- list(paste("Warning message:", x$warning, sep = "\n"))
                    }
                    else {
                        if (length(x$warning) > 50) {
                            warningstack <<- unlist(x$warning[1:50])
                            x$warning <- list("There were 50 or more warnings (use warnings() to see the first 50)")
                        }
                        else {
                            warningstack <<- unlist(x$warning)
                            x$warning <- list(paste(c("Warning messages:", paste(seq(length(warningstack)), warningstack, sep = ": ")), collapse = "\n"))
                        }
                    }
                }
                return(x)
            })
            globjs <- ls(envir = globalenv())
            added <- NULL
            if (length(globjs) >= 0) {
                for (i in globjs) { 
                    if (any(names(ev) == i) & is.data.frame(get(i, globalenv()))) {
                        added <- c(added, i)
                    }
                    assign(i, get(i, globalenv()), ev)
                }
                rm(list = ls(envir = globalenv()), envir = globalenv())
            }
            if (length(dev.list()) > 0) {
                testplot <- recordPlot()
                sapply(dev.list(), dev.off)
            }
            if (!identical(testplot, emptyplot) & !libcar) {
                grafic <<- testplot
                tosend$plot <- TRUE
                svg(filename = svgfile, width = plotsize[1], height = plotsize[2])
                replayPlot(grafic)
                dev.off()
            }
            sapply(dev.list(), dev.off)
            if (exists("SOME__BIG__NAME", envir = ev)) {
                rm("SOME__BIG__NAME", envir = ev)
            }
            if (file.exists("Rplots.pdf")) {
                file.remove("Rplots.pdf")
            }
            hashes_after <- lapply(ev, function(x) {
                numhash(x)
            })
            if (length(foo$calibrate) > 0) {
                calib_after <- numhash(ev[[foo$calibrate$dataset]][, foo$calibrate$x, drop = FALSE])
            }
            if (length(foo$xyplot) > 0) {
                xyplot_after <- numhash(ev[[foo$xyplot$dataset]][, c(foo$xyplot$x, foo$xyplot$y), drop = FALSE])
            }
            added <- c(added, setdiff(names(hashes_after), names(hashes_before)))
            deleted <- setdiff(names(hashes_before), names(hashes_after))
            common <- intersect(names(hashes_before), names(hashes_after))
            modified <- names(hashes_before)[!is.element(hashes_before[common], hashes_after[common])]
            if (length(modified) > 0) {
                tosend$modified <- infobjs(ev, modified, scrollvh)
            }
            if (length(added) > 0) {
                tosend$added <- infobjs(ev, added)
            }
            if (length(deleted) > 0) {
                tosend$deleted <- as.list(deleted)
            }
            if (!identical(calib_before, calib_after)) {
                xasp <- getXasp(ev[[foo$calibrate$dataset]][, foo$calibrate$x])
                tosend$poinths <- list(dataset = foo$calibrate$dataset,
                                       x = foo$calibrate$x,
                                       rownames = rownames(ev[[foo$calibrate$dataset]]),
                                       vals = unname(as.vector(ev[[foo$calibrate$dataset]][, foo$calibrate$x])),
                                       fuzzyvals = list(),
                                       prettyx = seq(xasp[1], xasp[2], length.out = xasp[3] + 1))
                if (recalibrate) {
                    checkit <- calibrateit(tocalibrate)
                    if (is.element("fuzzyvals", names(checkit))) {
                        tosend$poinths$fuzzyvals <- checkit$fuzzyvals
                    }
                }
                if (foo$calibrate$findth) {
                    tosend$poinths$thvals <- findTh(ev[[foo$calibrate$dataset]][, foo$calibrate$x], n = foo$calibrate$nth)
                    if (length(tosend$poinths$thvals) == 1) {
                        tosend$poinths$thvals <- as.list(tosend$poinths$thvals)
                    }
                    tosend$poinths$message <- "OK"
                }
            }
            if (!identical(xyplot_before, xyplot_after)) {
                tosend$xyplot <- getXYplot(foo$xyplot)
            }
            session$sendCustomMessage(type = "Rcommand", tosend)
        }
    })
    observe({ 
        foo <- input$changes
        if (!is.null(foo)) {
            session$sendCustomMessage(type = "getChanges", readLines(system.file("ChangeLog", package="QCA")))
        }
    })
    observe({ 
        foo <- input$help
        if (!is.null(foo)) {
            browseURL(file.path(path.package("QCA"), "staticdocs", "index.html"))
        }
    })
    observe({ 
        foo <- input$pingobj
        if (!is.null(foo)) {
            session$sendCustomMessage(type = "ping", paste("bar", foo))
        }
    })
    observe({ 
        foo <- input$closeplot
        if (!is.null(foo)) {
            if (file.exists(svgfile)) {
                file.remove(svgfile)
            }
            grafic <<- emptyplot
            svg(filename = svgfile)
            replayPlot(grafic)
            sapply(dev.list(), dev.off)
        }
    })
    observe({ 
        foo <- input$saveRplot
        if (!is.null(foo)) {
            filename <- paste(foo$filename, foo$type, sep=".")
            if (foo$type == "png") {
                png(filename, width = plotsize[1]*96, height = plotsize[2]*96)
            }
            else if (foo$type == "bmp") {
                bmp(filename, width = plotsize[1]*96, height = plotsize[2]*96)
            }
            else if (foo$type == "jpeg") {
                jpeg(filename, width = plotsize[1]*96, height = plotsize[2]*96)
            }
            else if (foo$type == "tiff") {
                tiff(filename, width = plotsize[1]*96, height = plotsize[2]*96)
            }
            else if (foo$type == "svg") {
                svg(filename, width = plotsize[1], height = plotsize[2])
            }
            else if (foo$type == "pdf") {
                pdf(filename, width = plotsize[1], height = plotsize[2])
            }
            replayPlot(grafic)
            dev.off()
        }
    })
    observe({ 
        foo <- input$plotsize
        if (!is.null(foo)) {
            if (!identical(grafic, emptyplot)) {
                plotsize <<- foo
                svg(filename = svgfile, width = plotsize[1], height = plotsize[2])
                replayPlot(grafic)
                dev.off()
                session$sendCustomMessage(type = "resizePlot", TRUE)
            }
            else {
                session$sendCustomMessage(type = "resizePlot", FALSE)
            }
        }
    })
    observe({ 
        foo <- input$quit
        if (!is.null(foo)) {
            if (file.exists(svgfile)) {
                file.remove(svgfile)
            }
            session$sendCustomMessage(type = "okquit", list(1))
        }
    })
    session$sendCustomMessage(type = "fullinfo", list(infobjs = infobjs(ev, ls(ev)), console = NULL))
    session$sendCustomMessage(type = "packages", getDatasets())
})
