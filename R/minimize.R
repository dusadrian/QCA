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

`minimize` <- function(
    input, include = "", dir.exp = NULL, details = FALSE, pi.cons = 0,
    sol.cons = 0, all.sol = FALSE, row.dom = FALSE, first.min = FALSE,
    max.comb = 0, use.labels = FALSE, method = "CCubes", ...
) {
    metacall <- match.call()
    dots <- as.list(substitute(list(...)))[-1]
    if (isTRUE(dots$categorical)) { 
        use.labels <- TRUE
        dots$categorical <- NULL
    }
    if (is.element("data", names(dots))) {
        input <- eval.parent(dots$data)
        dots$data <- NULL
    }
    min.pin <- isTRUE(dots$min.pin)
    enter <- if (is.element("enter", names(dots))) dots$enter else "\n" 
    if (missing(input)) {
        admisc::stopError(
            "The input (a truth table or a dataset) is missing.",
            ... = ...
        )
    }
    else {
        ttinput <- methods::is(input, "QCA_tt")
        if (is.matrix(input)) {
            if (is.null(colnames(input))) {
                admisc::stopError(
                    "The data should have column names.",
                    ... = ...
                )
            }
            if (any(duplicated(rownames(input)))) {
                rownames(input) <- seq(nrow(input))
            }
            input <- as.data.frame(input)
            for (i in seq(ncol(input))) {
                if (admisc::possibleNumeric(input[, i])) {
                    input[, i] <- admisc::asNumeric(input[, i])
                }
            }
        }
        if(!(is.data.frame(input) | ttinput)) {
            admisc::stopError(
                "The input should be a truth table or a dataset.",
                ... = ...
            )
        }
    }
    print.truth.table <- details & !ttinput
    if (ttinput) {
        nms <- colnames(input$recoded.data)[seq(length(input$noflevels))]
    }
    else {
        nms <- colnames(input)
    }
    if (length(dots) > 0) {
        for (i in seq(length(dots))) {
            tc <- admisc::tryCatchWEM(
                ieval <- eval(dots[[i]])
            )
            if (!is.null(tc$error)) {
                ieval <- admisc::recreate(dots[[i]], nms)
            }
            dots[[i]] <- ieval
        }
    }
    back.args <- c(
        "outcome", "conditions", "n.cut", "incl.cut", "complete", "show.cases",
        "sort.by" = "", "use.letters", "inf.test", "rowdom", "direxp", "neg.out",
        "data", "relation", "explain", "omit", "exclude"
    )
    check.args <- pmatch(names(dots), back.args)
    names(dots)[!is.na(check.args)] <- back.args[check.args[!is.na(check.args)]]
    explain     <- if (is.element("explain",     names(dots))) dots$explain      else "1"
    outcome     <- if (is.element("outcome",     names(dots))) dots$outcome      else ""
    conditions  <- if (is.element("conditions",  names(dots))) dots$conditions   else ""
    incl.cut    <- if (is.element("incl.cut",    names(dots))) dots$incl.cut     else 1
    n.cut       <- if (is.element("n.cut",       names(dots))) dots$n.cut        else 1
    complete    <- if (is.element("complete",    names(dots))) dots$complete     else FALSE
    show.cases  <- if (is.element("show.cases",  names(dots))) dots$show.cases   else FALSE
    dcc         <- if (is.element("dcc",         names(dots))) dots$dcc          else FALSE
    sort.by     <- if (is.element("sort.by",     names(dots))) dots$sort.by      else ""
    use.letters <- if (is.element("use.letters", names(dots))) dots$use.letters  else FALSE
    inf.test    <- if (is.element("inf.test",    names(dots))) dots$inf.test     else ""
    relation    <- if (is.element("relation",    names(dots))) dots$relation     else "sufficiency"
    neg.out     <- if (is.element("neg.out",     names(dots))) dots$neg.out      else FALSE
    pi.depth    <- if (is.element("pi.depth",    names(dots))) dots$pi.depth     else 0
    sol.cov     <- if (is.element("sol.cov",     names(dots))) dots$sol.cov      else 1
    sol.depth   <- if (is.element("sol.depth",   names(dots))) dots$sol.depth    else 0
    keep.trying <- if (is.element("keep.trying", names(dots))) dots$keep.trying  else FALSE
    if (is.element("omit", names(dots)) && !is.element("exclude", names(dots))) {
        dots$exclude <- dots$omit
    }
    if (is.null(include)) {
        admisc::stopError(
            "The <include> argument cannot be NULL.",
            ... = ...
        )
    }
    dir.exp <- admisc::recreate(substitute(dir.exp))
    row.dom     <- if (is.element("rowdom",      names(dots))) dots$rowdom       else row.dom
    dir.exp     <- if (is.element("direxp",      names(dots))) dots$direxp       else dir.exp
    if (identical(dir.exp, character(0))) {
        dir.exp <- NULL
    }
    else if (!is.null(dir.exp)) {
        if (length(dir.exp) == 1) {
            if (grepl(":", dir.exp)) {
                des <- unlist(strsplit(dir.exp, split = ":"))
                if (!all(is.element(des, nms))) {
                    admisc::stopError(
                        "Inexisting condition(s) in the sequence of directional expectations.",
                        ... = ...
                    )
                }
                dir.exp <- nms[seq(which(nms == des[1]), which(nms == des[2]))]
            }
            else {
                dir.exp <- admisc::splitstr(dir.exp)
            }
        }
    }
    if (identical(include, "")) {
        if (!is.null(dir.exp)) {
            admisc::stopError(
                "Directional expectations cannot be specified without including the remainders.",
                ... = ...
            )
        }
    }
    if (is.character(explain) & !identical(explain, "1")) {
        explain <- admisc::splitstr(explain)
    }
    if (is.character(include) & !identical(include, "")) {
        include <- admisc::splitstr(include)
    }
    if (ttinput) { 
        tt <- input
        ttargs <- setdiff(names(formals(truthTable)), c("show.cases", "use.labels"))
        if (any(is.element(ttargs, names(dots)))) {
            callist <- as.list(tt$call)
            common <- intersect(names(dots), ttargs)
            if (length(common) > 0) {
                for (i in seq(length(common))) {
                    callist[[common[i]]] <- dots[[common[i]]]
                }
            }
            dataname <- callist$data
            for (i in seq(2, length(callist))) {
                if (length(callist[[i]]) == 1) {
                    callist[[i]] <- admisc::recreate(callist[[i]])
                }
            }
            callist$data <- tt$initial.data
            tt <- do.call("truthTable", callist[-1])
            callist$data <- dataname
            tt$call <- as.call(callist)
        }
        if (isTRUE(use.labels)) {
            tt$options$use.labels <- TRUE
        }
        else if (isTRUE(tt$options$use.labels)) {
            use.labels <- TRUE
        }
    }
    else {
        if (identical(outcome, "")) {
            admisc::stopError(
                "Consider creating a truth table first, or formally specify the argument <outcome>.",
                ... = ...
            )
        }
        if (any(c(pi.cons, sol.cons) > 0) & incl.cut[1] == 1) {
            incl.cut[1] <- min(c(pi.cons, sol.cons))
        }
        if (is.character(outcome) & !identical(outcome, "")) {
            outcome <- admisc::splitstr(outcome)
        }
        if (length(outcome) > 1) {
            return(do.call("minimizeLoop", as.list(metacall)[-1], envir = parent.frame()))
        }
        outcome.copy <- outcome.name <- outcome
        indata <- input 
        testoutcome <- admisc::tryCatchWEM(
            trout <- admisc::translate(outcome, data = input)
        )
        if (is.element("error", names(testoutcome))) {
            admisc::stopError(
                "Incorrect outcome specification.",
                ... = ...
            )
        }
        testrout <- apply(trout, 2, function(x) {
            all(x != "-1")
        })
        if (sum(testrout) == 1) {
            outcome.name <- names(testrout)[testrout]
        }
        if (identical(conditions, "")) {
            conditions <- names(input)[-which(names(input) == outcome.name)]
        }
        else {
            conditions <- admisc::splitstr(conditions)
        }
        verify.data(input, outcome.name, conditions)
        if (length(conditions) == 1) {
            if (grepl(":", conditions)) {
                nms <- colnames(data)
                cs <- unlist(strsplit(conditions, split = ":"))
                conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
            }
        }
        input <- input[, unique(c(conditions, names(testrout)[testrout]))]
        verify.minimize(input, outcome.name, conditions, explain, include, use.letters)
        if (!is.element("incl.cut", names(dots))) {
            dots$incl.cut <- incl.cut
        }
        tt <- do.call("truthTable", c(list(data = input), dots))
    }
    curly <- grepl("\\{", tt$options$outcome)
    recdata <- tt$recoded.data
    conditions <- colnames(recdata)[seq(length(tt$noflevels))]
    outcome <- colnames(recdata)[ncol(recdata)]
    indata <- tt$initial.data[, match(colnames(recdata), colnames(tt$initial.data)), drop = FALSE]
    use.letters <- tt$options$use.letters
    show.cases <- show.cases | tt$options$show.cases 
    neg.out <- tt$options$neg.out
    output <- list()
    output$tt <- tt
    output$options$print.truth.table <- print.truth.table
    rowsNotMissing <- which(tt$tt$OUT != "?")
    if (any(tt$tt$OUT == "?")) {
        missings <- which(tt$tt$OUT == "?")
        tt$tt <- tt$tt[-missings, , drop = FALSE]
    }
    ttrownms <- admisc::asNumeric(rownames(tt$tt))
    noflevels <- tt$noflevels
    mbase <- as.integer(rev(c(1, cumprod(rev(noflevels))))[-1])
    mbaseplus <- rev(c(1, cumprod(rev(noflevels + 1))))[-1]
    alreadyletters <- sum(nchar(colnames(recdata)[-ncol(recdata)])) == ncol(recdata) - 1
    tt$tt[seq(length(conditions))] <- as.data.frame(lapply(tt$tt[seq(length(conditions))], function(x) {
        x[is.element(x, c("-", "dc"))] <- -1
        return(admisc::asNumeric(x))
    }))
    pos.incl <- unique(c(explain, include)) 
    subset.tt <- ttrownms[is.element(tt$tt[, "OUT"], explain)]
    subset.pos <- is.element(tt$tt[, "OUT"], pos.incl)
    pos.matrix <- as.matrix(tt$tt[subset.pos, seq(length(noflevels))])
    rownames(pos.matrix) <- drop(pos.matrix %*% mbase) + 1
    pos.matrix <- pos.matrix + 1
    neg.matrix <- as.matrix(tt$tt[!is.element(tt$tt[, "OUT"], pos.incl), seq(length(noflevels))])
    neg.matrix <- matrix(as.numeric(neg.matrix), ncol = length(noflevels)) + 1
    rownames(neg.matrix) <- drop((neg.matrix - 1) %*% mbase) + 1
    if (sum(subset.pos) == 0) {
        admisc::stopError(
            "None of the values in OUT is explained. Please check the truth table.",
            ... = ...
        )
    }
    inputt <- as.matrix(tt$tt[is.element(ttrownms, subset.tt), seq(length(noflevels)), drop = FALSE])
    rownames(inputt) <- drop(inputt %*% mbase) + 1
    inputt <- inputt + 1
    inputcases <- tt$cases[is.element(tt$indexes, subset.tt)]
    nofcases1 <- sum(tt$tt$n[tt$tt$OUT == 1])
    nofcases0 <- sum(tt$tt$n[tt$tt$OUT == 0])
    nofcasesC <- sum(tt$tt$n[tt$tt$OUT == "C"])
    excl.matrix <- matrix(nrow = 0, ncol = length(conditions))
    output$negatives <- sort(drop((neg.matrix - 1) %*% mbase) + 1)
    rownms <- rownames(inputt)
    if (nrow(pos.matrix) == 0) {
        admisc::stopError(
            "Nothing to explain. Please check the truth table.",
            ... = ...
        )
    }
    include <- admisc::trimstr(include)
    incl.rem <- is.element("?", include)
    if (
        nrow(neg.matrix) == 0 &
        incl.rem &
        !is.element("causalChain", names(dots))
    ) {
        admisc::stopError(
            paste(
                "All truth table configurations are used, all conditions are minimized.",
                "       Please check the truth table.",
                sep = "\n"
            ),
            ... = ...
        )
    }
    expressions <- pos.matrix
    recdata[, conditions] <- as.data.frame(lapply(recdata[, conditions, drop = FALSE], function(x) {
        x[is.element(x, c("-", "?", "dc"))] <- -1
        return(as.numeric(x))
    }))
    mv <- any(recdata[, seq(ncol(recdata) - 1)] > 1) | tt$multivalue
    collapse <- "*"
    changed <- FALSE
    if (use.letters & !alreadyletters) {
        colnms <- LETTERS[seq(ncol(inputt))]
        changed <- TRUE
    }
    else {
        colnms <- colnames(recdata[, seq(ncol(inputt)), drop = FALSE])
    }
    admisc::setColnames(expressions, colnms)
    admisc::setColnames(inputt, colnms)
    admisc::setColnames(pos.matrix, colnms)
    admisc::setColnames(neg.matrix, colnms)
    rownames(neg.matrix) <- (neg.matrix - 1) %*% mbase + 1
    output$initials <- admisc::writePrimeimp(
        impmat = inputt,
        mv = mv,
        collapse = collapse,
        curly = curly
    )
    expressions <- .Call("C_QMC", expressions, noflevels, PACKAGE = "QCA")
    if (is.element("simplify", names(dots))) {
        expressions <- admisc::sortExpressions(expressions)
    }
    callist <- list(
        expressions = expressions,
        mv = mv,
        collapse = collapse,
        inputt = inputt,
        row.dom = row.dom,
        initial = rownms,
        all.sol = all.sol,
        indata = indata,
        curly = curly,
        use.labels = tt$options$use.labels,
        enter = enter
    )
    callist <- c(callist, dots)
    if (!incl.rem || (!is.null(dir.exp) & !identical(include, ""))) {
        if (pi.cons > 0) {
            callist$outcome <- outcome
        }
        c.sol <- p.sol <- do.call("getSolution", callist)
    }
    if (incl.rem) {
        pos.matrix <- inputt
        if (method == "QMC") {
            expressions <- .Call("C_QMC", createMatrix(noflevels)[-output$negatives, , drop = FALSE] + 1, noflevels, PACKAGE = "QCA")
            admisc::setColnames(expressions, colnames(inputt))
        }
        else if (method == "eQMC") {
            if (nrow(neg.matrix) > 0) {
                expressions <- sort.int(setdiff(findSupersets(pos.matrix, noflevels + 1), findSupersets(neg.matrix, noflevels + 1)))
            }
            else {
                expressions <- sort.int(findSupersets(pos.matrix, noflevels + 1))
            }
            expressions <- .Call("C_removeRedundants", expressions, noflevels, mbaseplus, PACKAGE = "QCA")
            expressions <- admisc::sortExpressions(getRow(expressions, noflevels + 1))
            admisc::setColnames(expressions, colnames(inputt))
        }
        else { 
            extended.data <- as.matrix(tt$recoded.data)
            if (nrow(excl.matrix) > 0) {
                extended.data <- rbind(extended.data, cbind(excl.matrix, 0))
            }
            if (sol.cons > 0 & all.sol & sol.depth == 0) {
                sol.depth <- 7
            }
            expressions <- .Call("C_Cubes", list(
                            tt = cbind(rbind(pos.matrix, neg.matrix) - 1, rep(c(1, 0), c(nrow(pos.matrix), nrow(neg.matrix)))),
                            data = extended.data,
                            all.sol = all.sol,
                            row.dom = row.dom,
                            min.pin = min.pin,
                            pi.cons = pi.cons,
                            depth = as.integer(c(pi.depth, sol.depth)),
                            sol.cons = sol.cons,
                            sol.cov = sol.cov,
                            fs = tt$fs,
                            max.comb = max.comb,
                            first.min = first.min,
                            keep.trying = keep.trying,
                            gurobi = !isFALSE(dots$gurobi) && eval(parse(text = "requireNamespace('gurobi', quietly = TRUE)")),
                            solind = !isFALSE(dots$solind)),
                            PACKAGE = "QCA")
        }
        callist$expressions <- expressions
        p.sol <- do.call("getSolution", callist)
    }
    output$PIchart <- p.sol$mtrx
    class(output$PIchart) <- c("matrix", "QCA_pic")
    attr(output$PIchart, "PI") <- p.sol$expressions
    output$primes    <- p.sol$reduced$expressions
    output$solution  <- p.sol$solution.list[[1]]
    output$essential <- p.sol$solution.list[[2]]
    output$options$explain     <- explain
    output$options$include     <- include
    output$options$neg.out     <- neg.out
    output$options$details     <- details
    output$options$sol.cons    <- sol.cons
    output$options$sol.cov     <- sol.cov
    output$options$relation    <- relation
    output$options$show.cases  <- show.cases
    output$options$use.letters <- use.letters
    output$options$collapse    <- collapse
    output$options$curly       <- curly
    output$options$use.labels  <- use.labels
    expr.cases <- rep(NA, nrow(p.sol$reduced$expressions))
    tt.rows <- admisc::writePrimeimp(
        impmat = inputt,
        mv = mv,
        collapse = collapse
    )
    if (any(grepl("[*]", rownames(p.sol$reduced$expressions)))) {
        if (use.letters) {
            mtrxlines <- makeChart(
                primes = rownames(p.sol$reduced$expressions),
                snames = LETTERS[seq(length(conditions))],
                configs = tt.rows,
                mv = mv,
                noflevels = noflevels
            )
        }
        else {
            mtrxlines <- makeChart(
                primes = rownames(p.sol$reduced$expressions),
                snames = conditions,
                configs = tt.rows,
                mv = mv,
                noflevels = noflevels
            )
        }
    }
    else {
        if (use.letters) {
            mtrxlines <- makeChart(
                primes = rownames(p.sol$reduced$expressions),
                configs = tt.rows,
                snames = LETTERS[seq(length(conditions))],
                mv = mv,
                noflevels = noflevels
            )
        }
        else {
            mtrxlines <- makeChart(
                primes = rownames(p.sol$reduced$expressions),
                configs = tt.rows,
                snames = conditions,
                mv = mv,
                noflevels = noflevels
            )
        }
    }
    colnames(mtrxlines) <- colnames(p.sol$reduced$mtrx) <- rownms
    for (l in seq(length(expr.cases))) {
        expr.cases[l] <- paste(inputcases[mtrxlines[l, ]], collapse="; ")
    }
    output$inputcases <- inputcases
    conds <- conditions
    if (all(is.element(colnames(p.sol$reduced$expressions), conditions))) {
        conds <- colnames(p.sol$reduced$expressions)
    }
    else {
        if (all(is.element(colnames(p.sol$reduced$expressions), LETTERS))) {
            conds <- conditions[match(colnames(p.sol$reduced$expressions), LETTERS)]
        }
    }
    if (!is.element("simplify", names(dots))) {
        poflist <- list(
            setms = paste(
                rownames(p.sol$reduced$expressions),
                collapse = "+"
            ),
            outcome = tt$options$outcome,
            data = indata,
            relation = "sufficiency",
            use.labels = use.labels,
            neg.out = neg.out,
            minimize = TRUE,
            use.letters = tt$options$use.letters,
            show.cases = TRUE,
            cases = expr.cases
        )
        if (poflist$use.letters) {
            names(poflist$data)[seq(length(conditions))] <- LETTERS[seq(length(conditions))]
        }
        if (length(output$solution) > 1) {
            poflist$solution.list <- output$solution
            poflist$essential <- output$essential
        }
        poflist$use.labels <- output$tt$use.labels
        listIC <- do.call("pof", poflist)
        if (sol.cons > 0 & identical(include, "")) {
            error <- FALSE
            if (is.element("overall", names(listIC))) {
                lind <- length(listIC$individual)
                eligible <- logical(lind)
                for (i in seq(lind)) {
                    eligible[i] <- listIC$individual[[i]]$sol.incl.cov[1, 1] >= sol.cons
                }
                wel <- which(eligible)
                if (length(wel) == 0) {
                    error <- TRUE
                }
                else if (length(wel) == 1) {
                    listIC <- list(
                        incl.cov = listIC$individual[wel]$incl.cov,
                        pims = listIC$individual[wel]$pims,
                        sol.incl.cov = listIC$individual[wel]$sol.incl.cov,
                        options = listIC$options
                    )
                }
                else {
                    listIC$individual <- listIC$individual[wel]
                }
            }
            else {
                error <- listIC$sol.incl.cov[1, 1] < sol.cons
            }
            if (error) {
                admisc::stopError(
                    "There are no solutions, given these constraints.",
                    ... = ...
                )
            }
        }
        listIC$options$show.cases <- show.cases
        output$pims <- listIC$pims
        attr(output$pims, "conditions") <- conditions
        output$IC <- listIC
    }
    output$numbers <- c(
        OUT1 = nofcases1,
        OUT0 = nofcases0,
        OUTC = nofcasesC,
        Total = nofcases1 + nofcases0 + nofcasesC
    )
    mtrx <- p.sol$mtrx[p.sol$all.PIs, , drop = FALSE]
    SA <- TRUE
    if (is.element("SA", names(dots))) {
        SA <- dots$SA
    }
    if (SA & 3^length(noflevels) < 2^31) {
        conds <- conditions
        if (tt$options$use.letters) {
            conds <- LETTERS[seq(length(conditions))]
        }
        mbaseexpr <- rev(c(1, cumprod(rev(noflevels[is.element(conds, colnames(p.sol$reduced$expressions))] + 1))))[-1]
        output$SA <- lapply(p.sol$solution.list[[1]], function(x) {
            p.expressions <- p.sol$reduced$expressions[x, , drop = FALSE]
            temp <- apply(p.expressions, 1, function(pr) {
                indices <- rev(which(pr == 0))
                tempr <- NULL
                for (k in indices) {
                    if (is.null(tempr)) {
                        tempr <- drop(mbaseexpr %*% pr) + sum(mbaseexpr[pr == 0])
                        temp2 <- tempr
                    }
                    for (lev in seq(noflevels[k] - 1)) {
                        temp2 <- c(temp2, tempr + mbaseexpr[k]*lev)
                    }
                    tempr <- temp2
                }
                return(tempr)
            })
            if (all(is.null(temp))) {
                SAx <- matrix(nrow = 0, ncol = ncol(inputt))
            }
            else {
                temp <- sort(unique(as.vector(unlist(temp))))
                temp <- temp[!is.element(temp, drop(inputt %*% mbaseplus))]
                if (length(temp) == 0) {
                    SAx <- matrix(nrow = 0, ncol = ncol(inputt))
                }
                else {
                    SAx <- getRow(temp + 1, noflevels + 1) - 1L
                    rownames(SAx) <- drop(SAx %*% mbase) + 1
                }
            }
            colnames(SAx) <- colnames(inputt)
            return(SAx)
        })
        prettyNums <- formatC(
            seq(length(p.sol$solution.list[[1]])),
            digits = nchar(length(p.sol$solution.list[[1]])) - 1,
            flag = 0
        )
        if (!is.null(dir.exp) & !identical(include, "")) {
            if (!identical(c.sol$solution.list, NA)) {
            dir.exp <- verify.dir.exp(recdata, outcome, conditions, noflevels, dir.exp, enter)
            isoliscsol <- identical(dir.exp, matrix(0L, ncol = length(conditions)))
            if (all(unlist(lapply(output$SA, is.null))) | isoliscsol) {
                ECmat <- as.data.frame(matrix(ncol = length(conditions), nrow = 0))
                colnames(ECmat) <- colnames(inputt)
            }
            else {
                result <- .Call(
                    "C_getEC",
                    dir.exp,
                    c.sol$expressions,
                    c.sol$sol.matrix,
                    p.sol$expressions,
                    p.sol$sol.matrix,
                    output$SA,
                    as.integer(noflevels),
                    PACKAGE = "QCA"
                )
                EClist <- result[[1]]
                isols <- result[[2]]
                intsel <- result[[3]]
            }
            tt.rows <- admisc::writePrimeimp(
                impmat = inputt,
                mv = mv,
                collapse = collapse
            )
            i.sol <- vector("list", ncol(c.sol$sol.matrix)*ncol(p.sol$sol.matrix))
            index <- 1
            for (c.s in seq(ncol(c.sol$sol.matrix))) {
                for (p.s in seq(ncol(p.sol$sol.matrix))) {
                    names(i.sol)[index] <- paste("C", c.s, "P", p.s, sep = "")
                    if (all(unlist(lapply(output$SA, is.null))) | isoliscsol) {
                        i.sol[[index]]$EC <- ECmat
                        i.sol[[index]]$DC <- ECmat
                        i.sol.index <- c.sol
                        i.sol.index$solution.list[[1]] <- c.sol$solution.list[[1]][[c.s]]
                    }
                    else {
                        i.sol[[index]]$EC <- EClist[[index]]
                        i.sol[[index]]$DC <- output$SA[[p.s]][setdiff(rownames(output$SA[[p.s]]), rownames(EClist[[index]])), , drop = FALSE]
                        pos.matrix.i.sol <- unique(as.matrix(rbind(pos.matrix, EClist[[index]] + 1L)))
                        expressions <- .Call(
                            "C_QMC",
                            pos.matrix.i.sol,
                            noflevels,
                            PACKAGE = "QCA"
                        )
                        callist$expressions <- expressions
                        i.sol.index <- do.call("getSolution", callist)
                        i.sol.index$expressions <- i.sol.index$expressions[rowSums(i.sol.index$mtrx) > 0, , drop = FALSE]
                    }
                    i.sol[[index]]$solution       <- i.sol.index$solution.list[[1]]
                    i.sol[[index]]$essential      <- i.sol.index$solution.list[[2]]
                    i.sol[[index]]$primes         <- i.sol.index$reduced$expressions
                    i.sol[[index]]$PIchart        <- i.sol.index$mtrx
                    class(i.sol[[index]]$PIchart) <- c("matrix", "QCA_pic")
                    i.sol[[index]]$c.sol          <- c.sol$solution.list[[1]][[c.s]]
                    i.sol[[index]]$p.sol          <- p.sol$solution.list[[1]][[p.s]]
                    expr.cases <- rep(NA, nrow(i.sol.index$reduced$expressions))
                    if (use.letters) {
                        mtrxlines <- makeChart(
                            primes = rownames(i.sol.index$reduced$expressions),
                            configs = tt.rows,
                            snames = LETTERS[seq(length(conditions))],
                            mv = mv,
                            noflevels = noflevels
                        )
                    }
                    else {
                        mtrxlines <- makeChart(
                            primes = rownames(i.sol.index$reduced$expressions),
                            configs = tt.rows,
                            snames = conditions,
                            mv = mv,
                            noflevels = noflevels
                        )
                    }
                    for (l in seq(length(expr.cases))) {
                        expr.cases[l] <- paste(inputcases[which(mtrxlines[l, ])], collapse="; ")
                    }
                    poflist <- list(
                        setms = paste(
                            rownames(i.sol.index$reduced$expressions),
                            collapse = "+"
                        ),
                        outcome = tt$options$outcome,
                        data = indata,
                        relation = "sufficiency",
                        use.labels = use.labels,
                        neg.out = neg.out,
                        minimize = TRUE,
                        use.letters = tt$options$use.letters,
                        show.cases = TRUE,
                        cases = expr.cases
                    )
                    if (length(i.sol.index$solution.list[[1]]) > 1) {
                        poflist$solution.list <- i.sol.index$solution.list[[1]]
                        if (!identical(i.sol.index$solution.list[[1]], i.sol.index$solution.list[[2]])) {
                            poflist$essential <- i.sol.index$solution.list[[2]]
                        }
                    }
                    poflist$use.labels <- output$tt$use.labels
                    i.sol[[index]]$IC <- do.call("pof", poflist)
                    i.sol[[index]]$IC$options$show.cases <- show.cases
                    i.sol[[index]]$pims <- i.sol[[index]]$IC$pims
                    attr(i.sol[[index]]$pims, "conditions") <- conditions
                    i.sol[[index]]$IC$pims <- NULL
                    index <- index + 1
                }
            }
            output$i.sol <- i.sol
        }}
        names(output$SA) <- paste("M", prettyNums, sep = "")
        output$SA <- lapply(output$SA, as.data.frame)
    }
    if (any(names(output) == "i.sol")) {
        for (i in seq(length(output$i.sol))) {
            output$i.sol[[i]]$EC <- as.data.frame(output$i.sol[[i]]$EC)
        }
    }
    if (!methods::is(input, "QCA_tt")) {
        output$tt$options$outcome <- outcome.copy
    }
    output$complex <- p.sol$complex
    output$call <- metacall
    if (is.element("via.web", names(dots))) {
        output$via.web <- TRUE
    }
    if (mv & !grepl(mvregexp, output$tt$options$outcome)) {
        output$tt$options$outcome <- paste(output$tt$options$outcome, "[1]", sep = "")
    }
    return(structure(output, class = "QCA_min"))
}
`minimizeLoop` <-
function(...) {
    allargs <- list(...)
    verify.mqca(allargs)
    outcome <- admisc::splitstr(allargs$outcome)
    minimize.list <- lapply(outcome, function(x) {
        allargs[["outcome"]] <- x
        return(do.call("minimize", allargs))
    })
    names(minimize.list) <- outcome
    return(structure(minimize.list, class = "QCA_loopmin"))
}
`eqmcc` <- function(...) {
    .Deprecated(msg = "Function eqmcc() is deprecated, and has been renamed to minimize()\n")
    minimize(...)
}
