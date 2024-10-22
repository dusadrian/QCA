`truthTable` <- function(
    data, outcome = "", conditions = "", incl.cut = 1, n.cut = 1, pri.cut = 0,
    exclude = NULL, complete = FALSE, use.letters = FALSE, use.labels = FALSE,
    show.cases = FALSE, dcc = FALSE, sort.by = "", inf.test = "", ...
) {

    metacall <- match.call(expand.dots = TRUE)
    dots <- admisc::recreate(substitute(list(...)))

    back.args <- c("outcome", "conditions", "n.cut", "incl.cut", "complete",
                    "show.cases", "sort.by", "use.letters", "inf.test")
    check.args <- pmatch(names(dots), back.args)
    names(dots)[!is.na(check.args)] <- back.args[check.args[!is.na(check.args)]]

    ic0 <- 1

    if (is.character(incl.cut) & length(incl.cut) == 1) {
        incl.cut <- admisc::splitstr(incl.cut)
    }

    ic1 <- incl.cut[1]
    if (length(incl.cut) > 1) {
        ic0 <- incl.cut[2]
    }

    ###
    ### ### backwards compatibility
    ###
        neg.out <- isTRUE(dots$neg.out)

        if (is.element("incl.cut1", names(dots)) && identical(ic1, 1)) {
            ic1 <- dots$incl.cut1
            incl.cut[1] <- ic1
        }

        if (is.element("incl.cut0", names(dots)) && identical(ic0, 1)) {
            ic0 <- dots$incl.cut0
            incl.cut[2] <- ic0
        }

        if (isTRUE(dots$categorical)) { # backwards compatibility
            use.labels <- TRUE
            dots$categorical <- NULL
        }
    ###
    ### ### backwards compatibility
    ###

    initialcols <- colnames(data)
    outcome <- admisc::recreate(substitute(outcome), colnames(data))

    if (length(admisc::splitstr(outcome)) > 1) {
        admisc::stopError(
            "Only one outcome is allowed.",
            ... = ...
        )
    }

    # just in case it was a tibble, produced by any tidyverse package
    data <- as.data.frame(data)

    outcome.copy <- outcome
    initial.data <- data
    declared <- FALSE # just to initialise
    multivalue <- grepl(mvregexp, outcome)


    if (!identical(outcome, "")) {
        outcometest <- admisc::tryCatchWEM(
            admisc::translate(admisc::notilde(outcome), data = data)
        )
        if (!is.null(outcometest$error)) {
            if (
                grepl("does not match the set names", outcometest$error)
            ) {
                admisc::stopError(
                    "Incorrect outcome specification.",
                    ... = ...
                )
            }
            else {
                admisc::stopError(
                    gsub("\\n", "", outcometest$error),
                    ... = ...
                )
            }
        }

        if (grepl("\\+|\\*", outcome)) {
            initial.data[, outcome] <- data[, outcome] <- admisc::compute(outcome, data)
        }
        else if (multivalue) {
            if (any(grepl("\\{", outcome))) {
                outcome.value <- admisc::curlyBrackets(outcome)
                outcome <- admisc::curlyBrackets(outcome, outside = TRUE)
            }
            else {
                outcome.value <- admisc::squareBrackets(outcome)
                outcome <- admisc::squareBrackets(outcome, outside = TRUE)
            }

            outcometest <- data[, admisc::notilde(outcome)]
            declared <- inherits(outcometest, "declared")
            if (declared) {
                outcomelabels <- attr(outcometest, "labels", exact = TRUE)
            }

            data[, admisc::notilde(outcome)] <- is.element(
                data[, admisc::notilde(outcome)],
                admisc::splitstr(outcome.value)
            ) * 1
        }
        else {
            outcometest <- data[, admisc::notilde(outcome)]
            declared <- inherits(outcometest, "declared")
            if (declared) {
                outcomelabels <- attr(outcometest, "labels", exact = TRUE)
            }
            ### TODO: vezi ce se intampla cand are tilde, deci ceea ce se explica este absenta si atunci labels pentru valoarea 0 trebuie sa ajunga in categories
        }
    }



    conditions <- admisc::recreate(substitute(conditions), colnames(data))
    if (identical(conditions, "")) {
        conditions <- setdiff(
            colnames(data),
            ifelse(
                grepl("\\+|\\*", outcome),
                outcome,
                admisc::notilde(outcome)
            )
        )
    }
    else {
        # if (is.character(conditions) & length(conditions) == 1) {
            conditions <- admisc::splitstr(conditions)
        # }
    }

    if (length(conditions) > 30) {
        admisc::stopError(
            "Impossible to run a QCA analysis with so many conditions.",
            ... = ...
        )
    }

    if (is.character(sort.by) & length(sort.by) == 1 & !identical(sort.by, "")) {
        sort.by <- admisc::splitstr(sort.by)
    }

    decreasing <- TRUE # initialize a default value
    if (is.element("decreasing", names(dots))) {
        decreasing <- dots$decreasing
    }

    if (is.character(decreasing) & length(decreasing) == 1) {
        decreasing <- admisc::splitstr(decreasing)
    }


    if (!identical(inf.test, "")) {
        inf.test <- admisc::splitstr(inf.test)
    }


    if (is.matrix(data)) {
        if (is.null(colnames(data))) {
            admisc::stopError(
                "The data should have column names.",
                ... = ...
            )
        }

        if (any(duplicated(rownames(data)))) {
            rownames(data) <- seq(nrow(data))
        }

        data <- as.data.frame(data)
        data[] <- lapply(data, function(x) {
            if (admisc::possibleNumeric(x)) {
                x <- admisc::asNumeric(x)
            }
            return(x)
        })

        initial.data <- data
    }

    verify.tt(
        data,
        outcome,
        conditions,
        complete,
        show.cases,
        ic1,
        ic0,
        inf.test,
        ... = ...
    )

    if (length(conditions) == 1) {
        if (grepl(":", conditions)) {
            # they are already verified as existing columns, in verify.tt()
            nms <- colnames(data)
            cs <- unlist(strsplit(conditions, split = ":"))
            conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
        }
    }


    if (!grepl("\\+|\\*", outcome)) {
        data <- data[, c(conditions, admisc::notilde(outcome))]

        if (neg.out | admisc::tilde1st(outcome)) {
            data[, admisc::notilde(outcome)] <- 1 - data[, admisc::notilde(outcome)]
        }

        outcome <- admisc::notilde(outcome)

    }


    nofconditions <- length(conditions)

    infodata <- admisc::getInfo(data[, conditions, drop = FALSE])

    data[, conditions] <- infodata$data # having treated the "don't care" cells
    hastime <- infodata$hastime
    fuzzy.cc <- infodata$fuzzy.cc
    noflevels <- infodata$noflevels
    dc.code <- infodata$dc.code
    rownames(data) <- rownames(initial.data)

    categories <- infodata$categories
    if (declared & !is.null(categories)) {
        categories[[outcome]] <- names(outcomelabels)
    }

    condata <- data[, conditions, drop = FALSE]

    if (any(fuzzy.cc)) {
        if (any(data[, conditions[fuzzy.cc]] == 0.5)) {
            warning(
                paste0(
                    "\n",
                    "Fuzzy causal conditions should not have values of 0.5 in the data.",
                    "\n\n"
                )
            )
        }
        condata[, fuzzy.cc] <- lapply(
            condata[, fuzzy.cc, drop = FALSE],
            function(x) as.numeric(x > 0.5)
        )
    }

    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    line.data <- as.vector(as.matrix(condata) %*% mbase) + 1
    condata <- condata[order(line.data), , drop = FALSE]
    uniq <- which(!duplicated(condata))
    tt <- condata[uniq, , drop = FALSE]


    rownstt <- sort(line.data)[uniq]
    rownames(tt) <- rownstt

    # return(list(X = as.matrix(data[, conditions]), Y = as.matrix(tt), FUZ = as.numeric(fuzzy.cc), VO = data[, outcome]))
    ipc <- .Call(
        "C_truthTable",
        as.matrix(data[, conditions]),
        admisc::compute(outcome, data),
        as.matrix(tt),
        as.numeric(fuzzy.cc),
        PACKAGE = "QCA"
    )

    colnames(ipc) <- rownstt

    minmat <- ipc[seq(4, nrow(ipc)), , drop = FALSE]
    ipc <- ipc[1:3, , drop = FALSE]

    rownames(minmat) <- rownames(data)
    rownames(ipc) <- c("n", "incl", "PRI")

    # obremove means to remove observed configurations
    obremove <- ipc[1, ] < n.cut # | ipc[3, ] < pri.cut

    if (sum(!obremove) == 0) {
        admisc::stopError(
            "There are no configurations, using these cutoff values.",
            ... = ...
        )
    }

    tt$OUT <- "?"
    tt$OUT[!obremove] <- 1 * (
        admisc::agteb(ipc[2, !obremove], ic1) &
        admisc::agteb(ipc[3, !obremove], pri.cut)
    )

    tt$OUT[ipc[2, !obremove] < ic1 & admisc::agteb(ipc[2, !obremove], ic0)] <- "C"

    tt <- cbind(tt, t(ipc))

    # detect cases that have values equal to exactly 0.5
    zero.five <- apply(
        data[, conditions, drop = FALSE],
        1,
        function(x) any(admisc::aeqb(x, 0.5))
    )

    cases <- sapply(rownstt, function(x) {
        # and prevent allocating those cases to truth table rows
        paste(rownames(data)[line.data == x & !zero.five], collapse = ",")
    })

    DCC <- apply(minmat, 2, function(x) {
        paste(rownames(data)[x > 0.5 & data[, outcome] < 0.5], collapse = ",")
    })

    casesexcl <- cases[obremove]

    removed <- tt[obremove, , drop = FALSE]
    removed$OUT <- 1 * (
        admisc::agteb(ipc[2, obremove], ic1) &
        admisc::agteb(ipc[3, obremove], pri.cut)

    )

    removed$OUT[
        ipc[2, obremove] < ic1 & admisc::agteb(ipc[2, obremove], ic0)
    ] <- "C"

    frcallist <- NULL

    if (!is.null(exclude)) {

        if (admisc::possibleNumeric(exclude)) {
            exclude <- admisc::asNumeric(exclude)
        }

        if (is.character(exclude)) {
            frargs <- setdiff(names(formals(findRows)), "...")
            callist <- as.list(metacall)
            frcallist <- list(expression = exclude)
            callist$exclude <- NULL
            common <- intersect(names(dots), frargs)
            if (length(common) > 0) {
                for (i in seq(length(common))) {
                    frcallist[[common[i]]] <- dots[[common[i]]]
                    callist[[common[i]]] <- NULL
                }
            }
            # DNF expression
            # transform it into row numbers
            # exclude <- findRows(exclude, conditions, noflevels) # remainders only, by default
            exclude <- NULL
        }

        if (length(exclude) > 0) {
            # exclude <- setdiff(exclude, rownames(tt)[is.element(tt$tt$OUT, 0:1)])
            exclude <- exclude[exclude <= prod(noflevels)]
        }
    }

    if (length(exclude) == 0) {
        exclude <- NULL
    }


    if (length(conditions) > 7) {
        rownstt <- rownstt[!obremove]
        cases <- cases[!obremove]
        DCC <- DCC[!obremove]
        tt <- tt[!obremove, , drop = FALSE]

        if (!is.null(exclude)) {
            excl.matrix <- as.data.frame(getRow(exclude, noflevels))
            rownames(excl.matrix) <- exclude
            colnames(excl.matrix) <- conditions
            excl.matrix$OUT <- 0
            excl.matrix$n <- 0
            excl.matrix$incl <- "-"
            excl.matrix$PRI <- "-"

            tt <- rbind(tt, excl.matrix)
        }
    }
    else {
        ttc <- as.data.frame(matrix(nrow = prod(noflevels), ncol = ncol(tt)))
        colnames(ttc) <- colnames(tt)
        ttc[, seq(length(conditions))] <- createMatrix(noflevels)
        ttc$OUT   <- "?"
        ttc$n     <-  0
        ttc$incl  <- "-"

        # sometimes a causal condition may be named PRI (see Porter data)
        whichpri <- which(colnames(ttc) == "PRI")

        ttc[, whichpri[length(whichpri)]] <- "-"
        ttc[rownames(tt), ] <- tt

        if (!is.null(exclude)) {
            ttc$OUT[exclude] <- "0"
        }

        tt <- ttc
    }


    if (!identical(sort.by, "")) {

        if (is.logical(sort.by)) { # & !is.null(names(sort.by)) # if logical, it should _always_ have names
            decreasing <- as.vector(sort.by)
            sort.by <- names(sort.by)
        }
        else {
            # just to make sure we _do_ have a "decreasing" object
            if (missing(decreasing)) {
                decreasing <- rep(TRUE, length(sort.by))
            }
            else {
                if (is.logical(decreasing)) {
                    if (length(decreasing) == 1) {
                        decreasing <- rep(decreasing, length(sort.by))
                    }
                    else if (length(decreasing) < length(sort.by)) {
                        decreasing <- c(
                            decreasing,
                            rep(
                                TRUE,
                                length(sort.by) - length(decreasing)
                            )
                        )
                    }
                }
                else {
                    decreasing <- rep(TRUE, length(sort.by))
                }
            }
        }

        sort.by[sort.by == "out"] <- "OUT"

        decreasing <- decreasing[is.element(sort.by, names(tt))]
        sort.by <- sort.by[is.element(sort.by, names(tt))]


        rowsorder <- seq_len(nrow(tt))
        for (i in rev(seq(length(sort.by)))) {
            rowsorder <- rowsorder[
                order(
                    tt[rowsorder, sort.by[i]],
                    decreasing = decreasing[i]
                )
            ]
        }

        sortvector <- rep(1, nrow(tt))
        sortvector[tt[rowsorder, "OUT"] == "?"] <- 2
        rowsorder <- rowsorder[order(sortvector)]

    }

    if (any(hastime) && length(dc.code) > 1) {
        admisc::stopError("Multiple \"don't care\" codes found.")
    }

    # TO SOLVE: perhaps using categorical QCA, dc.code (or -1) does not perform well with CCubes
    for (i in seq(length(conditions))) {
        if (hastime[i]) {
            tt[, i][tt[, i] == max(tt[, i])] <- dc.code
            data[, i][data[, i] == max(data[, i])] <- -1 # was dc.code
            noflevels[i] <- noflevels[i] - 1
        }
    }

    statistical.testing <- FALSE

    if (inf.test[1] == "binom") {
        statistical.testing <- TRUE
        if (length(inf.test) > 1) {
            alpha <- as.numeric(inf.test[2])
            # already checked if a number between 0 and 1
        }
        else {
            alpha <- 0.05
        }

        observed <- which(tt$OUT != "?")
        success <- round(tt[observed, "n"] * as.numeric(tt[observed, "incl"]))

        tt$pval1 <- "-"
        if (length(incl.cut) > 1) {
            tt$pval0 <- "-"
        }
        tt[observed, "OUT"] <- 0

        for (i in seq(length(observed))) {

            pval1 <- tt[observed[i], "pval1"] <- binom.test(
                success[i],
                tt[observed[i], "n"],
                p = ic1,
                alternative = "greater"
            )$p.value

            if (length(incl.cut) > 1) {
                pval0 <- tt[observed[i], "pval0"] <- binom.test(
                    success[i],
                    tt[observed[i], "n"],
                    p = ic0,
                    alternative = "greater"
                )$p.value
            }

            if (pval1 < alpha) {
                tt[observed[i], "OUT"] <- 1
            }
            else if (length(incl.cut) > 1) {
                if (pval0 < alpha) {
                    tt[observed[i], "OUT"] <- "C"
                }
            }
        }
    }


    # deal with the show.cases in the print function
    # if (show.cases) {

        # apparently this is necessary, otherwise
        # the new column cases will be a factor
        tt$cases <- ""

        if (length(conditions) < 8) {
            tt$cases[rownstt] <- cases
        }
        else {
            if (nrow(tt) > length(cases)) {
                cases <- c(cases, rep("", nrow(tt) - length(cases)))
            }
            tt$cases <- cases
        }

    # }

    numerics <- unlist(lapply(initial.data, admisc::possibleNumeric))
    colnames(initial.data)[!numerics] <- initialcols[!numerics]

    multivalue <- multivalue | any(noflevels > 2)
    x <- list(
        tt = tt,
        indexes = rownstt,
        noflevels = as.vector(noflevels),
        initial.data = initial.data,
        recoded.data = data,
        cases = cases,
        DCC = DCC,
        minmat = minmat,
        categories = categories,
        multivalue = multivalue,
        options = list(
            outcome = outcome.copy,
            conditions = conditions,
            neg.out = neg.out,
            n.cut = n.cut,
            incl.cut = incl.cut,
            pri.cut = pri.cut,
            exclude = exclude,
            complete = complete,
            show.cases = show.cases,
            dcc = dcc,
            use.letters = use.letters,
            use.labels = use.labels,
            inf.test = statistical.testing
        )
    )

    if (any(obremove)) {
        removed$cases <- ""
        removed$cases <- casesexcl
        x$removed <- structure(
            list(
                tt = removed,
                categories = categories,
                multivalue = multivalue,
                options = list(
                    complete = FALSE,
                    show.cases = show.cases,
                    dcc = dcc,
                    removed = TRUE,
                    use.labels = use.labels
                )
            ),
            class = "QCA_tt"
        )
    }

    # also verify if not already letters
    if (use.letters & any(nchar(conditions) > 1)) {
        colnames(x$tt)[seq(nofconditions)] <- LETTERS[seq(nofconditions)]
    }

    if (!identical(sort.by, "")) {
        x$rowsorder <- rowsorder
    }

    x$fs <- unname(fuzzy.cc)
    x$call <- metacall
    x <- structure(x, class = "QCA_tt")

    if (!is.null(frcallist)) {
        tempcall <- callist
        for (i in seq(2, length(tempcall))) {
            tc <- tryCatch(eval.parent(tempcall[[i]]), error = function(e) e)
            if (is.list(tc) && identical(names(tc), c("message", "call"))) {
                tc <- as.character(tempcall[[i]])
            }

            tempcall[[i]] <- tc
        }
        x$call <- as.call(tempcall)
        frcallist$obj <- x
        tempcall$exclude <- do.call("findRows", frcallist)
        x <- do.call("truthTable", tempcall[-1])
        callist$exclude <- tempcall$exclude
        x$call <- as.call(callist)
    }

    return(x)
}
