`pof` <- function(
    setms = NULL, outcome = NULL, data = NULL, relation = "necessity",
    use.labels = FALSE, inf.test = "", incl.cut = c(0.75, 0.5), add = NULL, ...
) {

    # Objective: detect if an input (setms or outcome) was negated using 1 -
    # and present its name with a tilde in front in the output

    # this is an usually easy task via funargs, for instance if an outcome named "Y"
    # was negated, then funargs would show:
    # $outcome
    # [1] "1 - Y"

    # however, if the function pof() is called via do.call(), this is no longer the case
    # and funargs shows:
    # $outcome
    # [1] "c(0.47,0.7,0.03,0.03,0.03,0.27,0.15,0.44,0.47,0.36,"
    # [2] "0.6, 0.38)"

    ## Below an attempt to detect A -> B that is inverted into B <- A by the parser
    ## srcref is an attribute returned by the parser *only* from sourced files
    ## still to figure out how to deal with multiple expressions: c("A -> B", "A*B -> C")
    ## etc.
    # setms <- substitute(setms)
    # srcref <- as.character(attr(sys.call(), "srcref"))
    # if (length(srcref) > 0) {
    #     srcref <- unlist(strsplit(srcref, split = "pof\\("))[2]
    #     srcref <- unlist(strsplit(srcref, split = ","))[1]
    #     srcref <- gsub(paste0("\"|\'|[[:space:]|)|", "\u00a0", "]"), "", srcref)
    # }
    ## Then srcref should be passed to recreate() via the three dots ...

    outcome <- admisc::recreate(substitute(outcome), snames = names(data))

    if (is.null(setms) & !is.null(data)) {
        return(pofind(
            data = data,
            outcome = outcome,
            relation = relation,
            use.labels = use.labels,
            ... = ...
        ))
    }

    setms <- admisc::recreate(substitute(setms), snames = names(data))

    if (inherits(outcome, "declared")) {
        attributes(outcome) <- NULL
    }

    funargs <- lapply(
        lapply(match.call(), deparse)[-1],
        function(x) gsub(paste0("'|\"|[[:space:]|", "\u00a0", "]"), "", x)
    )

    dots <- list(...)

    if (isTRUE(dots$categorical)) { # backwards compatibility
        use.labels <- TRUE
        dots$categorical <- NULL
    }

    # because of the objective above, otherwise when doing
    # tilde1st(funargs$outcome) it would throw a warning
    funargs$outcome <- paste(funargs$outcome, collapse = "")

    if (is.null(setms)) {
        admisc::stopError(
            "The argument <setms> is missing.", ... = ...
        )
    }

    if (!(nec(relation) | suf(relation))) {
        admisc::stopError(
            "The relation should be either \"necessity\" or \"sufficiency\".",
            ... = ...
        )
    }

    # making sure the defaults are there
    ic1 <- 0.75
    ic0 <- 0.5

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

        if (is.element("incl.cut1", names(dots)) & identical(ic1, 0.75)) {
            ic1 <- dots$incl.cut1
        }

        if (is.element("incl.cut0", names(dots)) & identical(ic0, 0.5)) {
            ic0 <- dots$incl.cut0
        }
    ###
    ### ### backwards compatibility
    ###

    complete <- FALSE
    if (is.element("complete", names(dots))) {
        if (is.logical(dots$complete)) {
            complete <- dots$complete
        }
    }

    odata <- data
    infodata <- NULL
    categories <- list()

    if (is.element("categories", names(dots))) {
        categories <- dots$categories
        dots$categories <- NULL
    }

    if (!is.null(data)) {
        if (is.element("data.frame", class(data)) | is.matrix(data)) {
            # just in case it was a tibble, produced by package haven
            data <- as.data.frame(data)
        }

        infodata <- admisc::getInfo(data)
        categories <- infodata$categories

        data <- infodata$data

        if (is.element("minimize", names(dots))) {
            # the outcome is surely arranged in the last column
            if (is.element("use.letters", names(dots))) {
                if (dots$use.letters) {
                    colnames(data)[seq(1, ncol(data) - 1)] <- LETTERS[seq(1, ncol(data) - 1)]
                }
            }
        }
    }

    conditions <- outcomename <- ""
    condnegated <- outnegated <- FALSE


    ### temporary use of validateNames() until package admisc gets updated

    `extract` <- function(x, snames = "", data = NULL) {
        if (grepl("<=>|<->", x)) {
            admisc::stopError(
                "Incorrect expression: relation can be either necessity or sufficiency.",
                ... = ...
            )
        }

        multivalue <- grepl("\\{|\\}|\\[|\\]", x)
        relation <- ifelse(grepl("=|-", x), ifelse(grepl("=>|->", x), "suf", "nec"), NA)

        x <- gsub("<=|=>|<-|->", "@", gsub("[[:space:]]", "", x))
        x <- unlist(strsplit(x, split = "@"))
        xcopy <- x

        if (!multivalue) {
            x[1] <- mvSOP(x[1], snames = snames, data = data)
            if (!is.na(x[2])) {
                x[2] <- mvSOP(x[2], snames = snames, data = data)
            }
        }

        if (grepl("\\+|\\*", x[2])) {
            # if (grepl("\\+|\\*", x[1])) {
            #     admisc::stopError(
            #         "Incorrect output in the right hand side.", ... = ...
            #     )
            # }

            x <- rev(x)
            if (relation == "nec") {
                relation <- "suf"
            }
            else if (relation == "suf") {
                relation <- "nec"
            }
        }

        if (identical(snames, "") & !is.null(data)) {
            snames <- colnames(data)
        }

        if (identical(substring(x[1], 1, 2), "1-")) {
            x[1] <- negate(gsub("1-", "", x[1]), snames = snames)
        }

        if (identical(substring(x[2], 1, 2), "1-")) {
            x[2] <- negate(gsub("1-", "", x[2]), snames = snames)
        }

        outmtrx <- NA

        if (length(x) > 1) {
            outmtrx <- validateNames(x[2], snames = snames, data = data)
        }

        if (!is.na(outmtrx)) {
            if (!multivalue) {
                rownames(outmtrx) <- xcopy[2]
            }

            if (!is.null(data)) {
                # pof("SR~V + CLP~V + SCRP + SLRP => ~JSR", data = d.jobsecurity)
                # here, there are no * signs in the expression, and d.jobsecurity
                # contains both single letter conditions and multi-letter outcome
                # this decreases time when translating the expression, via admisc::validateNames()
                data <- data[, -which(is.element(colnames(data), colnames(outmtrx))), drop = FALSE]
            }
        }

        condmtrx <- validateNames(x[1], snames = snames, data = data)
        if (!multivalue & is.data.frame(condmtrx)) {
            rownames(condmtrx) <- admisc::trimstr(unlist(strsplit(xcopy[1], split = "\\+")))
            # x[1] <- xcopy[1]
        }

        return(
            list(
                condmtrx = condmtrx,
                outmtrx = outmtrx,
                expression = x[1],
                oexpr = xcopy[1],
                relation = relation,
                multivalue = multivalue
            )
        )
    }

    checkoutcome <- TRUE
    addexpression <- FALSE

    if (is.element("character", class(setms))) {
        # a character setms does not throw an evaluation error
        # it is not a formula, and not an R object

        if (missing(data)) {
            admisc::stopError(
                "The data argument is missing, with no default.", ... = ...
            )
        }

        if (length(setms) > 1) {
            admisc::stopError(
                "Only one expression allowed.", ... = ...
            )
        }

        toverify <- extract(setms, data = odata)

        if (!is.na(toverify$relation)) {
            relation <- toverify$relation
        }

        conditions <- colnames(toverify$condmtrx)
        if (is.na(toverify$outmtrx)) {

            if (missing(outcome)) {
                admisc::stopError(
                    "Expression without outcome.", ... = ...
                )
                # outcome <- "Y" # a generic name in case nothing else is found
            }

            temp <- subset(
                data,
                select = which(
                    is.element(
                        colnames(data),
                        conditions
                    )
                )
            )

            verify.qca(temp)

            # return(list(expression = toverify$expression, data = temp))
            setms <- admisc::compute(
                toverify$expression,
                data = temp,
                separate = TRUE
            )

            if (!toverify$multivalue & is.data.frame(setms)) {
                colnames(setms) <- admisc::trimstr(unlist(strsplit(toverify$oexpr, split = "\\+")))
            }

            # data <- data[, which(is.element(colnames(data), conditions)), drop = FALSE]
            funargs$setms <- toverify$expression

        }
        else {

            # here, the outcome was already verified by function extract()
            # that is part of the data
            outcomename <- colnames(toverify$outmtrx)

            temp <- subset(
                data,
                select = which(
                    is.element(
                        colnames(data),
                        c(conditions, outcomename)
                    )
                )
            )

            verify.qca(temp)

            # return(list(expression = toverify$expression, data = data[, -which(colnames(data) == outcomename)]))
            setms <- admisc::compute(toverify$expression, data = temp, separate = TRUE)

            if (!toverify$multivalue & is.data.frame(setms)) {
                colnames(setms) <- admisc::trimstr(unlist(strsplit(toverify$oexpr, split = "\\+")))
            }

            funargs$setms <- paste(
                paste(
                    unlist(toverify$expression),
                    collapse = "+"
                ),
                ifelse(
                    toverify$relation == "suf",
                    "->",
                    "<-"
                ),
                rownames(toverify$outmtrx)
            )

            # outcome <- as.numeric(is.element(data[, colnames(toverify$outcome)[1]], admisc::splitstr(toverify$outcome[1, 1])))
            outcome <- admisc::compute(
                rownames(toverify$outmtrx)[1], # [1] is redundant of course, but just in case
                data = temp
            )

            checkoutcome <- FALSE

        }


        if (is.vector(setms)) {
            setms <- data.frame(setms)
            colnames(setms) <- toverify$oexpr
        }

        rownames(setms) <- rownames(data)

        if (!is.element("minimize", names(dots)) & ncol(setms) > 1) {
            addexpression <- TRUE
        }
    }

    if (is.element("QCA_fuzzy", class(setms))) {
        # conditions <- attr(setms, "name")
        conditions <- "expression"
        setms <- data.frame(X = as.vector(setms))
        colnames(setms) <- conditions
    }


    if (checkoutcome) {
        if (missing(outcome)) {
            admisc::stopError(
                "Outcome is missing, with no default.", ... = ...
            )
        }

        if (is.element("character", class(outcome))) {

            # funargs$outcome and outcome are both character and contain the same thing
            # but funargs$outcome already has the [[:space:]] deleted

            if (grepl("\\+|\\*", outcome)) {
                outcomename <- outcome
            }
            else {
                if (admisc::tilde1st(gsub("1-", "", funargs$outcome))) {
                    outnegated <- !outnegated
                }

                oneminus <- identical(substr(funargs$outcome, 1, 2), "1-")

                if (oneminus) {
                    outnegated <- !outnegated
                    outcome <- gsub("1-", "", funargs$outcome)
                }

                outcome <- admisc::notilde(outcome)

                if (grepl("\\{", outcome)) {
                    outcomename <- admisc::curlyBrackets(outcome, outside = TRUE)
                }
                else {
                    outcomename <- admisc::squareBrackets(outcome, outside = TRUE)
                }
            }

            if (is.null(data)) {
                admisc::stopError(
                    "The data argument is missing, with no default.", ... = ...
                )
            }

            if (!is.element(outcomename, colnames(data))) {
                admisc::stopError(
                    "Outcome not found in the data.", ... = ...
                )
            }

            verify.qca(
                data[, which(colnames(data) == outcomename), drop = FALSE]
            )

            outcome <- admisc::compute(outcome, data = data)

            if (outnegated) {
                outcome <- 1 - outcome
            }

        }
        else if (is.vector(outcome)) {

            # Something like do.call("pof", list(1 - c(1,0,0,1,0), c(0,1,1,1,0)))
            # would NOT detect the negation of setms, which is OK
            # suppose we would have
            # ll <- list(1 - c(1,0,0,1,0), c(0,1,1,1,0))
            # do.call("pof", ll)
            # here, too, setms is negated and it would not be detected, hence OK
            # who uses a list instead of objects of conditions from a dataset, gets this

            if (admisc::tilde1st(gsub("1-", "", funargs$outcome))) {
                outnegated <- !outnegated
            }

            if (identical(substr(funargs$outcome, 1, 2), "1-")) {
                outnegated <- !outnegated
            }

            outcomename <- admisc::notilde(gsub("1-", "", funargs$outcome))

            if (identical(substr(outcomename, 1, 2), "c(")) {
                outcomename <- "Y"
            }
        }
    }


    if (is.vector(outcome)) {
        if (!is.numeric(outcome) & admisc::possibleNumeric(outcome)) {
            outcome <- admisc::asNumeric(outcome)
        }
        # if (identical(outcomename, "")) {
        # outcomename <- "Y"

        verify.qca(outcome)

        # if (!inherits(tc <- tryCatch(getName(funargs$outcome), error = function(e) e), "error")) {
        #     outcomename <- tc
        # }
        # }
    }
    else {
        admisc::stopError(
            paste(
                "The outcome should be either a column name in a dataset",
                "       or a vector of set membership values.",
                sep = "\n"
            ),
            ... = ...
        )
    }

    if (identical(substr(funargs$setms, 1, 2), "1-")) {
        condnegated <- !condnegated
    }

    if (is.vector(setms)) {
        setms <- data.frame(setms)

        conditions <- admisc::notilde(gsub("1-", "", funargs$setms))

        if (grepl("[$]", conditions)) {
            conditions <- tail(unlist(strsplit(conditions, split = "[$]")), 1)
        }
        else if (identical(substr(conditions, 1, 2), "c(")) {
            conditions <- "X"
        }

        colnames(setms) <- conditions
    }

    if (is.element("data.frame", class(setms))) {

        for (i in seq(ncol(setms))) {
            if (!is.numeric(setms[, i]) & admisc::possibleNumeric(setms[, i])) {
                setms[, i] <- admisc::asNumeric(setms[, i])
            }
        }

        verify.qca(setms)

        colnames(setms) <- gsub("[[:space:]]", "", colnames(setms))

        if (identical(conditions, "")) {
            conditions <- all.vars(parse(text = paste(colnames(setms), collapse = "+")))
        }


        if (condnegated) {
            conditions <- all.vars(parse(text = paste(conditions, collapse = "+")))

            # this is necessary when if the setms is a pims or coms component and it should be negated
            # in which case we need the set names, something like:
            # pof(1 - cCVF$pims, PROTEST, relation = "sufficiency")

            # assuming that setms is a pims or coms component, they should have a "conditions" attribute
            if (any(grepl("\\$coms|\\$pims", funargs$setms))) {
                toverify <- unlist(strsplit(admisc::notilde(gsub("1-", "", funargs$setms)), split = "\\$"))[1]
                if (grepl("pims", funargs$setms)) { # minimize()
                    tt <- eval.parent(
                        parse(text = sprintf("%s$tt", toverify)),
                        n = 1
                    )
                    if (tt$options$use.letters) {
                        conditions <- LETTERS[seq(length(conditions))]
                    }
                    else {
                        conditions <- tt$options$conditions
                    }
                }
                else {
                    conditions <- eval.parent(
                        parse(text = sprintf("%s$options$conditions", toverify)),
                        n = 1
                    )
                }
            }

            # in the meantime, conditions might have changed
            if (identical(conditions, "")) {
                colnames(setms) <- paste("~", colnames(setms), sep = "")
            }
            else {
                # return(list(input = colnames(setms), snames = conditions))
                colnames(setms) <- gsub("[[:space:]]", "", admisc::negate(colnames(setms), snames = conditions))
            }
        }
    }
    else {
        admisc::stopError(
            "The argument <setms> is not standard.", ... = ...
        )
    }

    if (any(na.omit(cbind(setms, outcome) > 1))) {
        admisc::stopError(
            "Set membership scores should be numbers between 0 and 1.",
            ... = ...
        )
    }

    notmiss <- apply(cbind(setms, outcome), 1, function(x) !any(is.na(x)))

    outcome <- outcome[notmiss]
    setms <- setms[notmiss, , drop = FALSE]

    if (neg.out) {
        outcome <- admisc::getLevels(outcome) - outcome - 1
    }

    result.list <- list()
    # return(list(as.matrix(cbind(setms, fuzzyor(setms))), outcome, nec(relation)))
    incl.cov <- .Call("C_pof", as.matrix(cbind(setms, fuzzyor(setms))), outcome, nec(relation), PACKAGE = "QCA")
    incl.cov[incl.cov < 0.00001] <- 0 # to take care of double precision numbers
    incl.cov <- as.data.frame(incl.cov)


    if (nec(relation)) {
        colnames(incl.cov) <- c("inclN", "RoN", "covN", "covU")
    }
    else {
        colnames(incl.cov) <- c("inclS", "PRI", "covS", "covU")
    }

    if (is.character(inf.test) & length(inf.test) == 1) {
        inf.test <- admisc::splitstr(inf.test)
    }

    if (!identical(inf.test, "")) {
        if (missing(data)) {
            data <- cbind(setms, outcome)
            colnames(data) <- c(conditions, outcomename)
        }
        verify.inf.test(inf.test, data)
    }

    if (identical(inf.test[1], "binom")) {

        statistical.testing <- TRUE

        if (length(inf.test) > 1) {
            alpha <- as.numeric(inf.test[2]) # already checked if a number between 0 and 1
        }
        else {
            alpha <- 0.05
        }

        if (nec(relation)) {
            nofcases <- rep(sum(outcome), ncol(setms) + 1)
        }
        else {
            nofcases <- c(colSums(setms), sum(fuzzyor(setms)))
        }

        success <- as.vector(round(nofcases * incl.cov[, which(grepl("incl", colnames(incl.cov)))[1]]))

        incl.cov$pval0 <- incl.cov$pval1 <- 0

        for (i in seq(length(success))) {
            incl.cov[i, "pval1"] <- binom.test(success[i], nofcases[i], p = ic1, alternative = "greater")$p.value
            incl.cov[i, "pval0"] <- binom.test(success[i], nofcases[i], p = ic0, alternative = "greater")$p.value
        }

    }

    result.list$incl.cov <- incl.cov

    if (nec(relation)) {
        result.list$incl.cov <- result.list$incl.cov[, -4]
    }
    else {
        result.list$incl.cov[nrow(incl.cov), 4] <- NA
    }

    colnms <- colnames(setms)

    if (addexpression) {
        colnms <- c(colnms, "expression")
    }
    else {
        result.list$incl.cov <- result.list$incl.cov[-nrow(incl.cov), , drop = FALSE]
        if (nrow(result.list$incl.cov) == 1 & suf(relation)) {
            result.list$incl.cov[1, 4] <- NA
        }
    }

    rownames(result.list$incl.cov) <- colnms


    if (is.element("show.cases", names(dots))) {
        if (dots$show.cases) {
            result.list$incl.cov <- cbind(result.list$incl.cov, cases = dots$cases, stringsAsFactors = FALSE)
        }
    }


    if (is.element("minimize", names(dots))) {
        result.list$pims <- as.data.frame(setms)
        result.list$sol.incl.cov <- incl.cov[nrow(incl.cov), 1:3]
    }


    if (is.element("solution.list", names(dots))) {

        solution.list <- dots$solution.list
        length.solution <- length(solution.list)
        individual <- vector("list", length = length.solution)

        for (i in seq(length.solution)) {
            individual[[i]] <- list()
            temp <- setms[, solution.list[[i]], drop = FALSE]

            incl.cov <- .Call("C_pof", as.matrix(cbind(temp, fuzzyor(temp))), outcome, nec(relation), PACKAGE = "QCA")
            incl.cov[incl.cov < 0.0001] <- 0
            incl.cov <- as.data.frame(incl.cov)

            rownames(incl.cov) <- c(colnames(temp), "expression")

            if (nec(relation)) {
                colnames(incl.cov) <- c("inclN", "RoN", "covN", "covU")
                incl.cov <- incl.cov[, -4]
            }
            else {
                colnames(incl.cov) <- c("inclS", "PRI", "covS", "covU")
                incl.cov[nrow(incl.cov), 4] <- NA
            }

            if (nrow(incl.cov) == 2 & suf(relation)) {
                incl.cov[1, 4] <- NA
            }

            individual[[i]]$incl.cov <- incl.cov[-nrow(incl.cov), ]
            individual[[i]]$sol.incl.cov <- incl.cov[nrow(incl.cov), 1:3]
            individual[[i]]$pims <- as.data.frame(temp)

        }

        # return(list(overall = result.list, individual = individual))
        return(structure(list(
            overall = result.list,
            individual = individual,
            essential = dots$essential,
            pims = as.data.frame(setms),
            relation = relation,
            categories = categories,
            options = c(
                list(
                    setms = setms,
                    outcome = outcome,
                    data = data,
                    relation = relation,
                    inf.test = inf.test,
                    incl.cut = incl.cut,
                    add = add,
                    use.labels = use.labels
                ),
                dots)
            ), class = "QCA_pof"))
    }

    result.list$categories <- categories

    if (!is.null(add)) {

        if (!(is.list(add) | is.function(add))) {
            admisc::stopError(
                "The argument <add> should be a function or a list of functions.",
                ... = ...
            )
        }

        if (is.list(add)) {
            if (!all(unlist(lapply(add, is.function)))) {
                admisc::stopError(
                    "Components from the list argument <add> should be functions.",
                    ... = ...
                )
            }

            toadd <- matrix(nrow = nrow(incl.cov), ncol = length(add))
            if (is.null(names(add))) {
                names(add) <- paste0("X", seq(length(add)))
            }

            if (any(duplicated(substr(names(add), 1, 5)))) {
                names(add) <- paste0("X", seq(length(add)))
            }

            colnames(toadd) <- substr(names(add), 1, 5)

            for (i in seq(length(add))) {
                coltoadd <- apply(
                    cbind(setms, fuzzyor(setms)),
                    2,
                    add[[i]],
                    outcome
                )
                if (ncol(setms) == 1) {
                    coltoadd <- coltoadd[1]
                }
                toadd[, i] <- coltoadd
            }
        }
        else {
            toadd <- matrix(nrow = nrow(incl.cov), ncol = 1)
            coltoadd <- apply(cbind(setms, fuzzyor(setms)), 2, add, outcome)
            if (ncol(setms) == 1) {
                coltoadd <- coltoadd[1]
            }
            toadd[, 1] <- coltoadd
            if (any(grepl("function", funargs$add))) {
                funargs$add <- "X"
            }
            colnames(toadd) <- substr(funargs$add, 1, 5)
        }

        result.list$incl.cov <- cbind(result.list$incl.cov, toadd)

    }

    result.list$options <- c(
        list(
            setms = setms,
            outcome = outcome,
            data = data,
            relation = relation,
            inf.test = inf.test,
            incl.cut = incl.cut,
            add = add,
            use.labels = use.labels
        ),
        dots)

    return(structure(result.list, class = "QCA_pof"))

}
