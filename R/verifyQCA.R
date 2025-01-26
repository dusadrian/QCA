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

`verify.data` <-
function(data, outcome = "", conditions = "") {
    if (!is.data.frame(data)) {
        admisc::stopError(
            "The input data should be a data frame."
        )
    }
    if (is.null(colnames(data))) {
        admisc::stopError(
            "Please specify the column names for your data."
        )
    }
    if (identical(outcome, "")) {
        admisc::stopError(
            "The outcome set is not specified."
        )
    }
    testoutcome <- admisc::tryCatchWEM(
        trout <- admisc::translate(outcome, colnames(data))
    )
    if (is.element("error", names(testoutcome))) {
        admisc::stopError(
            "Incorrect outcome specification."
        )
    }
    testrout <- apply(trout, 2, function(x) {
        all(x != "-1")
    })
    if (sum(testrout) == 1) {
        outcome <- names(testrout)[testrout]
    }
    if (!identical(conditions, "")) {
        if (any(grepl(":", conditions))) {
            if (length(conditions) > 1) {
                admisc::stopError(
                    "Only one sequence of conditions allowed."
                )
            }
            conditions <- unlist(strsplit(conditions, split = ":"))
            nms <- colnames(data)
            cs <- unlist(strsplit(conditions, split = ":"))
            conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
            if (is.element(outcome, conditions)) {
                admisc::stopError(
                    "Outcome found in the sequence of conditions."
                )
            }
        }
        if (is.element(outcome, conditions)) {
            admisc::stopError(
                paste0(
                    "\"",
                    outcome,
                    "\" cannot be both outcome _and_ condition."
                )
            )
        }
        if (!all(is.element(conditions, names(data)))) {
            admisc::stopError(
                "Conditions not found in the data."
            )
        }
        if (any(duplicated(conditions))) {
            admisc::stopError(
                "Duplicated conditions."
            )
        }
    }
    if (any(is.na(data))) {
        checked <- sapply(data, function(x) any(is.na(x)))
        admisc::stopError(
            paste0(
                "Missing values in the data are not allowed. ",
                "Please check columns:\n",
                paste(
                    names(checked)[checked],
                    collapse = ", "
                )
            )
        )
    }
}
`verify.qca` <-
function(data, ...) {
    dots <- list(...)
    simplify <- isTRUE(dots$simplify)
    if (is.data.frame(data)) {
        if (is.null(colnames(data))) {
            admisc::stopError(
                "The dataset doesn't have any columns names."
            )
        }
        checkNumUncal <- lapply(data, function(x) {
            is_x_factor <- is.factor(x)
            is_x_declared <- inherits(x, "declared")
            x <- setdiff(x, c("-", "dc", "?"))
            is_possible_numeric <- admisc::possibleNumeric(x)
            uncal <- mvuncal <- FALSE
            if (is_possible_numeric) {
                y <- na.omit(admisc::asNumeric(x))
                if (admisc::wholeNumeric(y)) {
                    mvuncal <- length(seq(0, max(y))) > 20 ||
                    (
                        !simplify && 
                        any(y > 1) && 
                        !setequal(
                            sort(unique(y)),
                            seq(0, max(y))
                        )
                    )
                }
                else {
                    uncal <- any(y > 1)
                }
            }
            return(c(is_possible_numeric, uncal, mvuncal, is_x_factor, is_x_declared))
        })
        checknumeric <- sapply(checkNumUncal, "[[", 1)
        checkuncal <- sapply(checkNumUncal, "[[", 2)
        checkmvuncal <- sapply(checkNumUncal, "[[", 3)
        checkfactor <- sapply(checkNumUncal, "[[", 4)
        checkdeclared <- sapply(checkNumUncal, "[[", 5)
        if (!all(checknumeric | checkfactor | checkdeclared)) {
            notnumeric <- colnames(data)[!checknumeric]
            errmessage <- paste(
                "The causal condition",
                ifelse(length(notnumeric) == 1, " ", "s "),
                paste(notnumeric, collapse=", "),
                ifelse(length(notnumeric) == 1, " is ", " are "),
                "not numeric or factor.",
                sep = ""
            )
            admisc::stopError(errmessage)
        }
        if (any(checkuncal)) {
            uncalibrated <- colnames(data)[checkuncal]
            errmessage <- paste0(
                "Uncalibrated data.\n       ",
                "Fuzzy sets should have values bound to the interval [0 , 1] ",
                "and all other sets should be crisp.\n       ",
                "Please check the following condition",
                ifelse(length(uncalibrated) == 1, "", "s"),
                ":\n",
                paste(uncalibrated, collapse = ", ")
            )
            admisc::stopError(errmessage)
        }
        if (any(checkmvuncal)) {
            admisc::stopError(paste(
                "Possibly uncalibrated multivalue conditions. Please check:\n",
                paste(colnames(data)[checkmvuncal], collapse = ", ")
            ))
        }
    }
    else if (is.vector(drop(data))) {
        if (!admisc::possibleNumeric(data)) {
            admisc::stopError(
                "Non numeric input."
            )
        }
    }
}
`verify.tt` <- function(
    data, outcome = "", conditions = "", complete = FALSE,
    show.cases = FALSE, ic1 = 1, ic0 = 1, inf.test, ...
) {
    if (!inherits(data, "data.frame")) {
        cls <- ifelse(methods::is(data, "QCA_sS"), "QCA_sS",
                ifelse(methods::is(data, "QCA_tt"), "QCA_tt",
                ifelse(methods::is(data, "QCA_pof"), "QCA_tt",
                paste(class(data), collapse = ", "))))
        errmessage <- paste0(
            "You have to provide a data frame, ",
            "the current \"data\" argument contains an object\n",
            "       of class ",
            class(data),
            ifelse(cls == "QCA_sS", ", created by superSubset()", ""),
            ifelse(cls == "QCA_tt", ", created by truthTable()", ""),
            ifelse(cls == "QCA_pof", ", created by pof()", ""),
            "."
        )
        admisc::stopError(errmessage)
    }
    if (methods::is(data, "QCA_tt")) {
        data <- data$initial.data
    }
    if (identical(outcome, "")) {
        admisc::stopError("Incorrect outcome specification.")
    }
    testoutcome <- admisc::tryCatchWEM(admisc::translate(outcome, colnames(data)))
    if (is.element("error", names(testoutcome))) {
        admisc::stopError("The outcome is not correct.")
    }
    if (!identical(conditions, "")) {
        if (length(conditions) == 1 & is.character(conditions)) {
            conditions <- admisc::splitstr(conditions)
            if (any(grepl(":", conditions)) & length(conditions) > 1) {
                admisc::stopError(
                    "Only one sequence of conditions allowed."
                )
            }
            conditions <- unlist(strsplit(conditions, split = ":"))
        }
        if (is.element(outcome, conditions)) {
            admisc::stopError(
                paste0(
                    "Variable \"",
                    outcome,
                    "\" cannot be both outcome _and_ condition!"
                )
            )
        }
        if (!all(is.element(conditions, names(data)))) {
            admisc::stopError(
                "Conditions not found in the data."
            )
        }
        if (any(duplicated(conditions))) {
            admisc::stopError(
                "Duplicated conditions."
            )
        }
    }
    else {
        conditions <- colnames(data)
        conditions <- setdiff(conditions, outcome)
    }
    if (any(is.na(data))) {
        checked <- sapply(data, function(x) any(is.na(x)))
        admisc::stopError(
            paste0(
                "Missing values in the data are not allowed. ",
                "Please check columns:\n",
                paste(
                    names(checked)[checked],
                    collapse = ", "
                )
            )
        )
    }
    if (any(c(ic1, ic0) < 0) | any(c(ic1, ic0) > 1)) {
        admisc::stopError(
            "The inclusion cut-off(s) should be bound to the interval [0, 1]."
        )
    }
    testoutcome <- admisc::tryCatchWEM(
        trout <- admisc::translate(outcome, colnames(data))
    )
    if (is.element("error", names(testoutcome))) {
        admisc::stopError("Incorrect outcome specification.")
    }
    testrout <- apply(trout, 2, function(x) {
        all(x != "-1")
    })
    data <- data[, unique(c(conditions, names(testrout)[testrout]))]
    data[] <- lapply(data, function(x) {
        if (!is.factor(x) & !inherits(x, "declared")) {
            x <- as.character(x)
            x[x %in% c("-", "dc", "?")] <- -1
            if (admisc::possibleNumeric(x)) {
                x <- admisc::asNumeric(x)
            }
        }
        return(x)
    })
    verify.qca(data, ... = ...)
    verify.inf.test(inf.test, data)
}
`verify.minimize` <-
function(data, outcome = "", conditions = "", explain = "",
         include = "", use.letters = FALSE) {
    if (all(explain == "")) {
        admisc::stopError(
            "You have not specified what to explain."
        )
    }
    if (any(explain == 0)) {
        admisc::stopError(
            "Negative output configurations cannot be explained."
        )
    }
    if (any(include == 0)) {
        admisc::stopError(
            paste(
                "Negative output configurations",
                "cannot be included in the minimization."
            )
        )
    }
    if (length(setdiff(explain, c(1, "C"))) > 0) {
        admisc::stopError(
            paste(
                "Only the positive output configurations",
                "and/or contradictions can be explained."
            )
        )
    }
    if (length(setdiff(include, c("?", "C", ""))) > 0) {
        admisc::stopError(
            paste(
                "Only the remainders and/or the contradictions",
                "can be included in the minimization."
            )
        )
    }
    if (is.element("C", explain) & is.element("C", include)) {
        admisc::stopError(
            "Contradictions are either explained or included, but not both."
        )
    }
    if (use.letters & ncol(data) > 27) {
        admisc::stopError(
            "Cannot use letters. There are more than 26 conditions."
        )
    }
    if (any(is.na(data))) {
        checked <- sapply(data, function(x) any(is.na(x)))
        admisc::stopError(
            paste0(
                "Missing values in the data are not allowed. ",
                "Please check columns:\n",
                paste(
                    names(checked)[checked],
                    collapse = ", "
                )
            )
        )
    }
}
`verify.dir.exp` <- function(
    data, outcome, conditions, noflevels, dir.exp = "", enter = NULL
) {
    if (is.null(enter)) enter <- "\n"
    if (is.null(dir.exp)) {
        return(dir.exp)
    }
    else {
        multivalue <- any(grepl(mvregexp, dir.exp))
        if (is.character(dir.exp)) {
            dir.exp <- gsub(paste(admisc::dashes(), collapse = "|"), "-", dir.exp)
        }
        if (identical(dir.exp, "")) {
            dir.exp <- paste(rep("-", length(conditions)), collapse = ",")
        }
        direxpsplit <- unlist(
            strsplit(
                gsub("[-|;|,|[:space:]]", "", dir.exp),
                split = ""
            )
        )
        oldway <- admisc::possibleNumeric(direxpsplit) | length(direxpsplit) == 0
        if (oldway) {
            if (length(dir.exp) == 1) {
                dir.exp <- admisc::splitstr(dir.exp)
            }
            expression <- NULL
            if (length(dir.exp) != length(conditions)) {
                admisc::stopError(
                    "Number of expectations does not match number of conditions."
                )
            }
            if (all(dir.exp == "-")) {
                return(matrix(0L, ncol = length(conditions)))
            }
            del <- strsplit(as.character(dir.exp), split = ";")
            if (is.null(names(dir.exp))) {
                names(del) <- conditions
            }
            else {
                if (length(names(dir.exp)) != length(conditions)) {
                    admisc::stopError(
                        "All directional expectations should have names, or none at all."
                    )
                }
                else if (length(setdiff(names(dir.exp), conditions)) > 0) {
                    admisc::stopError(
                        "Incorect names of the directional expectations."
                    )
                }
                names(del) <- names(dir.exp)
                del <- del[conditions]
            }
            for (i in seq(length(del))) {
                values <- del[[i]]
                if (any(values != "-")) {
                    values <- admisc::asNumeric(setdiff(values, "-"))
                    if (length(setdiff(values, seq(noflevels[i]) - 1)) > 0) {
                        errmessage <- paste0(
                            'Values specified in the directional expectations ',
                            'do not appear in the data, for condition "',
                            conditions[i],
                            '".'
                        )
                        admisc::stopError(errmessage)
                    }
                    else {
                        expression <- c(
                            expression,
                            paste0(conditions[i], "[", values, "]")
                        )
                    }
                }
            }
            multivalue <- TRUE
            dir.exp <- expression
        }
        else {
            if (length(dir.exp) == 1) {
                if (!grepl("[+]", dir.exp) &  grepl("[,]", dir.exp)) {
                    if (multivalue) {
                        values <- admisc::squareBrackets(dir.exp)
                        atvalues <- paste("@", seq(length(values)), sep = "")
                        for (i in seq(length(values))) {
                            dir.exp <- gsub(values[i], atvalues[i], dir.exp)
                        }
                        dir.exp <- gsub(",", "+", dir.exp)
                        for (i in seq(length(values))) {
                            dir.exp <- gsub(atvalues[i], values[i], dir.exp)
                        }
                    }
                    else {
                        dir.exp <- gsub(",", "+", dir.exp)
                    }
                }
            }
        }
        dir.exp <- paste(dir.exp, collapse = "+") 
        if (!multivalue) {
            if (any(noflevels > 2)) {
                admisc::stopError(
                    paste(
                        "For multivalue data, directional expectations",
                        "should be specified using square brackets."
                    )
                )
            }
        }
        if (!oldway) {
            dir.exp <- tryCatch(
                admisc::simplify(
                    expression = dir.exp,
                    snames = conditions,
                    noflevels = noflevels,
                    dir.exp = TRUE
                ),
                error = function(e) e,
                warning = function(w) w
            )
        }
        if (length(dir.exp) > 1) {
            admisc::stopError(
                "Ambiguous directional expectations."
            )
        }
        if (identical(dir.exp, "")) {
            admisc::stopError(
                "Directional expectations cancel each other out to an empty set."
            )
        }
        dir.exp <- admisc::translate(
            dir.exp,
            snames = conditions,
            noflevels = noflevels
        )
        return(matrix(as.integer(dir.exp) + 1L, ncol = ncol(dir.exp)))
    }
}
`verify.mqca` <-
function(allargs) {
    data <- allargs$input
    outcome <- admisc::splitstr(allargs$outcome)
    mvoutcome <- grepl(mvregexp, outcome) 
    if (any(mvoutcome)) {
        if (any(grepl("\\{", outcome))) {
            outcome.value <- admisc::curlyBrackets(outcome)
            outcome <- admisc::curlyBrackets(outcome, outside = TRUE)
        }
        else {
            outcome.value <- admisc::squareBrackets(outcome)
            outcome <- admisc::squareBrackets(outcome, outside = TRUE)
        }
        if (length(setdiff(outcome, names(data))) > 0) {
            outcome <- setdiff(outcome, names(data))
            errmessage <- paste0(
                'Outcome(s) not present in the data: "',
                paste(outcome, collapse = '", "'),
                '".'
            )
            admisc::stopError(errmessage)
        }
        for (i in seq(length(outcome))) {
            if (mvoutcome[i]) {
                mvnot <- setdiff(
                    admisc::splitstr(outcome.value[i]),
                    unique(data[, outcome[i]])
                )
                if (length(mvnot) > 0) {
                    admisc::stopError(
                        sprintf(
                            'Value(s) %s not found in the outcome "%s".',
                            paste(mvnot, collapse = ","),
                            outcome[i]
                        )
                    )
                }
            }
        }
    }
    else {
        if (length(setdiff(outcome, names(data))) > 0) {
            outcome <- setdiff(outcome, names(data))
            admisc::stopError(
                paste0(
                    'Outcome(s) not present in the data: "',
                    paste(outcome, collapse='", "'),
                    '".'
                )
            )
        }
        fuzzy.outcome <- apply(
            data[, outcome, drop=FALSE],
            2,
            function(x) any(x %% 1 > 0)
        )
        if (any(!fuzzy.outcome)) {
            outcome.copy <- outcome[!fuzzy.outcome]
            for (i in outcome.copy) {
                valents <- unique(data[, i])
                if (!all(valents %in% c(0, 1))) {
                    admisc::stopError(
                        paste0(
                            'Please specify the value of outcome variable "',
                            i,
                            '" to explain.'
                        )
                    )
                }
            }
        }
    }
    conditions <- allargs$conditions
    if (is.null(conditions)) {
        conditions <- names(data)
    }
    else {
        conditions <- admisc::splitstr(conditions)
    }
    if (length(setdiff(outcome, conditions)) > 0) {
        outcome <- setdiff(outcome, conditions)
        admisc::stopError(
            paste0(
                "Outcome(s) not present in the conditions' names: \"",
                paste(outcome, collapse = '", "'),
                '".'
            )
        )
    }
}
`verify.inf.test` <- function(inf.test, data) {
    if (all(inf.test != "")) {
        if (inf.test[1] != "binom") {
            admisc::stopError(
                "For the moment only \"binom\"ial testing for crisp data is allowed."
            )
        }
        else {
            fuzzy <- apply(data, 2, function(x) any(x %% 1 > 0))
            if (any(fuzzy)) {
                admisc::stopError(
                    "The binomial test only works with crisp data."
                )
            }
        }
        if (length(inf.test) > 1) {
            alpha <- as.numeric(inf.test[2])
            if (is.na(alpha) | alpha < 0 | alpha > 1) {
                admisc::stopError(
                    "The second value of inf.test should be a number between 0 and 1."
                )
            }
        }
    }
}
