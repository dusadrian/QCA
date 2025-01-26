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

`calibrate` <- function (
    x, type = "fuzzy", method = "direct", thresholds = NA,
    logistic = TRUE, idm = 0.95, ecdf = FALSE, below = 1, above = 1, ...
) {
    dots <- list(...)
    funargs <- lapply(
        lapply(match.call(), deparse)[-1],
        function(x) gsub(paste0("'|\"|[[:space:]|", "\u00a0", "]"), "", x)
    )
    if (is.element("q", names(dots))) {
        above <- dots$q
    }
    if (is.element("p", names(dots))) {
        below <- dots$p
    }
    if (admisc::possibleNumeric(x)) {
        x <- admisc::asNumeric(x) 
    }
    else {
        if (grepl("[$]", funargs$x) & is.null(x)) {
            x <- unlist(strsplit(funargs$x, split = "\\$"))
            admisc::stopError(
                sprintf(
                    "There is no column \"%s\" in the dataframe %s.",
                    x[2],
                    x[1]
                ),
                ... = ...
            )
        }
        admisc::stopError("The input is not numeric.", ... = ...)
    }
    if (!is.element(type, c("crisp", "fuzzy"))) {
        admisc::stopError("Incorrect calibration type.", ... = ...)
    }
    if (!is.element(method, c("direct", "indirect", "TFR"))) {
        admisc::stopError("Incorrect calibration method.", ... = ...)
    }
    if (method != "TFR") {
        if(all(is.na(thresholds))) {
            admisc::stopError("Threshold value(s) not specified.", ... = ...)
        }
        if (is.character(thresholds) & length(thresholds) == 1) {
            thresholds <- admisc::splitstr(thresholds)
        }
        if (admisc::possibleNumeric(thresholds)) {
            nmsths <- NULL
            if (!is.null(names(thresholds))) {
                nmsths <- names(thresholds)
            }
            thresholds <- admisc::asNumeric(thresholds)
            names(thresholds) <- nmsths
        }
        else {
            admisc::stopError("Thresholds must be numeric.", ... = ...)
        }
    }
    if (type == "crisp") {
        if (any(thresholds < min(x) | thresholds > max(x))) {
            admisc::stopError(
                "Threshold value(s) outside the range of x.", ... = ...
            )
        }
        if (!is.null(names(thresholds))) {
            admisc::stopError(
                "Named thresholds require fuzzy type calibration.", ... = ...
            )
        }
        thresholds <- sort(thresholds)
        return(findInterval(x, thresholds))
    }
    else if (type == "fuzzy") {
        check.equal <- function(x, y) {
            check.vector <- as.logical(unlist(lapply(x, all.equal, y)))
            check.vector[is.na(check.vector)] <- FALSE
            return(check.vector)
        }
        lth <- length(thresholds)
        nth <- names(thresholds)
        if (method == "direct") {
            if (lth != 3 & lth != 6) {
                admisc::stopError(
                    "For fuzzy direct calibration, there should be either 3 or 6 thresholds.",
                    ... = ...
                )
            }
            if (idm <= 0.5 | idm >= 1) {
                admisc::stopError(
                    "The inclusion degree of membership has to be bigger than 0.5 and less than 1.",
                    ... = ...
                )
            }
            if (lth == 3) {
                if (!is.null(names(thresholds))) {
                    if (length(unique(nth)) == sum(is.element(nth, c("e", "c", "i")))) {
                        thresholds <- thresholds[match(c("e", "c", "i"), nth)]
                    }
                }
                thresholds <- as.vector(thresholds)
                thEX <- thresholds[1]
                thCR <- thresholds[2]
                thIN <- thresholds[3]
                if (logistic) {
                    if (thresholds[1] > thresholds[3]) {
                        thEX <- thresholds[3]
                        thIN <- thresholds[1]
                    }
                    y <- (x < thCR) + 1
                    fs <- 1/(1 + exp(c(-1, 1)[y] * (x - thCR) * log(idm/(1 - idm)) / (c(thIN, thEX)[y] - thCR)))
                    if (thresholds[1] > thresholds[3]) {
                        fs <- 1 - fs
                    }
                }
                else {
                    if (any(table(c(thEX, thCR, thIN)) > 1)) {
                        enter <- "\n"
                        if (is.element("enter", names(dots))) {
                            enter <- ifelse(isFALSE(dots$enter), "", dots$enter)
                        }
                        cat(enter)
                        warning(
                            simpleWarning("Some thresholds equal, that should not be equal.\n\n"),
                            .call = FALSE
                        )
                    }
                    if (above <= 0 | below <= 0) {
                        admisc::stopError(
                            "Arguments <above> and <below> should be positive.",
                            ... = ...
                        )
                    }
                    increasing <- TRUE
                    if (thIN < thCR & thCR < thEX) {
                        increasing <- FALSE
                    }
                    if (ecdf) {
                        ecdfx <- x[-which(x < min(thresholds))]
                        ecdfx <- ecdfx[-which(ecdfx > max(thresholds))]
                        Fn <- ecdf(ecdfx)
                    }
                    fs <- rep(NA, length(x))
                    for (i in seq(length(x))) {
                        if (increasing) {
                            if (x[i] < thEX | check.equal(x[i], thEX)) {
                                fs[i] <- 0
                            }
                            else if (x[i] < thCR | check.equal(x[i], thCR)) {
                                fs[i] <- (((thEX - x[i])/(thEX - thCR))^below)/2
                                if (ecdf) {
                                    fs[i] <- (Fn(x[i])/Fn(thCR))/2
                                }
                            }
                            else if (x[i] < thIN | check.equal(x[i], thIN)) {
                                fs[i] <- 1 - (((thIN - x[i])/(thIN - thCR))^above)/2
                                if (ecdf) {
                                    fs[i] <- 1 - ((1 - Fn(x[i]))/(1 - Fn(thCR)))/2
                                }
                            }
                            else {
                                fs[i] <- 1
                            }
                        }
                        else {
                            if (x[i] < thIN | check.equal(x[i], thIN)) {
                                fs[i] <- 1
                            }
                            else if (x[i] < thCR | check.equal(x[i], thCR)) {
                                fs[i] <- 1 - (((thIN - x[i])/(thIN - thCR))^above)/2
                                if (ecdf) {
                                    fs[i] <- 1 - (Fn(x[i])/Fn(thCR))/2
                                }
                            }
                            else if (x[i] < thEX | check.equal(x[i], thEX)) {
                                fs[i] <- (((thEX - x[i])/(thEX - thCR))^below)/2
                                if (ecdf) {
                                    fs[i] <- ((1 - Fn(x[i]))/(1 - Fn(thCR)))/2
                                }
                            }
                            else {
                                fs[i] <- 0
                            }
                        }
                    }
                }
            }
            else { 
                if (!is.null(nth)) {
                    if (length(unique(nth)) == sum(
                        is.element(
                            nth,
                            c("e1", "c1", "i1", "i2", "c2", "e2")
                        )
                    )) {
                        thresholds <- thresholds[
                            match(
                                c("e1", "c1", "i1", "i2", "c2", "e2"),
                                nth
                            )
                        ]
                    }
                }
                thresholds <- as.vector(thresholds)
                thEX1 <- thresholds[1]
                thCR1 <- thresholds[2]
                thIN1 <- thresholds[3]
                thIN2 <- thresholds[4]
                thCR2 <- thresholds[5]
                thEX2 <- thresholds[6]
                if (thCR1 < min(thEX1, thIN1) | thCR1 > max(thEX1, thIN1)) {
                    admisc::stopError(
                        "First crossover threshold not between first exclusion and inclusion thresholds.",
                        ... = ...
                    )
                }
                if (thCR2 < min(thEX2, thIN2) | thCR2 > max(thEX2, thIN2)) {
                    admisc::stopError(
                        "Second crossover threshold not between second exclusion and inclusion thresholds.",
                        ... = ...
                    )
                }
                somequal <- FALSE
                if (
                    any(
                        table(c(thEX1, thCR1, thIN1)) > 1
                    ) |
                    any(
                        table(c(thIN2, thCR2, thEX2)) > 1
                    ) |
                    thCR1 == thCR2
                ) {
                    somequal <- TRUE
                }
                increasing <- TRUE
                if (
                    thIN1 < thCR1 &
                    thCR1 < thEX1 &
                    thEX1 <= thEX2 &
                    thEX2 < thCR2 &
                    thCR2 < thIN2
                ) {
                    increasing <- FALSE
                }
                if (increasing) {
                    if (thEX1 == thEX2) {
                        somequal <- TRUE
                    }
                }
                else {
                    if (thIN1 == thIN2) {
                        somequal <- TRUE
                    }
                }
                if (somequal) {
                    admisc::stopError(
                        "Some thresholds equal, that should not be equal.",
                        ... = ...
                    )
                }
                if (above <= 0 | below <= 0) {
                    admisc::stopError(
                        "Arguments <above> and <below> should be positive.",
                        ... = ...
                    )
                }
                fs <- rep(NA, length(x))
                for (i in seq(length(x))) {
                    if (increasing) {
                        if (x[i] < thEX1 | check.equal(x[i], thEX1)) {
                            fs[i] <- 0
                        }
                        else if (x[i] < thCR1 | check.equal(x[i], thCR1)) {
                            fs[i] <- (((thEX1 - x[i])/(thEX1 - thCR1))^below)/2
                        }
                        else if (x[i] < thIN1) {
                            fs[i] <- 1 - (((thIN1 - x[i])/(thIN1 - thCR1))^above)/2
                        }
                        else if (x[i] < thIN2 | check.equal(x[i], thIN2)) {
                            fs[i] <- 1
                        }
                        else if (x[i] < thCR2 | check.equal(x[i], thCR2)) {
                            fs[i] <- 1 - (((thIN2 - x[i])/(thIN2 - thCR2))^above)/2
                        }
                        else if (x[i] < thEX2 | check.equal(x[i], thEX2)) {
                            fs[i] <- (((thEX2 - x[i])/(thEX2 - thCR2))^below)/2
                        }
                        else {
                            fs[i] <- 0
                        }
                    }
                    else {
                        if (x[i] < thIN1 | check.equal(x[i], thIN1)) {
                            fs[i] <- 1
                        }
                        else if (x[i] < thCR1 | check.equal(x[i], thCR1)) {
                            fs[i] <- 1 - (((thIN1 - x[i])/(thIN1 - thCR1))^above)/2
                        }
                        else if (x[i] < thEX1) {
                            fs[i] <- (((thEX1 - x[i])/(thEX1 - thCR1))^below)/2
                        }
                        else if (x[i] < thEX2 | check.equal(x[i], thEX2)) {
                            fs[i] <- 0
                        }
                        else if (x[i] < thCR2 | check.equal(x[i], thCR2)) {
                            fs[i] <- (((thEX2 - x[i])/(thEX2 - thCR2))^below)/2
                        }
                        else if (x[i] < thIN2 | check.equal(x[i], thIN2)) {
                            fs[i] <- 1 - (((thIN2 - x[i])/(thIN2 - thCR2))^above)/2
                        }
                        else {
                            fs[i] <- 1
                        }
                    }
                }
            }
            fs[fs < 0.0001] <- 0
            fs[fs > 0.9999] <- 1
            return(fs)
        }
        else if (method == "indirect") {
            thresholds <- sort(thresholds)
            values <- round(seq(0, 1, by = 1 / length(thresholds)), 3)
            y <- rep(0, length(x))
            for (i in seq(length(thresholds))) {
                y[x > thresholds[i]] = values[i + 1]
            }
            x[x == 0] <- 0.00001
            fracpol <- glm(
                y ~ log(x) + I(x^(1/2)) + I(x^1) + I(x^2),
                family = quasibinomial(logit)
            )
            fs <- round(
                unname(
                    predict(fracpol, type = "response")
                ),
                6
            )
            fs[fs < 0.0001] <- 0
            fs[fs > 0.9999] <- 1
            return(fs)
        }
        else if (method == "TFR") {
            E <- ecdf(x)
            return(pmax(0, (E(x) - E(1)) / (1 - E(1))))
        }
    }
}
