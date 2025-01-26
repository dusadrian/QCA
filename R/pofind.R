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

`pofind` <- function(
    data = NULL, outcome = "", conditions = "", relation = "necessity",
    use.labels = FALSE, ...
) {
    if (missing(data)) {
        admisc::stopError(
            "Data is missing.", ... = ...
        )
    }
    dots <- list(...)
    if (isTRUE(dots$categorical)) { 
        use.labels <- TRUE
        dots$categorical <- NULL
    }
    funargs <- lapply(match.call(), deparse)
    outcome <- admisc::recreate(substitute(outcome), colnames(data))
    conditions <- admisc::recreate(substitute(conditions), colnames(data))
    enter <- if (is.element("enter", names(dots))) "" else "\n" 
    if (identical(conditions, "")) {
        conditions <- setdiff(colnames(data), admisc::notilde(outcome))
    }
    else {
        if (is.character(conditions) & length(conditions) == 1) {
            conditions <- admisc::splitstr(conditions)
            if (length(conditions) == 1) {
                if (grepl(":", conditions)) {
                    nms <- colnames(data)
                    cs <- unlist(strsplit(conditions, split = ":"))
                    if (!all(is.element(cs, nms))) {
                        admisc::stopError(
                            "Inexisting condition(s) in the sequence.", ... = ...
                        )
                    }
                    conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
                }
            }
        }
    }
    if (identical(outcome, "")) {
        admisc::stopError(
            "The outcome is missing.", ... = ...
        )
    }
    if (is.matrix(data)) {
        data <- as.data.frame(data)
    }
    verify.qca(data)
    for (i in seq(ncol(data))) {
        if (!is.numeric(data[, i])) {
            if (admisc::possibleNumeric(data[, i])) {
                data[, i] <- admisc::asNumeric(data[, i])
            }
        }
    }
    origoutcome <- outcome
    if (grepl("\\{", outcome)) {
        outcome <- admisc::curlyBrackets(outcome, outside = TRUE)
    }
    else {
        admisc::squareBrackets(outcome, outside = TRUE)
    }
    outcome <- admisc::notilde(outcome)
    if (!is.element(outcome, colnames(data))) {
        admisc::stopError(
            "Outcome not found in the data.", ... = ...
        )
    }
    if (identical(conditions, "")) {
        conditions <- setdiff(colnames(data), outcome)
    }
    else {
        conditions <- admisc::splitstr(conditions)
        verify.data(data, outcome, conditions)
        if (length(conditions) == 1) {
            if (grepl(":", conditions)) {
                nms <- colnames(data)
                cs <- unlist(strsplit(conditions, split = ":"))
                conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
            }
        }
    }
    data <- data[, c(conditions, outcome)]
    infodata <- admisc::getInfo(data[, conditions, drop = FALSE])
    noflevels <- infodata$noflevels
    if (any(noflevels > 2)) { 
        expression <- paste(unlist(lapply(seq(length(conditions)), function(x) {
            values <- sort(unique(data[, conditions[x]]))
            return(
                paste(conditions[x], "[", values, "]", sep = "")
            )
        })), collapse = "+")
    }
    else {
        negconditions <- paste("~", conditions, sep = "")
        expression <- paste(
            negconditions,
            conditions,
            sep = "+",
            collapse = "+"
        )
    }
    pofargs <- list(
        setms = expression,
        outcome = origoutcome,
        data = data,
        relation = relation,
        use.labels = use.labels,
        ... = ...
    )
    result <- do.call(pof, pofargs)
    result$incl.cov <- result$incl.cov[-nrow(result$incl.cov), , drop = FALSE]
    result$options$setms <- result$options$setms[, -ncol(result$options$setms), drop = FALSE]
    if (is.element("covU", colnames(result$incl.cov))) {
        result$incl.cov <- result$incl.cov[, setdiff(colnames(result$incl.cov), "covU")]
    }
    rownames(result$incl.cov)[seq(length(conditions))] <- paste(
        "",
        rownames(result$incl.cov)[
            seq(length(conditions))
        ]
    )
    result$options$data <- funargs$data
    return(result)
}
