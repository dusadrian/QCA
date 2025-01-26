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

`retention` <- function(
    data, outcome = "", conditions = "", incl.cut = 1, n.cut = 1,
    type = "corruption", dependent = TRUE, p.pert = 0.5, n.pert = 1, ...
) {
    outcome <- admisc::recreate(substitute(outcome))
    conditions <- admisc::recreate(substitute(conditions))
    type <- admisc::recreate(substitute(type))
    if (
        any(
            grepl(
                "\\{|\\[",
                c(outcome, conditions)
            )
        )
    ) {
        admisc::stopError(
            "Only binary data allowed.", ... = ...
        )
    }
    nms <- colnames(data)
    if (identical(conditions, "")) {
        conditions <- setdiff(nms, outcome)
    }
    else {
        conditions <- admisc::splitstr(conditions)
        if (length(conditions) == 1 & any(grepl(":", conditions))) {
            cs <- unlist(strsplit(conditions, split = ":"))
            if (!all(is.element(conditions, nms))) {
                admisc::stopError(
                    "Conditions from sequence not found in the data.", ... = ...
                )
            }
            conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
        }
        if (!all(is.element(conditions, nms))) {
            admisc::stopError(
                "Conditions not found in the data.", ... = ...
            )
        }
    }
    if (!is.element(outcome, nms)) {
        admisc::stopError(
            "Outcome not found in the data.", ... = ...
        )
    }
    data <- data[, c(conditions, outcome)]
    udata <- unique(data[, conditions])
    rownames(udata) <- seq(nrow(udata))
    cpos <- cneg <- rep(0, nrow(udata))
    for (i in seq(nrow(udata))) {
        for (j in seq(nrow(data))) {
            if (all(udata[i, ] == data[j, conditions])) {
                if (data[j, outcome] == 1) {
                    cpos[i] <- cpos[i] + 1
                }
                else if (data[j, outcome] == 0) {
                    cneg[i] <- cneg[i] + 1
                }
            }
        }
    }
    total <- cpos + cneg
    udata <- udata[total >= n.cut, , drop = FALSE]
    cpos  <- cpos[total >= n.cut]
    cneg  <- cneg[total >= n.cut]
    total <- total[total >= n.cut]
    calculatePairs <- function(x, n.pert, type = "deletion") {
        pairsxl <- admisc::combnk(nrow(udata), min(x, nrow(udata)))
        nofsetsxl <- 0
        for (j in seq(ncol(pairsxl))) {
            cposneg <- NULL
            for (l in seq(min(x, nrow(pairsxl)))) {
                cposneg <- c(cposneg, cbind(cpos, cneg)[pairsxl[l, j], ])
            }
            allpairs <- createMatrix(cposneg + 2) - 1
            allpairs <- allpairs[apply(allpairs, 1, function(y) all(y >= 0)), , drop = FALSE]
            linesubset <- rep(TRUE, nrow(allpairs))
            for (l in seq(min(x, nrow(pairsxl)))) {
                linesubset <- linesubset & rowSums(allpairs[, l*2 - c(1, 0)]) >= 1
            }
            allpairs <- allpairs[linesubset & rowSums(allpairs) <= n.pert, , drop = FALSE]
            for (i in seq(nrow(allpairs))) {
                lchanges <- rep(FALSE, min(x, nrow(pairsxl)))
                for (l in seq(min(x, nrow(pairsxl)))) {
                    initially <- cpos[pairsxl[l, j]]/total[pairsxl[l, j]]
                    if (type == "deletion") {
                        newtotaless <- total[pairsxl[l, j]] - allpairs[i, l*2 - 1]
                        after <- (cpos[pairsxl[l, j]] - allpairs[i, l*2 - 1])/newtotaless
                        lchanges[l] <- ((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                                       ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)
                    }
                    else if (type == "corruption") {
                        after <- (cpos[pairsxl[l, j]] + allpairs[i, l*2] - allpairs[i, l*2 - 1])/total[pairsxl[l, j]]
                        lchanges[l] <- (initially >= incl.cut & after <  incl.cut) |
                                       (initially <  incl.cut & after >= incl.cut)
                    }
                }
                if (all(lchanges)) {
                    combs <- 1
                    for (l in seq(min(x, nrow(pairsxl)))) {
                        combs <- combs*choose(cpos[pairsxl[l, j]], allpairs[i, l*2 - 1])
                        combs <- combs*choose(cneg[pairsxl[l, j]], allpairs[i, l*2])
                    }
                    nofsetsxl <- nofsetsxl + combs*choose(nrow(data) - sum(cposneg), n.pert - sum(allpairs[i, ]))
                }
            }
        }
        return(nofsetsxl)
    }
    if (dependent) {
        nofsets <- 0
        totalsets <- choose(nrow(data), n.pert)
        for (i in seq(n.pert)) {
            if (nofsets != totalsets) {
                nofsetsxl <- calculatePairs(i, n.pert, type = type)
                nofsets <- nofsets + ifelse(i %% 2 == 1, nofsetsxl, -1*nofsetsxl)
            }
        }
        return(as.vector(1 - nofsets/totalsets))
    }
    else {
        pfinal <- 1
        if (type == "deletion") {
            for (l in seq(nrow(udata))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq(nrow(allpairs))) {
                    newtotaless <- total[l] - allpairs[i, 1] - allpairs[i, 2]
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1])/newtotaless
                    if (((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                        ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        else if (type == "corruption") {
            for (l in seq(nrow(udata))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq(nrow(allpairs))) {
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1] + allpairs[i, 2])/total[l]
                    if ((initially >= incl.cut & after < incl.cut) | (initially < incl.cut & after >= incl.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        return(pfinal)
    }
}
