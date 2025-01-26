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

`nec` <- function(x) {
    !is.na(pmatch(x, "necessity"))
}
`suf` <- function(x) {
    !is.na(pmatch(x, "sufficiency"))
}
`colnms` <- function(mymat, rownms, tilde = FALSE) {
    apply(mymat, 1, function(x) {
        rownms1 <- rownms[x == 1]
        rownms[x == 1] <- if (tilde) paste0("~", rownms1) else tolower(rownms1)
        return(paste(rownms[x > 0], collapse = "*"))
    })
}
`colnms2` <- function(mymat, colnms, tilde = FALSE) {
    chars <- colnms[col(mymat)]
    lowerChars <- if (tilde) paste0("~", chars) else tolower(chars)
    chars <- ifelse(mymat==1L, lowerChars, chars)
    keep <- mymat > 0L
    charList <- split(chars[keep], row(chars)[keep])
    unlist(lapply(charList, paste, collapse = "*"))
}
`splitMainComponents2` <- function(expression) {
    expression <- gsub("[[:space:]]", "", expression)
    ind.char <- unlist(strsplit(expression, split=""))
    if (grepl("\\(", expression)) {
        open.brackets <- which(ind.char == "(")
        closed.brackets <- which(ind.char == ")")
        invalid <- ifelse(grepl("\\)", expression), length(open.brackets) != length(closed.brackets), FALSE)
        if (invalid) {
            cat("\n")
            stop("Invalid expression, open bracket \"(\" not closed with \")\".\n\n", call. = FALSE)
        }
        all.brackets <- sort(c(open.brackets, closed.brackets))
        if (length(all.brackets) > 2) {
            for (i in seq(3, length(all.brackets))) {
                if (all.brackets[i] - all.brackets[i - 1] == 1) {
                    open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                    closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                }
                if (all.brackets[i] - all.brackets[i - 1] == 2) {
                    if (ind.char[all.brackets[i] - 1] != "+") {
                        open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                        closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                    }
                }
            }
        }
        for (i in seq(length(open.brackets))) {
            plus.signs <- which(ind.char == "+")
            last.plus.sign <- plus.signs[plus.signs < open.brackets[i]]
            if (length(last.plus.sign) > 0) {
                open.brackets[i] <- max(last.plus.sign) + 1
            }
            else {
                if (1 == 1) { 
                    open.brackets[i] <- 1
                }
            }
            next.plus.sign <- plus.signs[plus.signs > closed.brackets[i]]
            if(length(next.plus.sign) > 0) {
                closed.brackets[i] <- min(next.plus.sign) - 1
            }
            else {
                closed.brackets[i] <- length(ind.char)
            }
        }
        big.list <- vector(mode="list", length = length(open.brackets) + 2)
        if (length(open.brackets) == 1) {
            if (open.brackets > 1) {
                big.list[[1]] <- paste(ind.char[seq(1, open.brackets - 2)], collapse = "")
            }
            nep <- min(which(unlist(lapply(big.list, is.null))))
            big.list[[nep]] <- paste(ind.char[seq(open.brackets, closed.brackets)], collapse = "")
            if (closed.brackets < length(ind.char)) {
                nep <- min(which(unlist(lapply(big.list, is.null))))
                big.list[[nep]] <- paste(ind.char[seq(closed.brackets + 2, length(ind.char))], collapse = "")
            }
        }
        else {
            for (i in seq(length(open.brackets))) {
                if (i == 1) {
                    if (open.brackets[1] > 1) {
                        big.list[[1]] <- paste(ind.char[seq(1, open.brackets[1] - 2)], collapse = "")
                    }
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                }
                else {
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    if (i == length(closed.brackets)) {
                        if (closed.brackets[i] < length(ind.char)) {
                            nep <- min(which(unlist(lapply(big.list, is.null))))
                            big.list[[nep]] <- paste(ind.char[seq(closed.brackets[i] + 2, length(ind.char))], collapse = "")
                        }
                    }
                }
            }
        }
        nulls <- unlist(lapply(big.list, is.null))
        if (any(nulls)) {
            big.list <- big.list[-which(nulls)]
        }
        big.list <- list(unlist(big.list))
    }
    else {
        big.list <- list(expression)
    }
    names(big.list) <- expression
    return(big.list)
}
`splitBrackets2` <- function(big.list) {
    big.list <- as.vector(unlist(big.list))
    result <- vector(mode="list", length = length(big.list))
    for (i in seq(length(big.list))) {
        result[[i]] <- trimstr(unlist(strsplit(unlist(strsplit(big.list[i], split="\\(")), split="\\)")), "*")
    }
    names(result) <- big.list
    return(result)
}
`splitPluses2` <- function(big.list) {
    return(lapply(big.list, function(x) {
        x2 <- lapply(x, function(y) {
            plus.split <- unlist(strsplit(y, "\\+"))
            return(plus.split[plus.split != ""])
        })
        names(x2) <- x
        return(x2)
    }))
}
`splitProducts` <- function(x, prod.split) {
    x <- as.vector(unlist(x))
    strsplit(x, split=prod.split)
}
`mvregexp` <- "\\[|\\]|\\{|\\}"
