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

`generate` <- function(
    expression = "", snames = "", noflevels = NULL, ...
) {
    expression <- admisc::recreate(substitute(expression))
    snames <- admisc::recreate(substitute(snames))
    suf <- grepl("=>|->", expression)
    if (grepl("<=|<-", expression) & !suf) {
        admisc::stopError(
            "Invalid expression, relation should (also) indicate sufficiency.",
            ... = ...
        )
    }
    if (!is.null(noflevels)) {
        if (is.character(noflevels) & length(noflevels) == 1) {
            noflevels <- splitstr(noflevels)
        }
    }
    outcome <- ""
    if (suf) {
        necsuf <- grepl("<=>|<->", expression)
        expression <- unlist(
            strsplit(
                expression,
                split = ifelse(
                    necsuf,
                    "<=>|<->",
                    "->|=>"
                )
            )
        )
        outcome <- trimstr(expression[2])
        expression <- trimstr(expression[1])
    }
    if (!identical(snames, "")) {
        snames <- splitstr(snames)
    }
    trexp <- translate(expression, snames = snames, noflevels = noflevels)
    snames <- colnames(trexp)
    if (is.null(noflevels)) {
        noflevels <- rep(2, length(snames)) 
    }
    tt <- as.data.frame(getMatrix(noflevels))
    pos <- expand(expression, snames = snames, noflevels = noflevels, implicants = TRUE) - 1
    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    posrownms <- as.vector(as.matrix(pos) %*% mbase) + 1
    tt$OUT <- 0
    tt$OUT[as.vector(as.matrix(pos) %*% mbase) + 1] <- 1
    if (identical(outcome, "")) {
        if (any(nchar(snames) > 1)) {
            if (!is.element("OUT", snames)) {
                outcome <- "OUT"
            }
            else {
                outname <- paste(sample(LETTERS, 10), collapse = "")
            }
        }
        else {
            outname <- setdiff(c("O", "X", "Y", "Z", LETTERS, letters), snames)[1]
        }
    }
    colnames(tt) <- c(snames, outcome)
    return(as.data.frame(tt))
}
