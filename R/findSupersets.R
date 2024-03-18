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

`findSupersets` <- function (
    input, noflevels = NULL, ...
) {
    dots <- list(...)
        if (is.element("input.combs", names(dots))) {
            input <- dots$input.combs
        }
    if (!is.matrix(input)) {
        if (!is.vector(drop(input))) {
            admisc::stopError(
                "input must be either an solution-space matrix or a vector of row numbers.",
                ... = ...
            )
        }
        else {
            if (any(input > prod(noflevels))) {
                admisc::stopError(
                    paste(
                        "Some line numbers do not belong in the solution-space for",
                        length(noflevels),
                        "causal conditions."
                    ),
                    ... = ...
                )
            }
            input <- getRow(input, noflevels)
        }
    }
    mbase <- rev(c(1, cumprod(rev(noflevels))))[-1]
    allcombn <- t(createMatrix(rep(2, length(noflevels)))[-1, ])
    primes <- sort.int(
        unique(
            as.vector(
                apply(input, 1, function(x) (x*mbase) %*% allcombn + 1)
            )
        )
    )
    if (primes[1] == 1) {
        return(primes[-1])
    }
    else {
        return(primes)
    }
}
