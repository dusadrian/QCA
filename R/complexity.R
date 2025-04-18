# Copyright (c) 2016 - 2025, Adrian Dusa
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

complexity <- function(
    n, layers = NULL, noflevels = NULL, ...
) {
    if (!is.numeric(n)) {
        admisc::stopError("Argument <n> should be numeric.", ... = ...)
    }
    if (length(n) != 1L) {
        admisc::stopError("Argument <n> should be a scalar of length 1.", ... = ...)
    }
    if (n < 0) {
        admisc::stopError("Argument <n> should be positive.", ... = ...)
    }
    if (is.null(noflevels)) noflevels <- rep(2, n)
    if (is.null(layers)) layers <- seq(n)
    if (any(layers > n)) {
        admisc::stopError("Argument <layers> cannot be greater than <n>.", ... = ...)
    }
    sumk <- .Call("C_omplexity", list(as.integer(n), as.integer(layers), as.integer(noflevels)), PACKAGE = "QCA")
    sumk[sumk < 0] <- Inf
    return(sumk)
    sumk <- rep(0, length(layers))
    for (i in seq(length(layers))) {
        sumk[i] <- sum(apply(admisc::combnk(n, layers[i]), 2, function(x) {
            prod(noflevels[x])
        }))
    }
}
