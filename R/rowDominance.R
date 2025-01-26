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

`rowDominance` <- function(
    mtrx
) {
    rownames(mtrx) <- seq(nrow(mtrx))
    mtrx.copy <- mtrx
    sums <- rowSums(mtrx)
    mtrx <- mtrx[order(sums, decreasing=TRUE), , drop=FALSE]
    sums <- sort(sums, decreasing=TRUE)
    line.no <- 1
    while(line.no < nrow(mtrx)) {
        less <- sums < sums[line.no]
        if (any(less)) {
            aa <- apply(mtrx[less, , drop=FALSE], 1, function(x) {all(mtrx[line.no, x])})
            mtrx <- rbind(mtrx[!less, , drop=FALSE], mtrx[less, , drop=FALSE][!aa, , drop=FALSE])
            sums <- rowSums(mtrx)
            line.no <- line.no + 1
        }
        else {
            break
        }
    }
    return(sort(match(rownames(mtrx), rownames(mtrx.copy))))
}
