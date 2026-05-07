# Copyright (c) 2016 - 2026, Adrian Dusa
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

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default) {
    prefix <- paste0("--", name, "=")
    match <- args[startsWith(args, prefix)]
    if (length(match) == 0) {
        return(default)
    }
    sub(prefix, "", match[[1]], fixed = TRUE)
}
nr <- as.integer(get_arg("rows", "20"))
nc <- as.integer(get_arg("cols", "1000"))
density <- as.numeric(get_arg("density", "0.03"))
seed <- as.integer(get_arg("seed", "1"))
set.seed(seed)
suppressPackageStartupMessages(library(QCA))
m <- matrix(runif(nr * nc) < density, nrow = nr)
for (r in seq_len(nr)) {
    if (!any(m[r, ])) {
        m[r, sample.int(nc, 1)] <- TRUE
    }
}
chart <- matrix(as.logical(m), nrow = nr)
.Call("C_resetScpProfile", PACKAGE = "QCA")
scp <- .Call("C_findminScpInternal", chart, PACKAGE = "QCA")
profile <- .Call("C_getScpProfile", PACKAGE = "QCA")
lp <- .Call("C_findminLpSolveInternal", chart, PACKAGE = "QCA")
cat(sprintf("scp=%d lp=%d\n", sum(scp), sum(lp)))
print(profile)
