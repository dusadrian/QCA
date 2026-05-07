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
ntests <- as.integer(get_arg("tests", "25"))
nr <- as.integer(get_arg("rows", "20"))
nc <- as.integer(get_arg("cols", "1000"))
density <- as.numeric(get_arg("density", "0.03"))
seed <- as.integer(get_arg("seed", "1"))
set.seed(seed)
suppressPackageStartupMessages(library(QCA))
for (k in seq_len(ntests)) {
    m <- matrix(runif(nr * nc) < density, nrow = nr)
    for (r in seq_len(nr)) {
        if (!any(m[r, ])) {
            m[r, sample.int(nc, 1)] <- TRUE
        }
    }
    chart <- matrix(as.logical(m), nrow = nr)
    scp_time <- system.time(
        scp <- .Call("C_findminScpInternal", chart, PACKAGE = "QCA")
    )[["elapsed"]]
    lp_time <- system.time(
        lp <- .Call("C_findminLpSolveInternal", chart, PACKAGE = "QCA")
    )[["elapsed"]]
    cat(sprintf(
        "test=%d rows=%d cols=%d density=%.4f scp=%d lp=%d scp_time=%.3f lp_time=%.3f\n",
        k, nr, nc, density, sum(scp), sum(lp), scp_time, lp_time
    ))
    if (sum(scp) != sum(lp)) {
        cat("Mismatch detected.\n")
        cat("scp_idx=", paste(which(scp > 0), collapse = ","), "\n", sep = "")
        cat("lp_idx=", paste(which(lp > 0), collapse = ","), "\n", sep = "")
        stop("SCP solver validation failed.", call. = FALSE)
    }
}
cat(sprintf(
    "validated=%d rows=%d cols=%d density=%.4f seed=%d\n",
    ntests, nr, nc, density, seed
))
