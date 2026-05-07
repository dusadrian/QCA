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
exact_matches <- 0L
accepted_gap1 <- 0L
lag_time_total <- 0
exp_time_total <- 0
scp_time_total <- 0
lp_time_total <- 0
for (k in seq_len(ntests)) {
    m <- matrix(runif(nr * nc) < density, nrow = nr)
    for (r in seq_len(nr)) {
        if (!any(m[r, ])) {
            m[r, sample.int(nc, 1)] <- TRUE
        }
    }
    chart <- matrix(as.logical(m), nrow = nr)
    lag_time <- system.time(
        lag <- .Call("C_findminLagrangianInfo", chart, PACKAGE = "QCA")
    )[["elapsed"]]
    exp_time <- system.time({
        exp_sol <- lag$solution
        exp_size <- lag$upper_bound
        lag_lb <- lag$lower_bound
        if (!is.na(lag_lb) && exp_size > 0 && exp_size <= ceiling(lag_lb - 1e-12) + 1L) {
            accepted_gap1 <<- accepted_gap1 + 1L
        }
        else {
            exp_sol <- .Call("C_findminScpInternal", chart, PACKAGE = "QCA")
            exp_size <- sum(exp_sol > 0)
        }
    })[["elapsed"]]
    scp_time <- system.time(
        scp <- .Call("C_findminScpInternal", chart, PACKAGE = "QCA")
    )[["elapsed"]]
    lp_time <- system.time(
        lp <- .Call("C_findminLpSolveInternal", chart, PACKAGE = "QCA")
    )[["elapsed"]]
    exp_exact <- sum(exp_sol > 0) == sum(lp > 0)
    if (exp_exact) {
        exact_matches <- exact_matches + 1L
    }
    lag_time_total <- lag_time_total + lag_time
    exp_time_total <- exp_time_total + exp_time
    scp_time_total <- scp_time_total + scp_time
    lp_time_total <- lp_time_total + lp_time
    cat(sprintf(
        paste0(
            "test=%d rows=%d cols=%d density=%.4f ",
            "lag_ub=%d lag_lb=%.3f exp=%d scp=%d lp=%d ",
            "exp_exact=%s lag_time=%.3f exp_time=%.3f scp_time=%.3f lp_time=%.3f\n"
        ),
        k, nr, nc, density,
        lag$upper_bound, lag$lower_bound, sum(exp_sol > 0), sum(scp > 0), sum(lp > 0),
        if (exp_exact) "TRUE" else "FALSE",
        lag_time, exp_time, scp_time, lp_time
    ))
}
cat(sprintf(
    paste0(
        "summary rows=%d cols=%d density=%.4f tests=%d seed=%d ",
        "exact_matches=%d accepted_gap1=%d ",
        "avg_lag_time=%.3f avg_exp_time=%.3f avg_scp_time=%.3f avg_lp_time=%.3f\n"
    ),
    nr, nc, density, ntests, seed,
    exact_matches, accepted_gap1,
    lag_time_total / ntests,
    exp_time_total / ntests,
    scp_time_total / ntests,
    lp_time_total / ntests
))
