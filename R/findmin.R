`findmin` <- function(chart, type = c("exact", "lagrangian"), ...) {

    dots <- list(...)
    verbose <- isTRUE(dots$verbose)

    if (missing(type) || is.null(type)) {
        type <- attr(chart, "type")
    }

    if (is.null(type)) {
        type <- "exact"
    }

    type <- match.arg(type, c("exact", "lagrangian"))

    if (!methods::is(chart, "QCA_pic")) {

        if (!is.matrix(chart) | (!is.logical(chart) & length(setdiff(chart, 0:1)) > 0)) {
            admisc::stopError(
                "Use a logical, TRUE/FALSE matrix. See makeChart()'s output.",
                ... = ...
            )
        }
    }

    cpi <- !is.null(attr(chart, "C_PI"))

    if (!cpi) {
        chart <- t(chart)
    }

    if (all(colSums(chart) > 0)) {

        just_minima <- !isTRUE(attr(chart, "solind")) && !isTRUE(dots$solind)
        gurobi <- !isFALSE(attr(chart, "gurobi")) && !isFALSE(dots$gurobi)
        solution <- NULL

        if (identical(type, "lagrangian")) {
            solution <- .Call(
                "C_findminLagrangian",
                matrix(as.logical(chart), nrow = nrow(chart)),
                PACKAGE = "QCA"
            )
        } else if (gurobi) {
            native_gurobi_available <- getOption("native.gurobi.available", NULL)
            if (is.null(native_gurobi_available)) {
                native_gurobi_available <- isTRUE(.Call("C_gurobiRuntimeAvailable", PACKAGE = "QCA"))
                options(native.gurobi.available = native_gurobi_available)
            }

            if (isTRUE(native_gurobi_available)) {
                solution <- .Call(
                    "C_findminExact",
                    matrix(as.logical(chart), nrow = nrow(chart)),
                    PACKAGE = "QCA"
                )
            } else if (verbose) {
                message("Gurobi not available, falling back to the internal lp_solve solver.")
            }
        }

        if (is.null(solution)) {
            solution <- .Call(
                "C_findminLpSolveInternal",
                matrix(as.logical(chart), nrow = nrow(chart)),
                PACKAGE = "QCA"
            )
        }
        

        if (just_minima) {
            solution <- as.integer(ceiling(sum(solution)))
        }
    }
    else {
        solution <- 0L
    }

    class(solution) <- c(class(solution), "QCA_findmin")
    return(solution)
}
