`findmin` <- function(chart, ...) {

    dots <- list(...)
    verbose <- isTRUE(dots$verbose)

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

        gurobi <- !isFALSE(attr(chart, "gurobi")) &&
                !isFALSE(dots$gurobi) &&
                eval(parse(
                    text = "requireNamespace('gurobi', quietly = TRUE)"
                ))

        just_minima <- !isTRUE(attr(chart, "solind")) && !isTRUE(dots$solind)

        if (gurobi) {
            chart <- matrix(as.numeric(chart), nrow = nrow(chart))

            model <- list(
                A = chart,
                obj = rep(1, ncol(chart)),
                modelsense = "min",
                rhs = rep(1, nrow(chart)),
                sense = rep(">=", nrow(chart)),
                vtype = "B"
            )

            params <- list(
                OutputFlag = verbose * 1,
                LogToConsole = verbose * 1
            )

            tc <- admisc::tryCatchWEM(
                solution <- eval(parse(text = "gurobi::gurobi(model, params)"))$x
            )
            # solution <- gurobi::gurobi(model, params)$x

            if (!is.null(tc$error)) {
                gurobi <- FALSE
                if (cpi) {
                    message(
                        sprintf(
                            "%s, using lpSolve instead.",
                            ifelse(
                                grepl("license", tc$error),
                                "No valid Gurobi license found",
                                "Gurobi threw an error"
                            )
                        )
                    )
                }
            }
            gurobi <- is.null(tc)
            # if FALSE, Gurobi threw an error (i.e. not finding a valid license).
        }

        if (!gurobi) {
            solution <- lpSolve::lp(
                "min",
                rep(1, ncol(chart)),
                chart,
                ">=",
                1,
                int.vec = seq(nrow(chart)),
                all.bin = TRUE
            )$solution
        }

        if (just_minima | any(solution > 0 & solution < 1)) {
            # for large PI charts, lpSolve may return fractions of solutions, not helpful for integer PI charts
            solution <- as.integer(ceiling(sum(solution)))
        }
    }
    else {
        solution <- 0L
    }

    class(solution) <- c(class(solution), "QCA_findmin")
    return(solution)
}
