`solveChart` <- function(
    chart, row.dom = FALSE, all.sol = FALSE, depth = NULL, max.comb = 0,
    first.min = FALSE, ...
) {
    
    if (!is.logical(chart) && length(setdiff(chart, 0:1)) > 0) {
        admisc::stopError(
            "Use a logical, T/F matrix. See makeChart()'s output.", ... = ...
        )
    }
    
    dots <- list(...)
    
    if (is.element("min.dis", names(dots))) {
        if (is.logical(dots$min.dis)) {
            all.sol <- !dots$min.dis
        }
    }

    if (all.sol) {
        row.dom <- FALSE
    }

    row.numbers <- seq(nrow(chart))
    
    if (row.dom) {
        row.numbers <- rowDominance(chart)
        chart <- chart[row.numbers, ]
    }


    # indexes <- which(chart, arr.ind = TRUE)

    
    foundm <- findmin(chart, ... = ...) # , quick = FALSE
    
    
    if (foundm == 0) {
        admisc::stopError(
            "The PI chart cannot be solved.", ... = ...
        )
    }
    
    if (is.null(depth)) depth <- 0L

    # return(list(matrix(as.logical(chart), nrow = nrow(chart)), all.sol, as.integer(depth), as.integer(foundm), as.integer(max.comb)))
    output <- .Call(
        "C_solveChart",
        matrix(as.logical(chart),
        nrow = nrow(chart)),
        all.sol,
        as.integer(depth),
        as.integer(foundm),
        max.comb,
        first.min,
        PACKAGE = "QCA"
    )
    

    if (output[[2]]) {
        warning(
            simpleWarning(
                "The PI chart is exceedingly complex, solution(s) not guaranteed to be exhaustive.\n\n"
            )
        )
    }

    output <- output[[1]]

    # just in case row dominance was applied
    output[output == 0] <- NA
    output <- matrix(as.integer(row.numbers[output]), nrow = nrow(output))
    output[is.na(output)] <- 0L
    
    return(output)
}
