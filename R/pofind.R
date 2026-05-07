`pofind` <- function(
    data = NULL, outcome = "", conditions = "", relation = "necessity",
    use.labels = FALSE, ...
) {
    
    if (missing(data)) {
        admisc::stopError(
            "Data is missing.", ... = ...
        )
    }

    dots <- list(...)

    if (isTRUE(dots$categorical)) { # backwards compatibility
        use.labels <- TRUE
        dots$categorical <- NULL
    }

    funargs <- lapply(match.call(), deparse)
    outcome <- admisc::recreate(substitute(outcome), colnames(data))
    conditions <- admisc::recreate(substitute(conditions), colnames(data))
    enter <- if (is.element("enter", names(dots))) "" else "\n" # internal
    
    if (identical(conditions, "")) {
        conditions <- setdiff(colnames(data), admisc::notilde(outcome))
    }
    else {
        if (is.character(conditions) & length(conditions) == 1) {
            conditions <- admisc::splitstr(conditions)
            if (length(conditions) == 1) {
                if (grepl(":", conditions)) {
                    # they are already verified as existing columns, in verify.tt()
                    nms <- colnames(data)
                    cs <- unlist(strsplit(conditions, split = ":"))
                    if (!all(is.element(cs, nms))) {
                        admisc::stopError(
                            "Inexisting condition(s) in the sequence.", ... = ...
                        )
                    }
                    conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
                }
            }
        }
    }
    
    if (identical(outcome, "")) {
        admisc::stopError(
            "The outcome is missing.", ... = ...
        )
    }
    
        
    if (is.matrix(data)) {
        data <- as.data.frame(data)
    }
    
    verify.qca(data)
    
    for (i in seq(ncol(data))) {
        if (!is.numeric(data[, i])) {
            if (admisc::possibleNumeric(data[, i])) {
                data[, i] <- admisc::asNumeric(data[, i])
            }
        }
    }
    origoutcome <- outcome

    if (grepl("\\{", outcome)) {
        outcome <- admisc::curlyBrackets(outcome, outside = TRUE)
    }
    else {
        admisc::squareBrackets(outcome, outside = TRUE)
    }

    outcome <- admisc::notilde(outcome)
    if (!is.element(outcome, colnames(data))) {
        admisc::stopError(
            "Outcome not found in the data.", ... = ...
        )
    }
    
    if (identical(conditions, "")) {
        conditions <- setdiff(colnames(data), outcome)
    }
    else {
        conditions <- admisc::splitstr(conditions)

        verify.data(data, outcome, conditions)
        
        if (length(conditions) == 1) {
            if (grepl(":", conditions)) {
                # they are already verified as existing columns
                nms <- colnames(data)
                cs <- unlist(strsplit(conditions, split = ":"))
                conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
            }
        }
    }
    
    data <- data[, c(conditions, outcome)]

    infodata <- admisc::getInfo(data[, conditions, drop = FALSE])
    noflevels <- infodata$noflevels

    if (any(noflevels > 2)) { # multivalue conditions
        expression <- paste(unlist(lapply(seq(length(conditions)), function(x) {
            values <- sort(unique(data[, conditions[x]]))
            return(
                paste(conditions[x], "[", values, "]", sep = "")
            )
        })), collapse = "+")
    }
    else {
        negconditions <- paste("~", conditions, sep = "")
        expression <- paste(
            negconditions,
            conditions,
            sep = "+",
            collapse = "+"
        )
    }

    # return(list(setms = paste(c(conditions, negconditions), collapse = "+"),
    #                 outcome = origoutcome,
    #                 data = data,
    #                 relation = relation,
    #                 ... = ...))
    
    pofargs <- list(
        setms = expression,
        outcome = origoutcome,
        data = data,
        relation = relation,
        use.labels = use.labels,
        ... = ...
    )
    
    
    result <- do.call(pof, pofargs)
    result$incl.cov <- result$incl.cov[-nrow(result$incl.cov), , drop = FALSE]
    result$options$setms <- result$options$setms[, -ncol(result$options$setms), drop = FALSE]
    if (is.element("covU", colnames(result$incl.cov))) {
        result$incl.cov <- result$incl.cov[, setdiff(colnames(result$incl.cov), "covU")]
    }
    
    rownames(result$incl.cov)[seq(length(conditions))] <- paste(
        "",
        rownames(result$incl.cov)[
            seq(length(conditions))
        ]
    )
    
    result$options$data <- funargs$data
    
    return(result)
}
