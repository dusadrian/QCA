`generate` <- function(
    expression = "", snames = "", noflevels = NULL, ...
) {
    
    expression <- admisc::recreate(substitute(expression))
    snames <- admisc::recreate(substitute(snames))

    suf <- grepl("=>|->", expression)
    if (grepl("<=|<-", expression) & !suf) {
        admisc::stopError(
            "Invalid expression, relation should (also) indicate sufficiency.",
            ... = ...
        )
    }

    if (!is.null(noflevels)) {
        if (is.character(noflevels) & length(noflevels) == 1) {
            noflevels <- splitstr(noflevels)
        }
    }
    
    outcome <- ""
    if (suf) {
        necsuf <- grepl("<=>|<->", expression)
        expression <- unlist(
            strsplit(
                expression,
                split = ifelse(
                    necsuf,
                    "<=>|<->",
                    "->|=>"
                )
            )
        )
        outcome <- trimstr(expression[2])
        expression <- trimstr(expression[1])
    }

    if (!identical(snames, "")) {
        snames <- splitstr(snames)
    }

    trexp <- translate(expression, snames = snames, noflevels = noflevels)
    snames <- colnames(trexp)
    
    if (is.null(noflevels)) {
        noflevels <- rep(2, length(snames)) # to count the outcome too
    }

    tt <- as.data.frame(getMatrix(noflevels))

    pos <- expand(expression, snames = snames, noflevels = noflevels, implicants = TRUE) - 1

    mbase <- c(rev(cumprod(rev(noflevels))), 1)[-1]
    posrownms <- as.vector(as.matrix(pos) %*% mbase) + 1

    tt$OUT <- 0
    tt$OUT[as.vector(as.matrix(pos) %*% mbase) + 1] <- 1

    if (identical(outcome, "")) {
        if (any(nchar(snames) > 1)) {
            if (!is.element("OUT", snames)) {
                outcome <- "OUT"
            }
            else {
                outname <- paste(sample(LETTERS, 10), collapse = "")
            }
        }
        else {
            outname <- setdiff(c("O", "X", "Y", "Z", LETTERS, letters), snames)[1]
            # perhaps capture error if this has length zero, but highly unlikely though
        }
    }

    colnames(tt) <- c(snames, outcome)
    return(as.data.frame(tt))
}
