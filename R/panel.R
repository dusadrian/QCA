`as.panel` <- function(
    x, row.names, ...
) {
    
    if (!missing(row.names)) {
        if (is.character(row.names)) {
            if (length(row.names) == 1L) {
                rowvar <- (1L:ncol(x))[match(colnames(x), row.names, 0L) == 1L]
                row.names <- x[[rowvar]]
                x <- x[-rowvar]
            }
        }
        else if (is.numeric(row.names) && length(row.names) == 1L) {
            rowvar <- row.names
            row.names <- x[[rowvar]]
            x <- x[-rowvar]
        }
        else {
            admisc::stopError(
                "invalid 'row.names' specification.", ... = ...
            )
        }
        
        if (is.object(row.names) || !(is.integer(row.names))) {
            row.names <- as.character(row.names)
        }
        
        if (anyNA(row.names)) {
            admisc::stopError(
                "missing values in 'row.names' are not allowed.", ... = ...
            )
        }
        
        attr(x, "row.names") <- row.names
    }
    
    if (!is.data.frame(x)) {
        x <- as.data.frame(x, stringsAsFactors = FALSE)
    }
    
    structure(x, class = c("panel", "data.frame"))
}



`[.panel` <- function(x, i, j, ...) {
    funargs <- unlist(lapply(match.call(), deparse)[-1])
    drop <- TRUE
    if (any(names(funargs) == "drop")) {
        drop <- as.logical(funargs["drop"])
    }
    
    rownms <- row.names(x)
    class(x) <- "data.frame"
    classes <- lapply(x, class)
    clevels <- lapply(x, levels)
    cordered <- lapply(x, is.ordered)
    x <- eval(
        parse(
            text = sprintf(
                "x[%s, %s, drop = %s]",
                funargs["i"],
                funargs["j"],
                drop
            )
        )
    )
    
    if (!is.null(dim(x))) {
        x <- as.matrix(x)
        rownms <- eval(
            parse(
                text = sprintf(
                    "rownms[%s]",
                    funargs["i"]
                )
            )
        )
        row.names(x) <- rownms
        x <- rebuild(
            as.data.frame(x, stringsAsFactors = FALSE),
            classes[colnames(x)],
            clevels[colnames(x)],
            cordered[colnames(x)]
        )
        class(x) <- c("QCA_panel", "data.frame")
    }
    
    return(x)
}



`row.names<-.panel` <- function(x, value) {
    
    classes <- lapply(x, class)
    clevels <- lapply(x, levels)
    cordered <- lapply(x, is.ordered)
    x <- as.matrix(x)
    admisc::setRownames(x, value)
    x <- rebuild(
        as.data.frame(x, stringsAsFactors = FALSE),
        classes,
        clevels,
        cordered
    )
    
    class(x) <- c("QCA_panel", "data.frame")
    return(x)
}



`rebuild` <- function(x, classes, clevels, cordered) {
    for (i in seq(ncol(x))) {
        x[[i]] <- if (is.element("factor", classes[[i]])){
            factor(x[[i]], levels = clevels[[i]], ordered = cordered[[i]])
        }
        else if (is.element("Date", classes[[i]])) {
            as.Date(x[[i]])
        }
        else if (is.element("POSIXct", classes[[i]])) {
            as.POSIXct(x[[i]])
        }
        else {
            methods::as(x[[i]], classes[[i]])
        }
    }
    
    return(x)
}
