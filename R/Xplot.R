`Xplot` <- function(
    x, jitter = FALSE, at = pretty(x), ...
) {

    dots <- list(...)
    funargs <- unlist(lapply(match.call(), deparse)[-1])

    xname <- admisc::getName(funargs[["x"]])

    linex <- 1.75
    jitfactor <- 0.5
    jitamount <- 0.5
    cexpoints <- 1
    cexaxis <- 0.8
    pch <- 21
    bgpoints <- NA

    if (length(testarg <- which(names(dots) == "line")) > 0) {
        linex <- dots$line
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "factor")) > 0) {
        jitfactor <- dots$factor
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "amount")) > 0) {
        jitamount <- dots$amount
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "cex")) > 0) {
        cexpoints <- dots$cex
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "cex.axis")) > 0) {
        cexaxis <- dots$cex.axis
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "pch")) > 0) {
        pch <- dots$pch
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "bg")) > 0) {
        bgpoints <- dots$bg
        dots <- dots[-testarg]
    }

    if (length(testarg <- which(names(dots) == "xlab")) > 0) {
        xname <- dots$xlab
        dots <- dots[-testarg]
    }


    y <- rep(1, length(x))

    if (jitter) {

        y <- jitter(y, jitfactor, jitamount)
    }

    toplot <- list(as.name("plot"), x, y)

    toplot$type <- "n"
    if (!is.null(at)) {
        toplot$xlim <- range(at)
    }
    toplot$ylim <- c(0, 2)
    toplot$xlab <- ""
    toplot$ylab <- ""
    toplot$axes <- FALSE

    if (length(dots) > 0) {
        toplot <- c(toplot, dots)
    }

    par(mar = c(ifelse(xname == "", 2, 3), 0.3, 0, 0))

    suppressWarnings(eval(as.call(toplot)))

    # box()
    axis(1, at = at, cex.axis = cexaxis)
    title(xlab = xname, cex.lab = cexaxis + 0.1, font.lab = 2, line = linex)

    plotpoints <- list(as.name("points"), x, y, pch = pch, cex = cexpoints, bg = bgpoints)

    suppressWarnings(eval(as.call(c(plotpoints, dots))))

}
