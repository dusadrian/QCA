`causalChain` <- function(
    data, ordering = NULL, strict = FALSE, pi.cons = 0, pi.depth = 0,
    sol.cons = 0, sol.cov = 1, sol.depth = 0, ...
) {

    metacall <- match.call(expand.dots = TRUE)
    # dots <- list(...)

    allargs <- as.list(metacall)[-1]
    allargs <- allargs[
        -which(
            is.element(
                c("data", "ordering", "strict"),
                names(allargs)
            )
        )
    ]

    if (
        is.element(c("pi.cons"), names(allargs)) && (
            is.null(pi.cons) ||
            !is.atomic(pi.cons) ||
            !is.numeric(pi.cons) ||
            length(pi.cons) != 1 ||
            is.na(pi.cons) ||
            pi.cons < 0 || pi.cons > 1
        )
    ) {
        admisc::stopError(
            "Argument <pi.cons> must be a single numeric value between 0 and 1.",
            ... = ...
        )
    }

    if (
        is.element(c("sol.cons"), names(allargs)) && (
            is.null(sol.cons) ||
            !is.atomic(sol.cons) ||
            !is.numeric(sol.cons) ||
            length(sol.cons) != 1 ||
            is.na(sol.cons) ||
            sol.cons < 0 || sol.cons > 1
        )
    ) {
        admisc::stopError(
            "Argument <sol.cons> must be a single numeric value between 0 and 1.",
            ... = ...
        )
    }

    if (
        is.element(c("sol.cov"), names(allargs)) && (
            is.null(sol.cov) ||
            !is.atomic(sol.cov) ||
            !is.numeric(sol.cov) ||
            length(sol.cov) != 1 ||
            is.na(sol.cov) ||
            sol.cov < 0 || sol.cov > 1
        )
    ) {
        admisc::stopError(
            "Argument <sol.cov> must be a single numeric value between 0 and 1.",
            ... = ...
        )
    }

    if (
        !is.element("incl.cut", names(allargs)) &&
        any(c(pi.cons, sol.cons) > 0)
    ) {
        allargs$incl.cut <- 0.5
    }

    verify.qca(data)
    noflevels  <- admisc::getInfo(data)$noflevels

    mv <- noflevels > 2
    names(noflevels) <- names(mv) <- colnames(data)


    if (methods::is(ordering, "character")) {
        
        ordering <- gsub("[[:space:]]", "", ordering)
        
        if (length(ordering) == 1) {
            ordering <- unlist(strsplit(ordering, split = "<"))
        }
        else {
            if (any(grepl("<", ordering))) {
                admisc::stopError(
                    "Causal ordering character \"<\" requires a single string.",
                    ... = ...
                )
            }
        }
        
        ordering <- lapply(ordering, admisc::splitstr)
    }

    if (length(allout <- unlist(ordering)) > 0) {
        if (length(setdiff(allout, colnames(data))) > 0) {
            admisc::stopError(
                "Some elements in the argument <ordering> not found in the data.",
                ... = ...
            )
        }
    }

    # allargs <- as.list(metacall)[-1]
    # allargs <- allargs[-which(is.element(names(allargs), c("ordering", "strict")))]

    allargs <- c(list(input = data), allargs)
    allargs$causalChain <- TRUE

    checkpos <- function(x, arg) {
        pos <- pmatch(names(allargs), arg)
        return(pos[!is.na(pos)])
    }

    pos <- checkpos(allargs, "include")
    if (length(pos) == 0) {
        allargs$include <- "?"
    }

    pos <- checkpos(allargs, "all.sol")
    if (length(pos) == 0) {
        allargs$all.sol <- TRUE
    }

    pos <- checkpos(allargs, "SA")
    if (length(pos) == 0) {
        allargs$SA <- FALSE
    }

    minimizeit <- function(allargs) {
        return(
            tryCatch(
                do.call("minimize", allargs),
                error = function(e) NA
            )
        )
        # tc <- tryCatch(do.call("minimize", allargs), error = function(e) e)
        # if (inherits(tc, "error")) return(NA)
        # return(tc)
    }

    allargs$enter <- ""

    minimize.list <- list()

    if (length(ordering) > 0) {
        
        if (any(table(unlist(ordering)) > 1)) {
            admisc::stopError(
                "Same condition(s) in multiple ordering levels.",
                ... = ...
            )
        }
        
        allcols <- colnames(data)
        if (length(restcols <- setdiff(allcols, unlist(ordering))) > 0) {
            ordering <- c(list(restcols), ordering)
        }
        
        for (i in seq(length(ordering))) {
            nextcols <- ordering[[i]]
            if (i == 1) {
                if (!strict & length(nextcols) > 1) {
                    for (j in seq(length(nextcols))) {
                        allargs$input <- subset(data, select = nextcols)
                        if (mv[nextcols[j]]) {
                            uniqv <- sort(unique(data[, nextcols[j]]))
                            for (v in seq(noflevels[nextcols[j]] - 1)) {
                                if (is.element(v, uniqv)) {
                                    allargs$outcome <- sprintf("%s[%s]", nextcols[j], v)
                                    # print(allargs$outcome)
                                    # return(allargs)
                                    minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                                }
                            }
                        }
                        else {
                            allargs$outcome <- nextcols[j]
                            # print(allargs$outcome)
                            # return(allargs)
                            minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                        }
                    }
                }
            }
            else {
                
                restcols <- unlist(ordering[seq(i - 1)])
                for (j in seq(length(nextcols))) {
                    
                    if (strict) {
                        allcols <- c(restcols, nextcols[j])
                    }
                    else {
                        allcols <- c(restcols, nextcols)
                    }
                    allcols <- allcols[order(match(allcols, colnames(data)))]
                    allargs$input <- subset(data, select = allcols)
                    
                    
                    if (mv[nextcols[j]]) {
                        uniqv <- sort(unique(data[, nextcols[j]]))
                        
                        for (v in seq(noflevels[nextcols[j]] - 1)) {
                            if (is.element(v, uniqv)) {
                                allargs$outcome <- sprintf("%s[%s]", nextcols[j], v)
                                # print(allargs$outcome)
                                # return(allargs)
                                minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                            }
                        }
                    }
                    else {
                        allargs$outcome <- nextcols[j]
                        # print(allargs$outcome)
                        # return(allargs)
                        minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                    }
                }
            }
        }
    }
    else {

        for (x in colnames(data)) {
            if (mv[x]) {
                uniqv <- sort(unique(data[, x]))
                
                for (v in seq(noflevels[x] - 1)) {
                    if (is.element(v, uniqv)) {
                        allargs$outcome <- sprintf("%s[%s]", x, v)
                        # print(allargs$outcome)
                        # return(allargs)
                        # if (FALSE) {
                        minimize.list[[allargs$outcome]] <- minimizeit(allargs)
                        # }
                    }
                }
            }
            else {
                allargs$outcome <- x
                # print(allargs$outcome)
                # return(allargs)
                minimize.list[[x]] <- minimizeit(allargs)
            }
            
        }
    }

    attr(minimize.list, "call") <- metacall
    return(structure(minimize.list, class = "QCA_chain"))
}
