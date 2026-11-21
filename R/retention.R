`retention` <- function(
    data, outcome = "", conditions = "", incl.cut = 1, n.cut = 1,
    type = "corruption", dependent = TRUE, p.pert = 0.5, n.pert = 1, ...
) {
    
    outcome <- admisc::recreate(substitute(outcome))
    conditions <- admisc::recreate(substitute(conditions))
    type <- admisc::recreate(substitute(type))

    if (
        any(
            grepl(
                "\\{|\\[",
                c(outcome, conditions)
            )
        )
    ) {
        admisc::stopError(
            "Only binary data allowed.", ... = ...
        )
    }

    nms <- colnames(data)

    if (identical(conditions, "")) {
        conditions <- setdiff(nms, outcome)
    }
    else {
        conditions <- admisc::splitstr(conditions)
        
        if (length(conditions) == 1 & any(grepl(":", conditions))) {
            cs <- unlist(strsplit(conditions, split = ":"))
            if (!all(is.element(cs, nms))) {
                admisc::stopError(
                    "Conditions from sequence not found in the data.", ... = ...
                )
            }
            conditions <- nms[seq(which(nms == cs[1]), which(nms == cs[2]))]
        }

        if (!all(is.element(conditions, nms))) {
            admisc::stopError(
                "Conditions not found in the data.", ... = ...
            )
        }
    }

    if (!is.element(outcome, nms)) {
        admisc::stopError(
            "Outcome not found in the data.", ... = ...
        )
    }
    
    if (!is.element(type, c("deletion", "corruption"))) {
        admisc::stopError("Unknown perturbation type.", ... = ...)
    }
    if (dependent) {
        if (
            !is.numeric(n.pert) || length(n.pert) != 1L || !is.finite(n.pert) ||
            n.pert != floor(n.pert) || n.pert < 0 || n.pert > nrow(data)
        ) {
            admisc::stopError("Invalid number of case perturbations.", ... = ...)
        }
        if (n.pert == 0) {
            return(1)
        }
    }

    data <- data[, c(conditions, outcome)]
    if (anyNA(data)) {
        admisc::stopError("Missing values are not allowed in retention analyses.", ... = ...)
    }

    # Use one integer code per condition value, including for labelled conditions.
    codes <- lapply(data[, conditions, drop = FALSE], function(x) match(x, unique(x)))
    keys <- do.call(paste, unname(codes))
    configurations <- unique(keys)
    groups <- match(keys, configurations)
    y <- data[[outcome]]
    cpos <- as.double(tabulate(groups[y == 1], length(configurations)))
    cneg <- as.double(tabulate(groups[y == 0], length(configurations)))

    total <- cpos + cneg
    cpos  <- cpos[total >= n.cut]
    cneg  <- cneg[total >= n.cut]
    total <- total[total >= n.cut]
    
    calculatePairs <- function(x, n.pert, type = "deletion") {
        pairsxl <- admisc::combnk(length(total), min(x, length(total)))
        nofsetsxl <- 0
        for (j in seq(ncol(pairsxl))) {
            cposneg <- NULL
            for (l in seq(min(x, nrow(pairsxl)))) {
                cposneg <- c(cposneg, cbind(cpos, cneg)[pairsxl[l, j], ])
            }
            allpairs <- createMatrix(cposneg + 2) - 1
            allpairs <- allpairs[apply(allpairs, 1, function(y) all(y >= 0)), , drop = FALSE]
            linesubset <- rep(TRUE, nrow(allpairs))
            for (l in seq(min(x, nrow(pairsxl)))) {
                linesubset <- linesubset & rowSums(allpairs[, l*2 - c(1, 0)]) >= 1
            }
            allpairs <- allpairs[linesubset & rowSums(allpairs) <= n.pert, , drop = FALSE]
            for (i in seq_len(nrow(allpairs))) {
                lchanges <- rep(FALSE, min(x, nrow(pairsxl)))
                for (l in seq(min(x, nrow(pairsxl)))) {
                    initially <- cpos[pairsxl[l, j]]/total[pairsxl[l, j]]
                    if (type == "deletion") {
                        newtotaless <- total[pairsxl[l, j]] - sum(allpairs[i, l*2 - c(1, 0)])
                        after <- (cpos[pairsxl[l, j]] - allpairs[i, l*2 - 1])/newtotaless
                        lchanges[l] <- ((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                                       ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)
                    }
                    else if (type == "corruption") {
                        after <- (cpos[pairsxl[l, j]] + allpairs[i, l*2] - allpairs[i, l*2 - 1])/total[pairsxl[l, j]]
                        lchanges[l] <- (initially >= incl.cut & after <  incl.cut) |
                                       (initially <  incl.cut & after >= incl.cut)
                    }
                }
                if (all(lchanges)) {
                    combs <- 1
                    for (l in seq(min(x, nrow(pairsxl)))) {
                        combs <- combs*choose(cpos[pairsxl[l, j]], allpairs[i, l*2 - 1])
                        combs <- combs*choose(cneg[pairsxl[l, j]], allpairs[i, l*2])
                    }
                    nofsetsxl <- nofsetsxl + combs*choose(nrow(data) - sum(cposneg), n.pert - sum(allpairs[i, ]))
                }
            }
        }
        return(nofsetsxl)
    }
    
    
    if (dependent) {
        nofsets <- 0
        totalsets <- choose(nrow(data), n.pert)
        for (i in seq_len(min(n.pert, length(total)))) {
            nofsetsxl <- calculatePairs(i, n.pert, type = type)
            if (i %% 2 == 1) {
                nofsets <- nofsets + nofsetsxl
            }
            else {
                nofsets <- nofsets - nofsetsxl
            }
        }
        
        return(as.vector(1 - nofsets/totalsets))
    }
    else {
        pfinal <- 1
        if (type == "deletion") {
            for (l in seq_len(length(total))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq_len(nrow(allpairs))) {
                    newtotaless <- total[l] - allpairs[i, 1] - allpairs[i, 2]
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1])/newtotaless
                    if (((initially >= incl.cut & after <  incl.cut) | newtotaless <  n.cut) |
                        ((initially <  incl.cut & after >= incl.cut) & newtotaless >= n.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        else if (type == "corruption") {
            for (l in seq_len(length(total))) {
                ptmp <- 1
                allpairs <- createMatrix(c(cpos[l], cneg[l]) + 2) - 1
                allpairs <- allpairs[apply(allpairs, 1, function(x) all(x >= 0)), , drop = FALSE]
                allpairs <- allpairs[rowSums(allpairs) >= 1, , drop = FALSE]
                for (i in seq_len(nrow(allpairs))) {
                    initially <- cpos[l]/total[l]
                    after <- (cpos[l] - allpairs[i, 1] + allpairs[i, 2])/total[l]
                    if ((initially >= incl.cut & after < incl.cut) | (initially < incl.cut & after >= incl.cut)) {
                           ptmp <- ptmp - dbinom(allpairs[i, 1], cpos[l], p.pert) * dbinom(allpairs[i, 2], cneg[l], p.pert)
                    }
                }
                pfinal <- pfinal*ptmp
            }
        }
        
        return(pfinal)
    }
}
