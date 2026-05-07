`replaceCategories` <- function(x, categories = NULL) {

    mv <- any(grepl("\\[|\\{", x)) && all(grepl("\\]|\\}", x))
    target <- replacement <- c()
    nms <- names(categories)
    
    for (i in seq(length(categories))) {
        if (mv) {
            values <- seq(length(categories[[i]])) - 1
            target <- c(target, paste0(nms[i], "[", values, "]"))
        }
        else {
            target <- c(target, paste0("~", nms[i]), nms[i])
        }

        replacement <- c(replacement, categories[[i]])
    }

    names(replacement) <- target

    for (i in seq(length(x))) {
        x[i] <- admisc::replaceText(
            x[i],
            replacement = replacement,
            target = target
        )
    }

    return(x)
}
