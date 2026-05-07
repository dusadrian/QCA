`nec` <- function(x) {
    !is.na(pmatch(x, "necessity"))
}


`suf` <- function(x) {
    !is.na(pmatch(x, "sufficiency"))
}


`colnms` <- function(mymat, rownms, tilde = FALSE) {
    apply(mymat, 1, function(x) {
        rownms1 <- rownms[x == 1]
        rownms[x == 1] <- if (tilde) paste0("~", rownms1) else tolower(rownms1)
        return(paste(rownms[x > 0], collapse = "*"))
    })
}


`colnms2` <- function(mymat, colnms, tilde = FALSE) {
    chars <- colnms[col(mymat)]
    lowerChars <- if (tilde) paste0("~", chars) else tolower(chars)
    chars <- ifelse(mymat==1L, lowerChars, chars)
    keep <- mymat > 0L
    charList <- split(chars[keep], row(chars)[keep])
    unlist(lapply(charList, paste, collapse = "*"))
}


`splitMainComponents2` <- function(expression) {
    
    expression <- gsub("[[:space:]]", "", expression)
    
    ind.char <- unlist(strsplit(expression, split=""))
    
    
    if (grepl("\\(", expression)) {
        # split the string in individual characters
    
        open.brackets <- which(ind.char == "(")
        closed.brackets <- which(ind.char == ")")
        
        invalid <- ifelse(grepl("\\)", expression), length(open.brackets) != length(closed.brackets), FALSE)
        
        if (invalid) {
            cat("\n")
            stop("Invalid expression, open bracket \"(\" not closed with \")\".\n\n", call. = FALSE)
        }
        
        
        all.brackets <- sort(c(open.brackets, closed.brackets))
        
        if (length(all.brackets) > 2) {
            for (i in seq(3, length(all.brackets))) {
                if (all.brackets[i] - all.brackets[i - 1] == 1) {
                    open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                    closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                }
                
                if (all.brackets[i] - all.brackets[i - 1] == 2) {
                    if (ind.char[all.brackets[i] - 1] != "+") {
                        open.brackets <- setdiff(open.brackets, all.brackets[seq(i - 1, i)])
                        closed.brackets <- setdiff(closed.brackets, all.brackets[seq(i - 1, i)])
                    }
                }
            }
        }
        
        for (i in seq(length(open.brackets))) {
            plus.signs <- which(ind.char == "+")
            last.plus.sign <- plus.signs[plus.signs < open.brackets[i]]
            if (length(last.plus.sign) > 0) {
                open.brackets[i] <- max(last.plus.sign) + 1
            }
            else {
                if (1 == 1) { # ????
                    open.brackets[i] <- 1
                }
            }
            next.plus.sign <- plus.signs[plus.signs > closed.brackets[i]]
            if(length(next.plus.sign) > 0) {
                closed.brackets[i] <- min(next.plus.sign) - 1
            }
            else {
                closed.brackets[i] <- length(ind.char)
            }
        }
                    
        # create an empty list with at least 3 times as many components as number of open brackets (just to make sure I have enough)
        big.list <- vector(mode="list", length = length(open.brackets) + 2)
        
        if (length(open.brackets) == 1) {
            # there is only one open bracket
            if (open.brackets > 1) {
                # there's something before that open bracket
                big.list[[1]] <- paste(ind.char[seq(1, open.brackets - 2)], collapse = "")
            }
            nep <- min(which(unlist(lapply(big.list, is.null))))
            big.list[[nep]] <- paste(ind.char[seq(open.brackets, closed.brackets)], collapse = "")
            if (closed.brackets < length(ind.char)) {
                # there is something beyond the closed bracket
                nep <- min(which(unlist(lapply(big.list, is.null))))
                big.list[[nep]] <- paste(ind.char[seq(closed.brackets + 2, length(ind.char))], collapse = "")
            }
        }
        else {
            for (i in seq(length(open.brackets))) {
                if (i == 1) {
                    # check if there's anything meaningful before the FIRST bracket
                    # i.e. containing a "+" sign, like "A + B(C + D)"
                    # before the first bracket is "A + B", but only B should be multiplied with "C + D"
                    
                    if (open.brackets[1] > 1) {
                        # there is something before the first bracket
                        big.list[[1]] <- paste(ind.char[seq(1, open.brackets[1] - 2)], collapse = "")
                    }
                    
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                }
                else {
                    nep <- min(which(unlist(lapply(big.list, is.null))))
                    big.list[[nep]] <- paste(ind.char[seq(open.brackets[i], closed.brackets[i])], collapse = "")
                    
                    if (i == length(closed.brackets)) {
                        if (closed.brackets[i] < length(ind.char)) {
                            # there is something beyond the last closed bracket
                            nep <- min(which(unlist(lapply(big.list, is.null))))
                    
                            big.list[[nep]] <- paste(ind.char[seq(closed.brackets[i] + 2, length(ind.char))], collapse = "")
                            
                        }
                    }
                    
                }
            }
        }
        
        nulls <- unlist(lapply(big.list, is.null))
        
        if (any(nulls)) {
            big.list <- big.list[-which(nulls)]
        }
        
        
        #### additional, to make a list containing a vector,
        #### rather than separate list components
        big.list <- list(unlist(big.list))
        
    }
    else {
        big.list <- list(expression)
    }
    
    names(big.list) <- expression
    
    return(big.list)
}




#####
# split each main component by separating brackets components
`splitBrackets2` <- function(big.list) {
    big.list <- as.vector(unlist(big.list))
    result <- vector(mode="list", length = length(big.list))
    for (i in seq(length(big.list))) {
        result[[i]] <- trimstr(unlist(strsplit(unlist(strsplit(big.list[i], split="\\(")), split="\\)")), "*")
    }
    names(result) <- big.list
    return(result)
    # return(lapply(big.list, function(x) {
    #     as.list(unlist(strsplit(unlist(strsplit(x, split="\\(")), split="\\)")))
    # }))
}




#####
# split by "+"
`splitPluses2` <- function(big.list) {
    return(lapply(big.list, function(x) {
        x2 <- lapply(x, function(y) {
            plus.split <- unlist(strsplit(y, "\\+"))
            return(plus.split[plus.split != ""])
        })
        names(x2) <- x
        return(x2)
    }))
}


#####
`splitProducts` <- function(x, prod.split) {
    x <- as.vector(unlist(x))
    strsplit(x, split=prod.split)
}



`mvregexp` <- "\\[|\\]|\\{|\\}"




# remove any letters
# > gsub("[a-zA-Z]", "", "A{1} + B{1}*c{0}")
# [1] "{1} + {1}*{0}"

# remove all digits
# > gsub("\\d", "", "A{1} + B{1}*c{0}")
# [1] "A{} + B{}*c{}"

# remove any non-alphanumeric
# > gsub("\\W", "", "A{1} + B{1}*c{0}")
# [1] "A1B1c0"

# remove anything but preserve letters
# > gsub("[^a-zA-Z]", "", "A{1} + B{1}*c{0}")
# [1] "ABc"

# preserve only the brackets (to check if they match):
# > gsub("[^\\{}]", "", "A{1} + B{1}*c{0}")
# [1] "{}{}{}"
