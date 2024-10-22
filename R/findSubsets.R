`findSubsets` <- function(
    input, noflevels = NULL, stop = NULL, ...
) {
    
    dots <- list(...)
    
    ### backwards compatibility 
        if (is.element("row.no", names(dots)) & missing(input)) {
            input <- dots$row.no
        }
        
        if (is.element("maximum", names(dots))) {
            stop <- dots$maximum
        }
    ### 
    
    stop <- ifelse(missing(stop), prod(noflevels), stop)
    
    result <- lapply(input, function(x) {
        .Call(
            "C_findSubsets",
            x,
            noflevels - 1,
            rev(c(1, cumprod(rev(noflevels))))[-1],
            stop,
            PACKAGE = "QCA"
        )
    })

    return(sort(unique(unlist(result))))
}
