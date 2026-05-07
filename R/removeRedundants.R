`removeRedundants` <- function(
    implicants, noflevels
) {
    .Call(
        "C_removeRedundants",
        implicants,
        noflevels - 1,
        rev(c(1, cumprod(rev(noflevels))))[-1], # mbase
        PACKAGE = "QCA"
    )
}
