`combint` <- function (
    n, k, ogte = 0, zerobased = FALSE
) {
    n <- as.integer(n)
    k <- as.integer(k)
    ogte <- as.integer(ogte)
    zerobased <- as.integer(zerobased)
    
    # the order n, k, ogte and zerobased is important for C
    .Call("C_ombnk", list(n = n, k = k, ogte = ogte, zerobased = zerobased), PACKAGE = "QCA")
}
