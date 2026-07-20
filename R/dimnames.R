`setColnames` <- function(
    matrix, colnames
) {
    invisible(.Call("C_setColnames", matrix, colnames))
}

`setRownames` <- function(
    matrix, rownames
) {
    invisible(.Call("C_setRownames", matrix, rownames))
}


`setDimnames` <- function(
    matrix, nameslist
) {
    invisible(.Call("C_setDimnames", matrix, nameslist))
}
