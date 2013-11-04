##
##  Beta Function
##


sp.beta <- function(p, q) {
    if (length(p) != 1 || length(q) != 1 ||
        !is.numeric(p) || !is.numeric(q) || p <= 0 || q <= 0)
        stop("Arguments 'p' and 'q' must be positive real numbers.")

    bt <- 0
    V <- .Fortran("beta", as.numeric(p), as.numeric(q), bt = as.numeric(bt),
                    PACKAGE = 'specfun')
    return(V$bt)
}
