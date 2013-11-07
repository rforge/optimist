##
##  Beta Functions
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


sp.betainc <- function(x, p, q) {
    if (length(x) != 1 || length(p) != 1 || length(q) != 1 || 
        !is.numeric(x) || !is.numeric(p) || !is.numeric(q))
        stop("All arguments 'x', 'p', 'q' must be single real numbers.")
    if (p <= 0 || q <= 0)
        stop("Arguments 'p' and 'q' must be greater than zero.")
    if (x < 0 || x > 1)
        stop("Argument 'x' must be a real number between 0 and 1.")

    bix <- 0
    V <- .Fortran("incob", as.numeric(p), as.numeric(q),
                    as.numeric(x), bix = as.numeric(bix),
                    PACKAGE = 'specfun')
    return(V$bix)
}