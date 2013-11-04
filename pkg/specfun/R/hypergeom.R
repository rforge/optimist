##
##  Hypergeometric Function
##


sp.hypergeom <- function(a, b, c, z) {
    stopifnot(is.numeric(a), is.numeric(b), is.numeric(c))

    if (abs(z) > 1)
        stop("Argument 'z' must fulfill abs(z) <= 1.")
    if (length(z) != 1)
        stop("Argument 'z' must be a single real or complex number.")

    if (is.numeric(z)) {
        if (floor(z) == ceiling(z) && z < 0) return(NaN)
        hf <- 0
        V <- .Fortran("hygfx", as.numeric(a), as.numeric(b), as.numeric(c),
                        as.numeric(z), hf = as.numeric(hf),
                        PACKAGE = 'specfun')
        return(V$hf)

    } else if (is.complex(z)) {
        zhf <- 0 + 0i
        V <- .Fortran("hygfz", as.numeric(a), as.numeric(b), as.numeric(c),
                        as.complex(z), zhf = as.complex(zhf),
                        PACKAGE = 'specfun')
        return(V$zhf)

    } else {
        stop("Argument 'z' must be a single real or complex number.")
    }
}
