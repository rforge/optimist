##
##  Beta Function
##


sp.psi <- function(z) {
    if (length(z) != 1)
        stop("Argument 'z' must be a single real or complex number.")

    if (is.numeric(z)) {
        if (floor(z) == ceiling(z) && z <= 0) return(NaN)
        ps <- 0
        V <- .Fortran("psi", as.numeric(z), ps = as.numeric(ps),
                        PACKAGE = 'specfun')
        return(V$ps)

    } else if (is.complex(z)) {
        x <- Re(z); y <- Im(z)
        if (y == 0 && floor(x) == ceiling(x) && x <= 0) return(NaN)
        psr <- 0; psi <- 0
        V <- .Fortran("cpsi",
                        as.numeric(x), as.numeric(y),
                        psr = as.numeric(psr), psi = as.numeric(psi),
                        PACKAGE = 'specfun')
        return(V$psr + V$psi * 1i)

    } else {
        stop("Argument 'z' must be a single real or complex number.")
    }
}
