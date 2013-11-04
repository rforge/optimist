##
##  Gamma Function
##


sp.gamma <- function(z) {
    if (length(z) != 1)
        stop("Argument 'z' must be a single real or complex number.")

    if (is.numeric(z)) {
        if (floor(z) == ceiling(z) && z <= 0) return(NaN)
        ga <- 0
        V <- .Fortran("gamma", as.numeric(z), ga = as.numeric(ga),
                        PACKAGE = 'specfun')
        return(V$ga)

    } else if (is.complex(z)) {
        x <- Re(z); y <- Im(z)
        if (y == 0 && floor(x) == ceiling(x) && x <= 0) return(NaN)
        gr <- 0; gi <- 0
        kf <- 1L
        V <- .Fortran("cgamma",
                        as.numeric(x), as.numeric(y), as.integer(kf),
                        gr = as.numeric(gr), gi = as.numeric(gi),
                        PACKAGE = 'specfun')
        return(V$gr + V$gi * 1i)

    } else {
        stop("Argument 'z' must be a single real or complex number.")
    }
}


sp.lgamma <- function(z) {
    if (length(z) != 1)
        stop("Argument 'z' must be a single real or complex number.")

    if (is.numeric(z)) {
        if (z <= 0)
            stop("As real number, argument 'z' must be greater than 0.")
        gl <- 0
        kf <- 0L
        V <- .Fortran("lgamma", as.integer(kf),
                        as.numeric(z), gl = as.numeric(gl),
                        PACKAGE = 'specfun')
        return(V$gl)

    } else if (is.complex(z)) {
        x <- Re(z); y <- Im(z)
        if (y == 0 && floor(x) == ceiling(x) && x <= 0) return(NaN)
        gr <- 0; gi <- 0
        kf <- 0L
        S <- .Fortran("cgamma",
                        as.numeric(x), as.numeric(y), as.integer(kf),
                        gr = as.numeric(gr), gi = as.numeric(gi),
                        PACKAGE = 'specfun')
        return(S$gr + S$gi * 1i)

    } else {
        stop("Argument 'z' must be a single real or complex number.")
    }
}
