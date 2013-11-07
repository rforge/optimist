##
##  Gamma Functions
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
        V <- .Fortran("cgamma",
                        as.numeric(x), as.numeric(y), as.integer(kf),
                        gr = as.numeric(gr), gi = as.numeric(gi),
                        PACKAGE = 'specfun')
        return(V$gr + V$gi * 1i)

    } else {
        stop("Argument 'z' must be a single real or complex number.")
    }
}


sp.gammainc <- function(s, x) {
    if(length(s) != 1 || length(x) != 1 ||
       !is.numeric(s) || !is.numeric(x) ||
       s <= 0 || x < 0)
        stop("Arguments 'a' and 'x' must be single real numbers > 0.")
    if (-x + s*log(x) > 700 || s > 170)
        stop("Argument(s) 's' and/or 'x' are too large.")

    gin <- 0; gim <- 0; gip <- 0
    V <- .Fortran("incog", as.numeric(s), as.numeric(x), gin = as.numeric(gin),
                    gim = as.numeric(gim), gip = as.numeric(gip),
                    PACKAGE = 'specfun')

    return(list(gin = V$gin, gim = V$gim, gip = V$gip))
}
