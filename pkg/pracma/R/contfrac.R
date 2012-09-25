##
##  c o n t f r a c . R  Continuous Fractions
##


rat <- function(x, tol = 1e-6) {
    if (length(x) == 0)
        return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    xs <- c(x)

    n <- length(xs)
    R <- character(n)
    for (i in 1:n) {
        x <- xs[i]
        k <- contfrac(x, tol = tol)$cf
        if (length(k) >= 1) {
            cf <- paste("[ ", k[1], sep="")
        }
        if (length(k) >= 2) {
            cf <- paste(cf, "; ", k[2], sep="")
        }
        if (length(k) >= 3) {
            for (j in 3:length(k)) {
                cf <- paste(cf, ", ", k[j], sep="")
            }
            cf <- paste(cf, "]", sep="")
        }
        R[i] <- cf
    }
    return(R)
}

rats <- function(x, tol = 1e-6) {
    if (length(x) == 0)
        return(c())
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    xs <- c(x)

    n <- length(xs)
    R <- numeric(n)
    for (i in 1:n) {
        x <- xs[i]
        k <- contfrac(x, tol = tol)$rat
        cf <- paste(k[1], "/", k[2], sep="")
        cat(cf, "\n")
        R[i] <- k[1]/k[2]
    }
    invisible(R)    
}
