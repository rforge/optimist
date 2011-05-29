##
##  g r a d i e n t . R  Discrete Derivatives
##


gradient <- function(f, x = rep(1, length(f))) {
    if (length(f) == 0 || length(x) == 0)
        return(c())
    if (length(f) != length(x))
        stop("Length of vectors 'f' and 'x' must be the same.")
    if (is.unsorted(x) || any(diff(x) == 0))
        stop("Argument 'x' must be sorted strongly increasing.")

    n <- length(f)
    if (n == 1) return(NA)

    g <- numeric(n)
    g[1] <- (f[2] - f[1]) / (x[2] - x[1])
    g[n] <- (f[n] - f[n-1]) / (x[n] - x[n-1])

    if (n > 2)
        g[2:(n-1)] <- (f[3:n] - f[1:(n-2)]) / (x[3:n] - x[1:(n-2)])

    return(g)
}