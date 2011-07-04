##
##  ne v i l l e . R  Neville and Newton Interpolation
##


newtonpoly <- function(x, y, xs = c()) {
    stopifnot(is.numeric(x), is.numeric(y))
    if (length(xs) != 0 && !is.numeric(xs))
        stop("Argument 'xs' must be empty or a numeric vector.")
    x <- c(x); y <- c(y)
    n <- length(x)
    if (length(y) != n)
        stop("Vectors 'x' and 'y' must be of the same length.")

    # Newton's polynomial interpolation formula
    p <- y
    for (k in 2:n) {
        p[k:n] <- (p[k:n] - p[k-1]) / (x[k:n] - x[k-1])
    }
    if (length(xs) == 0) return(p)

    # Evaluating Newton's interpolation formula
    ys <- rep(p[n], length(xs))
    for (k in 1:(n-1)) {
        ys <- p[n-k] + (xs - x[n-k]) * ys
    }
    return(ys)
}


neville <- function(x, y, xs) {
    n <- length(x)
    ys <- y
    for (k in 1:(n-1)) {
        y[1:(n-k)] <- ((xs - x[(k+1):n]) * y[1:(n-k)] +
                       (x[1:(n-k)] - xs) * y[2:(n-k+1)]) /
                       (x[1:(n-k)] - x[(k+1):n])
    }
    ys <- y[1]
    return(ys)
}
