##
##  q u a d 2 d . R  and  q u a d v .R
##


quad2d <- function(f, xa, xb, ya, yb, n = 32, ...) {
    stopifnot(is.numeric(xa), length(xa) == 1, is.numeric(ya), length(ya) == 1,
              is.numeric(xb), length(xb) == 1, is.numeric(yb), length(yb) == 1)

    fun <- match.fun(f)
    f <- function(x, y) fun(x, y, ...)

    # Get Gauss-Legendre nodes and weights in x- and y-direction.
    cx <- gaussLegendre(n, xa, xb)
    x  <- cx$x
    wx <- cx$w
    cy <- gaussLegendre(n, ya, yb)
    y  <- cy$x
    wy <- cy$w

    # Compute function f at all nodes in x- and y-direction
    mgrid <- meshgrid(x, y)
    Z <- f(mgrid$X, mgrid$Y)

    Q <- wx %*% Z %*% as.matrix(wy)
    return(Q)
}


quadv <- function(f, a, b, n = 32, ...) {
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1)

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    # Get Gauss_legendre nodes and weights
    cc <- gaussLegendre(n, a, b)
    x <- cc$x
    w <- cc$w

    Q <- 0
    for (i in 1:n) Q <- Q + w[i] * f(x[i])

    return(Q)
}
