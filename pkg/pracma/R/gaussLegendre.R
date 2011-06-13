##
##  g a u s s L e g e n d r e . R  Gauss-Legendre Quadrature Formula
##


gaussLegendre <- function(n, a, b) {
    stopifnot(is.numeric(a), length(a) == 1, is.numeric(b), length(b) == 1,
              is.numeric(n), length(n) == 1, n >= 2)

    i <- seq(1, n-1, by = 1)
    d <- i / sqrt(4*i^2 - 1)

    E <- eigen(mdiag(d, 1) + mdiag(d, -1), symmetric = TRUE)
    L <- E$values
    V <- E$vectors

    inds <- order(L)
    x <- L[inds]
    V <- t(V[, inds])
    w <- 2 * V[, 1]^2

    x <- 0.5 * ((b-a)*x + a+b)
    w <- -0.5 * (a-b)*w

    return(list(x = x, w = w))
}
