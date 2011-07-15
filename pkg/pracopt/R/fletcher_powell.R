##
##  f l e t c h e r _ p o w e l l . R  Davidon-Fletcher-Powell Method
##


fletcher_powell <- function(x0, f, g = NULL,
                      kmax = 100, tol = .Machine$double.eps^(1/2)) {
    eps <- .Machine$double.eps
    if (tol < eps) tol <- eps
    if (!is.numeric(kmax) || length(kmax) > 1 || kmax < 1)
        stop("Argument 'kmax' must be a positive integer.")
    kmax <- floor(kmax)
    if (! is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")
    n <- length(x0)
    if (n == 1)
        stop("Function 'f' is univariate; use some other optimization method.")
    
    # User provided or numerical gradient
    f <- match.fun(f)
    if (! is.null(g)) {
        g <- match.fun(g)
    } else {
        g <- function(x) fpgrad(f, x)   
    }     

    x <- x0  # column vector?
    H <- diag(1, n, n)
    f0 <- f(x)
    g0 <- g(x)  # row vector !

    for (k in 1:kmax) {
        s <- - H %*% as.matrix(g0)  # downhill direction
        # Find minimal f(x + a*s) along the direction s
        f1 <- f0
        z <- sqrt(sum(s^2))
        if (z == 0) return(list(xmin = x, fmin = f1, ninter = k))

        s <- c(s / z)
        a1 <- 0
        a3 <- 1; f3 <- f(x + a3*s)
        while (f3 >= f1) {
            a3 <- a3/2; f3 <- f(x + a3*s)
            if (a3 < tol/2) return(list(xmin = x, fmin = f1, niter = k))
        }

        a2 <- a3/2; f2 <- f(x + a2*s)
        h1 <- (f2 - f1)/a2
        h2 <- (f3 -f2)/(a3 - a2)
        h3 <- (h2 - h1)/a3
        a0 <- 0.5*(a2 - h1/h3); f0 <- f(x + a0*s)
        
        if (f0 < f3) a <- a0
        else         a <- a3

        d <- a * s; dp <- as.matrix(d)
        xnew <- x + d
        fnew <- f(xnew)
        gnew <- g(xnew)
        y <- gnew - g0; yp <- as.matrix(y)
        A <- (dp %*% d) / sum(d * y)
        B <- (H %*% yp) %*% t(H %*% yp) / c(y %*% H %*% yp)
        Hnew <- H + A - B
        if (max(abs(d)) < tol) break

        # Prepare for next iteration
        H <- Hnew
        f0 <- fnew
        g0 <- gnew
        x <- xnew
    }
    return(list(xmin = x, fmin = f(x), niter = k))
}


fpgrad <- function(f, x0, h = 1e-4, ...) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric value.")

    fun <- match.fun(f)
    f <- function(x) fun(x, ...)

    if (length(f(x0)) != 1)
        stop("Function 'f' must be a univariate function of 2 variables.")
    n <- length(x0)

    hh <- rep(0, n)
    gr <- numeric(n)
    for (i in 1:n) {
        hh[i] <- h
        gr[i] <- (-f(x0+2*hh) + 8*f(x0+hh) - 8*f(x0-hh) + f(x0-2*hh)) / (12*h)
        hh[i] <- 0
    }
    return(gr)
}



















