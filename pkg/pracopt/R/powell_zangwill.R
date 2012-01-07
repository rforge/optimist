powell_zangwill <- function(x0, f, g = NULL, kmax = 100,
                   tol1 = 1e-4, tol2 = 0.1) {
    if (is.null(g))
        g <- function(x) grad(f, x)

    n <- length(x0)
    xk <- x0
    D <- diag(n)
    dtk <- 1
    al <- numeric(n)
    xkw <- xk
    for (i in 1:n) {
        dkw <- D[, i]
        a0 <- softline(xkw,  dkw, f, g)
        a1 <- softline(xkw, -dkw, f, g)
        if (a0 >= a1) ak <- a0 else ak <- -a1
        al[i] <- ak
        xkw <- xkw + ak*dkw
    }
    m <- which.max(al)
    akm <- al[m]

    dk <- xkw - xk
    lamk <- vnorm(dk)
    a0 <- softline(xkw,  dk, f, g)
    a1 <- softline(xkw, -dk, f, g)
    if (a0 >= a1) ak <- a0 else ak <- -a1
    adk <- ak*dk
    err <- vnorm(adk)

    k <- 1
    while (err >= tol1 && k < kmax) {
        xk <- xkw + adk
        rk <- akm * dtk / lamk
        if (rk > tol2) {
            D[, m] <- dk  #  D <- cbind(D[, 2:n], dk)
            dtk <- rk
        }
        xkw <- xk
        for (i in 1:n) {
            dkw <- D[, i]
            a0 <- softline(xkw,  dkw, f, g)
            a1 <- softline(xkw, -dkw, f, g)
            if (a0 >= a1) ak <- a0 else ak <- -a1
            al[i] <- ak
            xkw <- xkw + ak*dkw
         }
         m <- which.max(al)
         akm <- al[m]

         dk <- xkw - xk
         lamk <- vnorm(dk)
         a0 <- softline(xkw,  dk, f, g)
         a1 <- softline(xkw, -dk, f, g)
         if (a0 >= a1) ak <- a0 else ak <- -a1
         adk <- ak*dk
         err <- vnorm(adk)
         k <- k + 1
    }

    xs <- xkw + adk
    fs <- f(xs)
    
    return(list(xs = xs, fs = fs, kiter = k))
}
