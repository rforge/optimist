##
##  t e s t f u n c t i o n s . R  Test Functions
##


#-- Rosenbrock's non-convex performance test function ------------------------
fnRosenbrock <- function(x) {
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n-1)]
    sum(100*(x1 - x2^2)^2 + (1 - x2)^2)
}

grRosenbrock <- function(x) {
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- 2*(x[1] - 1) + 400*x[1] * (x[1]^2 - x[2])
    if (n > 2) {
      for (i in 2:(n-1)) {
        g[i] <- -2*(1-x[i]) + 400*x[i]*(x[i]^2-x[i+1]) + 200*(x[i]-x[i-1]^2)
      }
    }
    g[n] <- 200*(x[n] - x[n-1]^2)
    g
}

#-- Nesterov's smooth Chebyshev-Rosenbrock function --------------------------
fnNesterov <- function(x) {
    n <- length(x)
    f <- (1 - x[1])^2 / 4
    for (i in 1:(n-1)) {
        f <- f + (1 + x[i+1] - 2*x[i]^2)^2
    }
    f
}

grNesterov <- function(x) {
    n <- length(x)
    g <- rep(NA, n)
    g[1] <- (x[1] - 1) / 2
    for (i in 1:(n-1)) {
        r = 1 + x[i+1] - 2*x[i]^2
        g[i+1] <- g[i+1] + 2*r
        g[i] <- g[i] - 8*x[i]*r
    }
    g
}

#-- Rastrigin's function -----------------------------------------------------
fnRastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10*cos(2*pi*x))
}

grRastrigin <- function(x) {
    2 * x + 20 * pi * sin(2 * pi * x)
}

#-- Hald's function ----------------------------------------------------------
# xmin = (0.99987763,  0.25358844, -0.74660757,  0.24520150, -0.03749029 )
# fmin = 0.0001223713
initHald <- function(x) {
    # if (length(x) != 5) stop("Parameter vector for Hald must have length 5.")
    i <- 1:21
    t <- -1 + (i - 1)/10
    (x[1] + x[2] * t) / ( 1 + x[3]*t + x[4]*t^2 + x[5]*t^3 ) - exp(t)
}

fnHald <- function(x) {
    f <- initHald(x)
    max(abs(f))
}

grHald <- function(x) {
    g <- rep(NA, 5)
    f <- initHald(x)
    k <- which.max(abs(f))
    s1 <- sign(f[k])
    t <- -1.0 + (k-1) * 0.1
    a <-  1 + x[3]*t + x[4]*t^2 + x[5]*t^3
    b <- x[1] + x[2] * t
    s2 <- sign(t * a * b)
    if (a == 0) return(g)
    g[1:2] <- s1 * c(1, t) / a
    g[3:5] <- s2 * c(t, t^2, t^3) * (x[1] + x[2] * t) / a^2
    g
}

