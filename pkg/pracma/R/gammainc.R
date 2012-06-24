gammainc <- function(x, a) {
    if (!is.numeric(a) || !is.numeric(x))
        stop("All arguments must be real numbers.")
    if (length(a) > 1 || length(x) > 1)
        stop("Arguments must be of length 1; function is not vectorized.")
    if (x == 0 && a == 0) return(1)
    if (a < 0)
        stop("Argument 'a' must be real and nonnegative.")

    if (x > 0)  xam <- -x + a*log(x)
    else        xam <- -x + a*log(x + 0i)
    if (abs(xam) > 700.0 || abs(a) > 170.0) {
        warning("Arguments 'x' and/or 'a' are too large.")
        return(NA)
    }

    # Auxiliary computation of the gamma function
    gamma_ <- function(x, ga) {
        g  <- c(1.0,
                0.5772156649015329e0,-0.6558780715202538,-0.420026350340952e-1,
                0.1665386113822915e0,-0.421977345555443e-1,-0.96219715278770e-2,
                0.72189432466630e-2,-0.11651675918591e-2,-0.2152416741149e-3,
                0.1280502823882e-3,-0.201348547807e-4,-0.12504934821e-5,
                0.11330272320e-5,-0.2056338417e-6,0.61160950e-8,
                0.50020075e-8,-0.11812746e-8,0.1043427e-9,0.77823e-11,
                -0.36968e-11,0.51e-12,-0.206e-13,-0.54e-14,0.14e-14,0.1e-15)

        if (floor(x) == ceiling(x)) {
            if (x > 0.0) {
                ga <- 1.0
                m1 <- x-1
                for (k in 2:m1) ga <- ga * k
            } else {
                ga <- 1.0e+300
            }
        } else {
            if (abs(x) > 1.0) {
                z <- abs(x)
                m <- fix(z)
                r <- 1.0
                for (k in 1:m) {
                    r <- r * (z-k)
                }
                z <- z - m
            } else {
                z <- x
            }
            gr <- g[26]
            for (k in 25:1) {
                gr <- gr * z + g[k]
            }
            ga <- 1.0/(gr*z)
            if (abs(x) > 1.0) {
                ga <- ga*r
                if (x < 0.0) ga <- -pi/(x*ga*sin(pi*x))
            }
        }
        return(list(a=x, ga=ga))
    }

    # Computation of the incomplete gamma function
    gin <- gim <- gip <- 0

    if (x == 0.0) {
        aga <- gamma_(a, ga)
        a <- aga$a; ga <- aga$ga
        gim <- ga
        gip <- 0.0
    } else if (x <= 1.0 + a) {
        s <- 1/a
        r <- s
        for  (k in 1:60) {
            r <- r * x/(a+k);
            s <- s+r;
            if (abs(r/s) < 1e-15) break
        }
        gin <- exp(xam) * s
        aga <- gamma_(a, ga)
        a <- aga$a; ga <- aga$ga
        gip <- gin/ga
        gim <- ga - gin
    } else if (x > 1.0 + a) {
        t0 <- 0
        for  (k in 60:1) {
            t0 <- (k-a)/(1 + k/(x+t0))
        }
        gim <- exp(xam)/(x+t0)
        aga <- gamma_(a, ga)
        a <- aga$a; ga <- aga$ga
        gin <- ga - gim
        gip <- 1 - gim/ga
    }
    return(c(lowinc = Re(gin), uppinc = Re(gim), reginc = Re(gip)))
}
