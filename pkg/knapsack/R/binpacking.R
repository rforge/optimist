##
##  b i n p a c k i n g . R  Bin Packing Problem
##


binpacking <- function(weights, cap, back = -1, jck = 0, lb = 0) {
    stopifnot(is.numeric(weights), is.numeric(cap))
    if (length(cap) != 1)
        stop("Argument 'cap' must be a scalar (i.e., of length 1).")
    if (floor(cap) != ceiling(cap) || cap >= 2^31)
        stop("'cap' must be a positive imteger smaller than 2^31.")

    if (any(weights <= 0))
        stop("'weights' must be a vector of positive numbers.")
    if (any(floor(weights) != ceiling(weights)) || any(weights >= 2^31))
        stop("All 'weights' must be positive integers < 2^31 !")

    if (is.unsorted(rev(weights))) {
        stop("'weights' must be a vector of decreasing numbers.")
    }
    if (any(weights > cap)) {
        stop("All 'weights' must be smaller than 'cap'.")
    }

    n <- length(weights)
    if (n == 1) {
        return(list(nbins = 1, xbins = c(1)))
    }

    z <- -1                     # solution found: 0
    xstar <- numeric(n)         # xstar[j]: which bin item j belongs to
    jdim <- as.integer(n)       # dimension of all dummy arrays
    # back <- as.integer(back)  # exact solution: -1, else no. of backtracks 
    # jck <- as.integer(jck)    # no input check:  0, else 1
    # lb <- as.integer(lb)      # lower bound, lb >= sum(weights)/cap
    alb <- ceiling(sum(weights)/cap)
    if (lb < alb) lb <- alb     # lower bound on solution value

    # dummy arrays
    wr <- xstarr <- dum <- res <- rel <- x <- r <- wa <- vector("integer", jdim)
    wb <- kfix <- fixit <- xred <- ls <- lsb <- xheu <- vector("integer", jdim)

    R <- .Fortran("mtp",
                  as.integer(n), as.integer(weights), as.integer(cap),
                  z = as.integer(z), xstar = as.integer(xstar),
                  as.integer(jdim), as.integer(back),
                  as.integer(jck), as.integer(lb),
                  wr, xstarr, dum, res, rel, x, r, wa,
                  wb, kfix, fixit, xred, ls, lsb, xheu)

    if (R$z <= 0) {
        warning("Bin packing routine failed to find solution.", call. = FALSE)
        return(list(nbins = R$z, xbins = NA))
    } else {
        return(list(nbins = R$z, xbins = R$xstar))
    }
}