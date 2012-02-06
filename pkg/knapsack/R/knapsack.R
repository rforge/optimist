##
##  k n a p s a c k . R  Knapsack Problems
##


knapsack <- function(profits, weights, capacity, bounds = NULL,
                     check = TRUE) {
    stopifnot(is.numeric(profits), is.numeric(weights), is.numeric(capacity))
    if (length(capacity) != 1)
        stop("Argument 'capacity' must be a scalar (i.e., of length 1).")

    if (any(weights <= 0))
        stop("'weights' must be a vector of positive numbers.")
    if (any(profits <= 0))
        stop("'profits' must be a vector of positive numbers.")

    n <- length(weights)
    if (n <= 1)
        stop("Length of 'weights' and 'profits' must be greater than 1.")
    m <- n + 3L
    if (length(profits) != n)
        stop("Vectors 'profits' and 'weights' must have the same length.")

    # profits, weights must be integer vextors, all elements < 2^31
    if (any(floor(profits) != ceiling(profits)) ||
        any(floor(weights) != ceiling(weights)) ||
        any(profits >= 2^31) || any(weights >= 2^31))
        stop("All inputs must be positive integers < 2^31 !")

    # profits/weights must be in decreasing order
    o <- order(profits/weights, decreasing = TRUE)
    p <- profits[o]
    w <- weights[o]

    if (is.null(bounds)) {
        R <- knapsack01(p, w, capacity, check = check)
        return(list(indices = o[R$indices],
                    profit = R$profit, capacity = R$capacity))

    } else if (length(bounds) >= 2) {
        stopifnot(is.numeric(bounds))
        if (any(bounds <= 0))
            stop("'bounds' must be a vector of positive numbers.")
        if (length(bounds) != n)
            stop("'profits', 'weights', 'bounds' must have the same length.")
        if (any(floor(bounds) != ceiling(bounds)) ||
            any(bounds >= 2^31))
            stop("All inputs must be positive integers < 2^31 !")

        if (any(bounds > capacity/weights)) {
            warning("Some bounds exceed the capacity -- they will be reduced.")
            bounds <- pmin(bounds, floor(capacity/weights))
        }

        b <- bounds[o]
        R <- knapsackbd(p, w, b, capacity, check = check)
        return(list(indices = o[R$indices], nitems = R$nitems,
                    profit = R$profit, capacity = R$capacity))

        } else if (length(bounds) == 1) {
            if (is.infinite(bounds)) {
                R <- knapsackub(p, w, capacity, check = check)
                return(list(indices = o[R$indices], nitems = R$nitems,
                            profit = R$profit, capacity = R$capacity))
                
            } else {
                bounds <- rep(bounds, n)
                knapsack(profits, weights, capacity, bounds, check = check)
            }

    } else {
        stop("Unknown 'bounds' type. Inform the maintainer.")
    }
}

#-- --------------------------------------------------------------------
knapsack01 <- function(profits, weights, capacity, check = TRUE)
{
    n  <- length(profits)
    m  <- n + 1
    xs <- numeric(m)            # solution vector of 0/1
    zz <- as.integer(0)         # value of optimal solution
    ck <- 1                     # ck <- if (check) 1 else 0

    # dummy parameters, used within Fortran code
    xx <- mn <- ps <- ws <- zs <- vector("integer", m)

    S <- .Fortran("mt1", as.integer(n), as.integer(profits),
                         as.integer(weights), as.integer(capacity),
                         zz = as.integer(zz), xs = as.integer(xs),
                         as.integer(m), as.integer(ck),
                         as.integer(xx), as.integer(mn),
                         as.integer(ps), as.integer(ws), as.integer(zs),
                         PACKAGE = 'knapsack')

    zz <- S$zz
    if (zz > 0) {
        inds <- which(S$xs == 1)
        capa <- sum(weights[inds])
        prof <- sum(profits[inds])  # S$zz
        return(list(indices = inds, profit = prof, capacity = capa))

    } else {
        if (zz == -1)      stop("Not satisfied: '2 <= n <= m-1'.")
        else if (zz == -2) stop("Not satisfied: 'All input must be positive integers'.")
        else if (zz == -3) stop("Not satisfied: 'Maximum weight smaller than capacity'.")
        else if (zz == -4) stop("Not satisfied: 'Sum of weights greater than capacity'.")
        else if (zz == -5) stop("Not satisfied: 'Profits/weights is not decreasing'.")
        else               stop("Unknown error: Inform the package maintainer.")
    }
}

#-- --------------------------------------------------------------------
knapsackbd <- function(profits, weights, bounds, capacity, check = TRUE)
{
    n <- length(profits)
    m <- n + ceiling(sum(log2(bounds))) + 3

    xs <- numeric(n)            # solution vector of 0/1
    zz <- as.integer(0)         # value of optimal solution
    jck <- 1L                   # ck <- if (check) 1 else 0
    jfo <- 1L                   # exact solution required
    jfs <- 0L                   # items not sorted
    jub <- 0L                   # upper bound, needed if jfo == 0

    # dummy parameters, used within Fortran code
    id1 <- id2 <- id3 <- id4 <- id5 <- id6 <- id7 <- vector("integer", m)
    rd8 <- numeric(m)

    S <- .Fortran("mtb2", as.integer(n), as.integer(profits),
                          as.integer(weights), as.integer(bounds),
                          as.integer(capacity),
                          zz = as.integer(zz), xs = as.integer(xs),
                          as.integer(m), as.integer(m),
                          as.integer(jfo), as.integer(jfs),
                          as.integer(jck), as.integer(jub),
                          id1, id2, id3, id4, id5, id6, id7, rd8,
                          PACKAGE = 'knapsack')

    zz <- S$zz
    xs <- S$xs
    if (zz > 0) {
        inds <- (1:n)[xs > 0]
        nits <- xs[xs > 0]
        capa <- sum(xs * weights)
        prof <- sum(xs * profits)
        return(list(indices = inds, nitems = nits, profit = prof, capacity = capa))

    } else {
        if (zz == -1)      stop("Not satisfied: '2 <= n <= m-1'.")
        else if (zz == -2) stop("Not satisfied: 'All input must be positive integers'.")
        else if (zz == -3) stop("Not satisfied: 'Maximum weight smaller than capacity'.")
        else if (zz == -4) stop("Not satisfied: 'Sum of weights greater than capacity'.")
        else if (zz == -5) stop("Not satisfied: 'n + log2(b) <= m-3'.")
        else if (zz == -6) stop("Not satisfied: 'Profits/weights in decreasing order'.")
        else               stop("Unknown error: Inform the maintainer.")
    }
}

#-- --------------------------------------------------------------------
knapsackub <- function(profits, weights, capacity, check = TRUE)
{
    n <- length(profits)
    m <- n + 1L

    ps <- c(profits, 0)
    ws <- c(weights, 0)
    xs <- numeric(m)            # solution vector of 0/1
    zz <- as.integer(0)         # value of optimal solution

    jfo <- 1L                   # exact solution required
    jck <- 1L                   # ck <- if (check) 1 else 0
    jub <- 0L                   # upper bound, needed if jfo == 0

    # dummy parameters, used within Fortran code
    p0 <- w0 <- pp <- vector("integer", m)
    x0 <- rr <- numeric(m)

    S <- .Fortran("mtu2", as.integer(n), as.integer(ps),
                          as.integer(ws), as.integer(capacity),
                          zz = as.integer(zz), xs = as.integer(xs),
                          as.integer(m),
                          as.integer(jfo), as.integer(jck),
                          jub = as.integer(jub),
                          p0, w0, x0, rr, pp,
                          PACKAGE = 'knapsack')

    zz <- S$zz
    xs <- S$xs
    if (zz > 0) {
        inds <- (1:n)[xs > 0]
        nits <- xs[xs > 0]
        capa <- sum(xs[1:n] * weights)
        prof <- sum(xs[1:n] * profits)
        return(list(indices = inds, nitems = nits, profit = prof, capacity = capa))

    } else {
        if (zz == -1)      stop("Not satisfied: '2 <= n <= m-1'.")
        else if (zz == -2) stop("Not satisfied: 'All input must be positive integers'.")
        else if (zz == -3) stop("Not satisfied: 'Maximum weight smaller than capacity'.")
        else if (zz == -4) stop("Not satisfied: 'Sum of weights greater than capacity'.")
        else if (zz == -5) stop("Not satisfied: 'n + log2(b) <= m-3'.")
        else if (zz == -6) stop("Not satisfied: 'Profits/weights in decreasing order'.")
        else               stop("Unknown error: Inform the maintainer.")
    }
}
#-- --------------------------------------------------------------------

