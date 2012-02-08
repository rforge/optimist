##
##  s u b s e t s u m . R  Subset Sum Problems
##


subsetsum <- function(weights, target, check = TRUE) {
    stopifnot(is.numeric(weights), is.numeric(target))
    if (length(target) != 1)
        stop("Argument 'target' must be a scalar (i.e., of length 1).")
    if (floor(target) != ceiling(target) || target >= 2^31)
        stop("'target' must be a positive imteger smaller than 2^31.")

    if (any(weights <= 0))
        stop("'weights' must be a vector of positive numbers.")
    if (any(floor(weights) != ceiling(weights)) || any(weights >= 2^31))
        stop("All 'weights' must be positive integers < 2^31 !")

    n <- length(weights)
    if (n <= 1)
        stop("Length of 'weights' must be greater than 1.")

    if (sum(weights) <= target) {
        return(list(ssum = sum(weights), inds = 1:n))
    }

    if (target %in% weights) {
        return(list(ssum = target, inds = min(which(weights == target))))
    } else if (any(weights > target)) {
        stop("Argument 'target' must be larger than all 'weights'.")
    }

    z <- 0L
    jdn <- n + 1L
    w <- c(weights, 0)
    x <- vector("integer", jdn)
    jdd <- 5000    # max. length of programmming list; suggested 5000
    itmm <- jdd    # max no. in the core problem + 1; suggested jdd (or 91)
    jck <- 1       # whether to check input data

    # dummy variables
    wo <- ind <- xx <- ws<- zs <- sm <- vector("integer", itmm)
    td1 <- td2 <- td3 <- matrix(0L, nrow = jdd, ncol = 2)

    R <- .Fortran("mtsl",
                  as.integer(n), as.integer(weights), as.integer(target),
                  z = as.integer(z), x = as.integer(x),
                  as.integer(jdn), as.integer(jdd),
                  as.integer(itmm), as.integer(jck),
                  wo, ind, xx, ws, zs, sm, td1, td2, td3)

    inds <- which(R$x == 1)
    return(list(ssum = R$z, inds = inds))
}
