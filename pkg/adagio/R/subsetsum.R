##
##  s u b s e t s u m . R  Subset Sum Algorithm
##


# Assume S decreasing, no elements > t, total sum >= t
subsetsum <- function(S, t) {
    stopifnot(is.numeric(S), is.numeric(t))
    if (length(t) > 1) {
        t <- t[1]
        warning("Length of 't' must be 1; will take only first component.")
    }
    if (any(floor(S) != ceiling(S)) || any(S <= 0) ||
        floor(t) != ceiling(t)      || t <= 0) 
        stop("Arguments 'S' and 't' must be positive integer vectors.")

    if (any(S >= t))
        stop("No element of 'S' shall be greater or equal to 't'.")
    if (sum(S) < t) {
        warning("Total sum of 'S' is smaller than 't'; no solution possible.")
        return(NA)
    }

    n <- length(S)
    inds <- NULL

    L <- c(0)
    for (i in 1:n) {
        L <- unique(c(L, L+S[i]))
        L <- L[L <= t]
        if (max(L) == t) {
            inds <- c(i)
            t <- t - S[i]
            while (t > 0) {
                K <- c(0)
                for (j in 1:n) {
                    K <- unique(c(K, K+S[j]))
                    K <- K[K <= t]
                    if (max(K) == t) break
                }
                inds <- c(inds, j)
                t <- t - S[j]
            }
            break
        }
    }
    return(inds)
}
