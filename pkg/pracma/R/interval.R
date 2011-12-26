##
##  i n t e r v a l . R  Interval Functions
##


interv_union <- function(M) {
    if (is.null(M)) return(c())
    if (is.vector(M)) {
        if (length(M) %% 2 == 0) M <- as.matrix(M, ncol = 2)
        else stop("Length of 'M' must be an even number.")
    }
    if (!is.numeric(M) || !is.matrix(M) || ncol(M) != 2)
        stop("Argument 'M' must be a matrix with two columns.")
    if (nrow(M) == 1) return(c(M))
    
    # sorted according to first and second col
    o <- order(M[, 1], M[, 2])
    L <- M[o, 1]; R <- M[o, 2]
    if (any(L > R))
        stop("Left endpoint cannot be larger than right endpoint.")

    k <- 1
    Mnew <- matrix(c(L[k], R[k]), 1, 2)
    for (i in 2:nrow(M)) {
        if (L[i] <= Mnew[k, 2]) {
            Mnew[k, 2] <- max(R[i], Mnew[k, 2])
        } else {
            k <- k+1
            Mnew <- rbind(Mnew, c(L[i], R[i]))
        }
    }
    return(Mnew)
}

interv_intersect <- function(M) {
    if (is.null(M)) return(c())
    if (is.vector(M)) {
        if (length(M) %% 2 == 0) M <- as.matrix(M, ncol = 2)
        else stop("Length of 'M' must be an even number.")
    }
    if (!is.numeric(M) || !is.matrix(M) || ncol(M) != 2)
        stop("Argument 'M' must be a matrix with two columns.")
    if (any(M[, 1] > M[, 2]))
        stop("Left endpoint cannot be larger than right endpoint.")
    if (nrow(M) == 1) return(c(M))

    L <- max(M[, 1]); R <- min(M[, 2])
    return(if (L <= R) c(L, R) else c())
}
