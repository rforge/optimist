###
### ROOTS.R  Matlab ROOTS Function
###

roots <- function(p) {
    if (is.null(p) || length(p) == 0) return(matrix(0, nrow=0, ncol=0))
    if ( !is.vector(p, mode="numeric") &&
         !is.vector(p, mode="complex") )
        stop("Argument p must be a vector of real or complex numbers.")
    if (length(p) == 1) return(matrix(0, nrow=0, ncol=0))

    # Find non-zero entries in p
    inz <- which(p != 0)
    nnz <- length(inz)
    if (nnz == 0) return(c())

    # Strip leading and trailing zeros, but remember the trailing zeros
    q <- p[inz[1]:inz[nnz]]
    r <- rep(0, length(p) - inz[nnz])

    A <- compan(q)
    return(c(r, eig(A)))
}


mroots <- function(p) {
    q  <- polyder(p);      d <- polygcf(p, q)
    pr <- deconv(p, d)$q; rr <- roots(pr)
    r <- c()
    for (i in 1:length(rr)) {
        r <- c(r, rep(rr[i], mroot(p, rr[i])))
    }
    if (length(p)-1 != length(r))
        warning("Not all roots be determined accurately enough.")
    return(r)
}

