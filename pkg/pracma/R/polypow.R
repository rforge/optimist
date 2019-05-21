##
##  p ol y p o w . R  Polynomial Powers
##


polypow <- function(p, n){
    if ( !is.vector(p, mode="numeric") && !is.vector(p, mode="complex") )
        stop("Arguments 'p' must be a real or complex vector.")
    if ( !is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 0 )
        stop("Argument 'n' must be a non-negative integer.")

    pp <- c(1)
    while (n > 0) {
    	pp <- polymul(pp, p)
    	n <- n - 1
    }

    return(pp)
}


polygcf <- function(p, q, tol=1e-12) {
    if ( !is.vector(p, mode="numeric") && !is.vector(p, mode="complex") )
        stop("Arguments 'p' must be a real or complex vector.")
    if ( !is.vector(q, mode="numeric") && !is.vector(q, mode="complex") )
        stop("Arguments 'q' must be a real or complex vector.")

    np <- Norm(p)
    pd <- polydiv(p,q)
    a <- pd$d; r0 <- pd$r
    if (Norm(r0) > np*tol) {
        pd <- polydiv(q,r0)
        a <- pd$d; r1 <- pd$r
        if (Norm(r1) > np*tol) {
            rn <- 1
            while (Norm(rn) > np*tol) {
                pd <- polydiv(r0,r1)
                a <- pd$d; rn <- pd$r
                r0 <- r1
                r1 <- rn
            }
            g <- r0
        } else {
            g <- r0
        }
    } else {
       g <- q
    }
    # g <- g / g[1]
    return(g)
}


polytrans <- function(p, q){
    if ( (!is.vector(p, mode="numeric") && !is.vector(p, mode="complex")) ||
         (!is.vector(q, mode="numeric") && !is.vector(q, mode="complex")) )
        stop("Arguments 'p' and 'q' must be real or complex vectors.")

    n <- length(p)
    if (length(p) == 1)
        return(p)

    pt <- 0
    for (i in 1:n) {
    	pt <- polyadd(pt, p[i]*polypow(q, n-i))
    }

    return(pt)
}
