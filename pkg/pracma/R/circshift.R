circshift <- function(a, sz) {
    if (is.null(a)) return(a)
    
    if (is.vector(a) && length(sz) == 1) {
        n <- length(a)
        s <- sz %% n
        if (s > 0) a <- a[c((s+1):n, 1:s)]

    } else if (is.matrix(a) && length(sz) == 2) {
        n <- nrow(a); m <- ncol(a)
	    s1 <- sz[1] %% n
	    s2 <- sz[2] %% m
        i1 <- if (s1 > 0) c((s1+1):n, 1:s1) else c(1:n)
        i2 <- if (s2 > 0) c((s2+1):m, 1:s2) else c(1:m)
	    a <- a[i1, i2]

    } else
        stop("Length of 'sz' must be equal to the no. of dimensions of 'a'.")

	return(a)
}
