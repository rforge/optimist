lsqsep <- function(flist, p0, xdata, ydata, const = TRUE) {
    stopifnot(is.numeric(xdata), is.numeric(ydata), is.numeric(p0))
    n <- length(xdata)
    if (length(ydata) != n)
        stop("Numeric arguments 'xdata', 'ydata' must have the same length.")
    m <- length(flist)
    # lapply is.function

    fapply <- function(b) {
        M <- matrix(1, nrow = n, ncol = m + 1)
        for (i in 1:m) {
            fi <- flist[[i]]
            xi <- fi(b[i], xdata)
            M[, i+1] <- xi
        }
        if (!const) M <- M[, 2:ncol(M)]
        a <- qr.solve(M, ydata)         # sum((M %*% a - ydata)^2)
        M %*% a - ydata                 # for lsqnonlin
    }

    # Find the function parameters b
    Lsq <- lsqnonlin(fapply, p0)
    b <- Lsq$x

    # Find the linear parameters a
    M <- matrix(1, nrow = n, ncol = m + 1)
    for (i in 1:m) {
        fi <- flist[[i]]
        xi <- fi(b[i], xdata)
        M[, i+1] <- xi
    }
    if (!const) M <- M[, 2:ncol(M)]
    a <- qr.solve(M, ydata)

    if (const) {
        a0 <- a[1]; a <- a[2:length(a)]
    } else {
        a0 <- 0
    }
    return(list(a0 = a0, a = a, b = b, ssq = NA))
}
