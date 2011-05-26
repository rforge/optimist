deval <- function(x, y, xp, idx = NULL) {
    stopifnot(is.vector(x, mode = "numeric"), is.numeric(y),
              is.vector(xp, mode= "numeric"))
    if (is.vector(y)) y <- as.matrix(y)
    if (length(x) != nrow(y))
        stop("Length of 'x' must be equal to the number of rows in 'y'.")
    if (is.null(idx)) idx <- 1:ncol(y)
    if (! all(idx %in% 1:ncol(y)))
        stop("Indices 'idx' must be between 1 and no. of columns of 'y'.")

    fint <- findInterval(xp, x)
    flen <- length(fint)

    yp <- matrix(NA, nrow = flen, ncol = length(idx))

    for (i in 1:flen) {
        fi <- fint[i]
        if (fi == 0 || fi == length(x)) next

        yp[i, ] <- y[fi, idx] + 
                   (xp[i] - x[fi])/(x[fi+1] - x[fi]) * (y[fi+1, idx] - y[fi, idx])
    }

    if (flen == 1) yp <- drop(yp)
    return(yp)
}


deeve <- function(x, y, yv = 0, idx = NULL){
    stopifnot(is.vector(x, mode = "numeric"), is.numeric(y),
              is.numeric(yv), length(yv) == 1)
    if (is.vector(y)) y <- as.matrix(y)
    if (length(x) != nrow(y))
        stop("Length of 'x' must be equal to the number of rows in 'y'.")
    if (is.null(idx)) idx <- ncol(y)
    else if (length(idx) > 1) {
        idx <- idx[1]
        warning("Several indices found; only accepting the first one.")
    }

    y <- y[, idx]
    if (yv < min(y) || yv > max(y))
        return(NA)

    # findInterval() needs nondecreasingly sorted vector
    i1 <- which(yv == y)
    y1 <- (y - yv)[1:(length(y)-1)]
    y2 <- (y - yv)[2:length(y)]
    i2 <- which(y1 * y2 < 0)
    fint <- sort(c(i1, i2))
    flen <- length(fint)

    x0 <- numeric(flen)
    for (i in 1:flen) {
        fi <- fint[i]
        x0[i] <- (yv - y[fi]) / (y[fi+1] - y[fi]) * (x[fi+1] - x[fi]) + x[fi]
    }

    return(x0)
}
