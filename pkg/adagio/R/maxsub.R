##
##  m a x s u b . R  Maximal Sum Subaarray
##


maxsub <- function(x, inds = TRUE, compiled = TRUE) {
    if (!is.numeric(x))
        stop("Argument 'x' must be a numeric vector.")
    n <- length(x)
    i1 <- 0; i2 <- 0
    s <- 0.0

    if (compiled) {
        S <- .Fortran("maxsubf", x = as.numeric(x), n = as.integer(n),
                                 s = as.numeric(s),
                                 i1 = as.integer(i1), i2 = as.integer(i2),
                                 PACKAGE = "adagio")
        if (inds)
            return(list(sum = S$s, inds = c(S$i1, S$i2)))
        else
            return(S$s)

    } else {
        if (!inds) {
            m1 <- m2 <- 0.0
            for (i in 1:n) {
                m2 <- max(m2 + x[i], 0.0)
                m1 <- max(m1, m2)
            }
            return(m1)
        } else {
            m1 <- m2 <- 0
            p1 <- p2 <- 0
            q1 <- q2 <- 1

            for (i in 1:n) {
                if (m2 > -x[i]) {
                    m2 <- m2 + x[i]
                    q2 <- i
                    if (m2 > m1) {
                        m1 <- m2
                        p1 <- q1; p2 <- q2
                    }
                } else {
                    m2 <- 0
                    q1 <- q2 <- i+1
                }
            }
            return(list(sum = m1, inds = c(p1, p2)))
        }
    }
}
