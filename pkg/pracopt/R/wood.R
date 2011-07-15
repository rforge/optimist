##
##  Wood Test Function
##


wood <- function(x) {
    x0 <- x[1];  x1 <- x[2];  x2 <- x[3];  x3 <- x[4]
    s1 <- x1 - x0^2; s2 <- 1 - x0; s3 <- x1 - 1
    t1 <- x3 - x2^2; t2 <- 1 - x2; t3 <- x3 - 1
    t4 <- s3 + t3; t5 <- s3 - t3
    100*s1^2 + s2^2 + 90*t1^2 + t2^2 + 10*t4^2 + t5^2/10
}


##  Example
x0 <- c(-3, -1, -3, -1)
wood(x0)

