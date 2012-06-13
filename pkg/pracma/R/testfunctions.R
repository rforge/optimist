##
##  t e s t f u n c t i o n s . R  Test Functions
##


##  Rastrigin function
# No. of Vars.:  n >= 1
# Bounds:  -5.12 <= xi <= 5.12
# Local minima:  many
# Minimum:  0.0
# Solution:  xi = 0, i=1:n
rastrigin <- function(x) {
    n <- length(x)
    10*n + sum(x^2 - 10*cos(2*pi*x))
}

rosenbrock <- function(x) {
    n <- length(x)
    x1 <- x[2:n]
    x2 <- x[1:(n-1)]
    sum(100*(x1-x2^2)^2 + (1-x2)^2)
}

wood <- function(x) {
    x0 <- x[1];  x1 <- x[2];  x2 <- x[3];  x3 <- x[4]
    s1 <- x1 - x0^2; s2 <- 1 - x0; s3 <- x1 - 1
    t1 <- x3 - x2^2; t2 <- 1 - x2; t3 <- x3 - 1
    t4 <- s3 + t3; t5 <- s3 - t3
    100*s1^2 + s2^2 + 90*t1^2 + t2^2 + 10*t4^2 + t5^2/10
}
