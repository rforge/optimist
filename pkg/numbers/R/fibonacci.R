##
##  f i b o n a c c i . R  Fibonacci Sequence
##


fibonacci <- function(n, sequence = FALSE) {
    if (!is.numeric(n) || length(n) != 1 || floor(n) != ceiling(n) || n < 0)
        stop("Argument 'n' must be a single integer >= 0.")
    if (n <= 1) return(c(1))

    if (sequence) {
        if (n == 2) return(c(1, 2))
        fib <- numeric(n)
        fib[1:2] <- c(1, 2)
        for (k in 3:n) {
            fib[k] <- fib[k-1] + fib[k-2]
        }
    } else {
        if (n <= 1) {
            return(1)
        } else {
            fib = fibonacci(n-1) + fibonacci(n-2)
        }
    }
    return(fib)
}

##  Examples
fibonacci(0)                         # 1
fibonacci(2)                         # 2
fibonacci(2, sequence = TRUE)        # 1 2

# Golden ratio
F <- fibonacci(25, sequence = TRUE)  # ... 75025 121393
f25 <- F[25]/F[24]                   #     1.618033989
phi <- (sqrt(5) + 1)/2
abs(f25 - phi)                       # 7.945178e-11

# Compare recursive with iterative approach
# system.time(F30 <- fibonacci(30))                       # user: 17.075 s
# system.time(F30 <- fibonacci(30, sequence = TRUE)[30])  # user:  0.006 s
