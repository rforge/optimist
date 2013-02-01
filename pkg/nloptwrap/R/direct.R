direct <-
function(fn, lower, upper, nl.info = FALSE, control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_GN_DIRECT"

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    x0 <- (lower + upper) / 2

    S0 <- nloptr(x0,
                eval_f = fn,
                lb = lower,
                ub = upper,
                opts = opts)

    if (nl.info) print(S0)
    S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
                convergence = S0$status, message = S0$message)
    return(S1)
}


directL <-
function(fn, lower, upper, nl.info = FALSE, control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_GN_DIRECT_L"

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    x0 <- (lower + upper) / 2

    S0 <- nloptr(x0,
                eval_f = fn,
                lb = lower,
                ub = upper,
                opts = opts)

    if (nl.info) print(S0)
    S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
                convergence = S0$status, message = S0$message)
    return(S1)
}
