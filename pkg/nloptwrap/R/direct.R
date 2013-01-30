direct <-
function(fn, lower, upper, control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_GN_DIRECT"

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    x0 <- (lower + upper) / 2

    S <- nloptr(x0,
                eval_f = fn,
                lb = lower,
                ub = upper,
                opts = opts)

    print(S)
    return(S)
}


directL <-
function(fn, lower, upper, control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_GN_DIRECT_L"

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    x0 <- (lower + upper) / 2

    S <- nloptr(x0,
                eval_f = fn,
                lb = lower,
                ub = upper,
                opts = opts)

    print(S)
    return(S)
}
