##
##  n m . R  Nelder-Mead and Subplex
##


neldermead <-
function(x0, fn, lower = NULL, upper = NULL,
                 control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_LN_NELDERMEAD"

    fun <- match.fun(fn)
    fn <- function(x) fun(x, ...)

    S <- nloptr(x0, fn, lb = lower, ub = upper,
                opts = opts)

    print(S)
    return(S)
}


sbplx <-
function(x0, fn, lower = NULL, upper = NULL,
                 control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_LN_SBPLX"

    fun <- match.fun(fn)
    fn <- function(x) fun(x, ...)

    S <- nloptr(x0, fn, lb = lower, ub = upper,
                opts = opts)

    print(S)
    return(S)    
}
