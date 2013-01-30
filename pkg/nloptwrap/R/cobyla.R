##
##  c o b y l a . R  COBYLA, BOBYQA, and NEWUOA 
##


cobyla <-
function(x0, fn, lower = NULL, upper = NULL,
            hin = NULL, heq = NULL,
            control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_LN_COBYLA"

    f1 <- match.fun(fn)
    fn <- function(x) f1(x, ...)

    if (!is.null(hin)) {
        f2  <- match.fun(hin)
        hin <- function(x) (-1)*f2(x)      # NLOPT expects hin <= 0
    }

    if (!is.null(heq)) {
        f3  <- match.fun(heq)
        heq <- function(x) f3(x)
    }

    S <- nloptr(x0,
                eval_f = fn,
                lb = lower, ub = upper,
                eval_g_ineq = hin,
                eval_g_eq = heq,
                opts = opts)

    print(S)
    return(S)
}


bobyqa <-
function(x0, fn, lower = NULL, upper = NULL,
                 control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithms"] <- "NLOPT_LN_BOBYQA"

    fun <- match.fun(fn)
    fn <- function(x) fun(x, ...)

    S <- nloptr(x0, fn, lb = lower, ub = upper,
                opts = opts)

    print(S)
    return(S)
}


newuoa <-
function(x0, fn, control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithms"] <- "NLOPT_LN_NEWUOA"

    fun <- match.fun(fn)
    fn <- function(x) fun(x, ...)

    S <- nloptr(x0, fn, opts = opts)

    print(S)
    return(S)
}