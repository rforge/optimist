##
##  l b f g s . R  Low-storage BFGS
##


lbfgs <-
function(x0, fn, gr = NULL, lower = NULL, upper = NULL,
            control = list(), ...)
{
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_LD_LBFGS"

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    if (is.null(gr)) {
        gr <- function(x) nl.grad(x, fn)
    }

    S <- nloptr(x0,
                eval_f = fn,
                eval_grad_f = gr,
                lb = lower,
                ub = upper,
                opts = opts)

    print(S)
    return(S)
}
