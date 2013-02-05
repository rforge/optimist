##
##  a u g l a g 2 . R  Augmented Lagrangian
##


auglag2 <-
function(x0, fn, gr = NULL, lower = NULL, upper = NULL,
            hin = NULL, hinjac = NULL, heq = NULL, heqjac = NULL,
            nl.info = FALSE, control = list(), ...)
{
    local_opts <- list(algorithm = "NLOPT_LD_LBFGS",
                       xtol_rel  = 1e-4)
    opts <- nl.opts(control)
    opts["algorithm"] <- "NLOPT_AUGLAG"
    opts[["local_opts"]] <- local_opts

    fun <- match.fun(fn)
    fn  <- function(x) fun(x, ...)

    if (is.null(gr)) {
        gr <- function(x) nl.grad(x, fn)
    }

    if (!is.null(hin)) {
        .hin <- match.fun(hin)
        hin <- function(x) (-1) * .hin(x)   # change  hin >= 0  to  hin <= 0 !
        if (is.null(hinjac)) {
            hinjac <- function(x) nl.jacobian(x, hin)
        } else {
            .hinjac <- match.fun(hinjac)
            hinjac <- function(x) (-1) * .hinjac(x)
        }
    }

    S0 <- nloptr(x0,
                eval_f = fn,
                eval_grad_f = gr,
                lb = lower,
                ub = upper,
                eval_g_ineq = hin,
                eval_jac_g_ineq = hinjac,
                eval_g_eq = heq,
                eval_jac_g_eq = heqjac,
                opts = opts)

    if (nl.info) print(S0)
    S1 <- list(par = S0$solution, value = S0$objective, iter = S0$iterations,
                convergence = S0$status, message = S0$message)
    return(S1)
}
