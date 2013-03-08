##
##  i n t e g r a l . R  Numerical Integration
##


integral <- function(fun, xmin, xmax,
                method = c("Kronrod","Richardson","Clenshaw","Simpson","Romberg"),
                tol = 1e-8, ...)
{
    stopifnot(is.numeric(xmin), length(xmin) == 1,
              is.numeric(xmax), length(xmax) == 1)
    fun <- match.fun(fun)
    f <- function(x) fun(x, ...)
    g <- function(x) (1/x^2) * f(1/x)

    if (xmin == xmax) return(0)

    method <- match.arg(method)
    if (is.finite(xmin) && is.finite(xmax)) {
        Q <- switch(method,
                "Clenshaw"   = clenshaw_curtis(f, xmin, xmax),
                "Kronrod"    = quadgk(f, xmin, xmax, tol = tol),
                "Richardson" = quadgr(f, xmin, xmax, tol = tol)$value,
                "Romberg"    = romberg(f, xmin, xmax, tol = tol)$value,
                "Simpson"    = simpadpt(f, xmin, xmax, tol = tol)
                )
    } else {
        Q <- quadinf(f, xmin, xmax, tol = tol, method = method)
    }

    return(Q)
}
