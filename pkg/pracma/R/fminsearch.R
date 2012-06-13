##
##  f m i n s e a r c h . R
##


fminsearch <- function(f, x0, ..., minimize = TRUE,
                                   tol = .Machine$double.eps^(2/3)) {
    if (!is.numeric(x0))
        stop("Argument 'x0' must be a numeric vector.")

    scl <- if(minimize) 1 else -1

    fopt <- optim(x0, f, ..., method = "Nelder-Mead",
                  control = list(fnscale = scl, reltol = tol))

    return(list(x = fopt$par, fval = fopt$value))
}
