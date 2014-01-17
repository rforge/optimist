##
##  f z e r o . R  Brent-Dekker Algorithm
##


fzero <- function(f, x, ..., maxiter = 100, tol = .Machine$double.eps^(1/2)) {
    if (!is.numeric(x) || length(x) > 2)
        stop("Argument 'x' must be a scalar or a vector of length 2.")

    err <- try(fun <- match.fun(f), silent = TRUE)
    if (class(err) == "try-error") {
        stop("Argument function 'f' not known in parent environment.")
    } else {
        f <- function(x) fun(x, ...)
    }

    if (length(x) == 2) {
        zero <- brent_dekker(f, x[1], x[2], maxiter = maxiter, tol = tol)

    } else {
        if (x != 0) dx <- x/50
        else        dx <- 1/50
        sqrt2 <- sqrt(2)

        a <- b <- x
        fa <- fb <- f(x)
        if (fa == 0) return(list(x = a, fval = fa))

        iter <- 0
        while (fa * fb > 0 && iter < maxiter) {
            iter <- iter + 1
            dx <- sqrt2 * dx
            a  <- a - dx
            fa <- f(a)
            if (fa * fb <= 0) break
            b  <- b + dx
            fb <- f(b)
        }
        if (iter == maxiter) {
            warning("Maximum number of iterations exceeded; no zero found.")
            return(list(x = NA, fval = NA))
        }
        zero <- brent_dekker(f, a, b, maxiter = maxiter, tol = tol)
    }

    return(list(x = zero$root, fval = zero$f.root))
}


brent_dekker <- function(f, a, b,
                        maxiter = 100, tol = .Machine$double.eps^0.5)
# Brent and Dekker's root finding method,
# based on bisection, secant method and quadratic interpolation
{
    stopifnot(is.numeric(a), is.numeric(b),
              length(a) == 1, length(b) == 1)

	x1 <- a; f1 <- f(x1)
	if (f1 == 0) return(list(root = a, f.root = 0, f.calls = 1, estim.prec = 0))
	x2 <- b; f2 <- f(x2)
	if (f2 == 0) return(list(root = b, f.root = 0, f.calls = 1, estim.prec = 0))
	if (f1*f2 > 0.0)
	    stop("Brent-Dekker: Root is not bracketed in [a, b].")

	x3 <- 0.5*(a+b)
	# Beginning of iterative loop
	niter <- 1
	while (niter <= maxiter) {
		f3 <- f(x3)
		if (abs(f3) < tol) {
		    x0 <- x3
		    break
		}

		# Tighten brackets [a, b] on the root
		if (f1*f3 < 0.0) b <- x3 else a <- x3
		if ( (b-a) < tol*max(abs(b), 1.0) ) {
		    x0 <- 0.5*(a + b)
		    break
	    }

		# Try quadratic interpolation
		denom <- (f2 - f1)*(f3 - f1)*(f2 - f3)
		numer <- x3*(f1 - f2)*(f2 - f3 + f1) + f2*x1*(f2 - f3) + f1*x2*(f3 - f1)
		# if denom==0, push x out of bracket to force bisection
		if (denom == 0) {
			dx <- b - a
		} else {
			dx <- f3*numer/denom
		}

		x <- x3 + dx
		# If interpolation goes out of bracket, use bisection.
		if ((b - x)*(x - a) < 0.0) {
			dx <- 0.5*(b - a)
			x  <- a + dx;
		}

		# Let x3 <-- x & choose new x1, x2 so that x1 < x3 < x2.
		if (x1 < x3) {
			x2 <- x3; f2 <- f3
		} else {
			x1 <- x3; f1 <- f3
		}

		niter <- niter + 1
		if (abs(x - x3) < tol) {
		    x0 <- x
		    break
	    }
		x3 <- x;
	}

    if (niter > maxiter)
        warning("Maximum numer of iterations, 'maxiter', has been reached.")

    prec <- min(abs(x1-x3), abs(x2-x3))
    return(list(root = x0, f.root = f(x0),
                f.calls = niter+2, estim.prec = prec))
}
