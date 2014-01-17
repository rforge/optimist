\name{fzero}
\alias{fzero}
\alias{brent_dekker}
\title{
  Root Finding
}
\description{
  Find root of continuous function of one variable.
}
\usage{
fzero(f, x, ..., maxiter = 100, tol = .Machine$double.eps^(1/2))

brent_dekker(f, a, b, maxiter = 100, tol = .Machine$double.eps^0.5)
}
\arguments{
  \item{f}{function whose root is to be found.}
  \item{x}{a point near the root or an interval giving end points.}
  \item{a, b}{left and right end points of an interval.}
  \item{maxiter}{maximum number of iterations.}
  \item{tol}{absolute tolerance.}
  \item{...}{additional arguments to be passed to the function.}
}
\details{
  \code{fzero} tries to find a zero of \code{f} near \code{x}, if \code{x}
  is a scalar. Expands the interval until different signs are found at the
  endpoints or the maximum number of iterations is exceeded.

  If \code{x} is a vector of length two, \code{fzero} assumes \code{x} is
  an interval where the sign of \code{x[1]} differs from the sign of
  \code{x[1]}. An error occurs if this is not true.

  \code{brent_dekker} implements a version of the Brent-Dekker algorithm,
  a well known root finding algorithms for real, univariate, continuous
  functions. The Brent-Dekker approach is a clever combination of secant
  and bisection with quadratic interpolation.

  This approach will not find zeroes of quadratic order.
}
\value{
  \code{fzero} returns a list with
  \item{x}{location of the root.}
  \item{fval}{function value at the root.}
  while \code{brent_dekker} returns a list with
  \item{root}{} 
  \item{f.root}{} 
  \item{f.calls}{} 
  \item{estim.prec}{} 
}
\references{
  Quarteroni, A., R. Sacco, and F. Saleri (2007). Numerical Mathematics.
  Second Edition, Springer-Verlag, Berlin Heidelberg.
}
\note{
  \code{fzero} mimics the Matlab function of the same name, while \code{fzero} 
  in Octave uses a more sophisticated approach with cubic interpolation.
}
\seealso{
  \code{\link{uniroot}}, \code{\link{ridder}}, \code{\link{newtonRaphson}}
}
\examples{
fzero(sin, 3)                    # 3.141593
fzero(cos,c(1, 2))               # 1.570796
fzero(function(x) x^3-2*x-5, 2)  # 2.094551

# Legendre polynomial of degree 5
lp5 <- c(63, 0, -70, 0, 15, 0)/8
f <- function(x) polyval(lp5, x)
brent_dekker(f, 0.6, 1)         # 0.9061798459 correct to 10 places
}
\keyword{ math }
