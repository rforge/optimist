\name{bobyqa}
\alias{bobyqa}
\title{
  Bound Optimization by Quadratic Approximation
}
\description{
  BOBYQA performs derivative-free bound-constrained optimization using an
  iteratively constructed quadratic approximation for the objective function.
}
\usage{
bobyqa(x0, fn, lower = NULL, upper = NULL,
        control = list(), ...)
}
\arguments{
  \item{x0}{starting point for searching the optimum.}
  \item{fn}{objective function that is to be minimized.}
  \item{lower, upper}{lower and upper bound constraints.}
  \item{control}{list of options, see \code{nl.opts} for help.}
  \item{\ldots}{additional arguments passed to the function.}
}
\details{
  This is an algorithm derived from the BOBYQA Fortran subroutine of Powell,
  converted to C and modified for the NLOPT stopping criteria.
}
\value{
  Prints the same information as the called NLOPT function and returns
  a long list of information about the optimum found.
}
\references{
  M. J. D. Powell. ``The BOBYQA algorithm for bound constrained optimization
  without derivatives,'' Department of Applied Mathematics and Theoretical
  Physics, Cambridge England, technical reportNA2009/06 (2009).
}
\note{
  Because BOBYQA constructs a quadratic approximation of the objective,
  it may perform poorly for objective functions that are not
  twice-differentiable.
}
\seealso{
  \code{\link{cobyla}}, \code{\link{newuoa}}
}
\examples{
## Will be provided later.
}
