\name{newuoa}
\alias{newuoa}
\title{
  New Unconstrained Optimization with quadratic Approximation
}
\description{
  NEWUOA solves quadratic subproblems in a spherical trust regionvia a
  truncated conjugate-gradient algorithm. For bound-constrained problems,
  BOBYQA shold be used instead, as Powell developed it as an enhancement
  thereof for bound constraints.
}
\usage{
newuoa(x0, fn, control = list(), ...)
}
\arguments{
  \item{x0}{starting point for searching the optimum.}
  \item{fn}{objective function that is to be minimized.}
  \item{control}{list of options, see \code{nl.opts} for help.}
  \item{\ldots}{additional arguments passed to the function.}
}
\details{
  This is an algorithm derived from the NEWUOA Fortran subroutine of Powell,
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
  NEWUOA may be largely superseded by BOBYQA.
}
\seealso{
  \code{\link{bobyqa}}, \code{\link{cobyla}}
}
\examples{
## Will be provided later.
}
