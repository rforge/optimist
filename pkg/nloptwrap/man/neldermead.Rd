\name{neldermead}
\alias{neldermead}
\title{
  Nelder-Mead Simplex
}
\description{
  An implementation of almost the original Nelder-Mead simplex algorithm.
}
\usage{
neldermead(x0, fn, lower = NULL, upper = NULL,
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
  Provides xplicit support for bound constraints, using essentially the method
  proposed in [Box].
  Whenever a new point would lie outside the bound constraints the point is
  moved back exactly onto the constraint.
}
\value{
  Prints the same information as the called NLOPT function and returns
  a long list of information about the optimum found.
}
\references{
  J. A. Nelder and R. Mead, ``A simplex method for function minimization,''
  The Computer Journal 7, p. 308-313 (1965).

  M. J. Box, ``A new method of constrained optimization and a comparison with
  other methods,'' Computer J. 8 (1), 42-52 (1965).
}
\note{
  The author of NLopt would tend to recommend the Subplex method instead.
}
\seealso{
  \code{dfoptim::nmk}
}
\examples{
# Will be filled later
}
