\name{sbplx}
\alias{sbplx}
\title{
  Subplex Algorithm
}
\description{
  Subplex is a variant of Nelder-Mead that uses Nelder-Mead on a sequence of
  subspaces.
}
\usage{
sbplx(x0, fn, lower = NULL, upper = NULL,
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
  SUBPLEX is claimed to be much more efficient and robust than the original
  Nelder-Mead, while retaining the latter's facility with discontinuous
  objectives.

  This implementation has explicit support for bound constraints (via the
  method in the Box paper as described on the \code{neldermead} help page).
}
\value{
  Prints the same information as the called NLOPT function and returns
  a long list of information about the optimum found.
}
\references{
  T. Rowan, ``Functional Stability Analysis of Numerical Algorithms'', Ph.D.
  thesis, Department of Computer Sciences, University of Texas at Austin, 1990.
}
\note{
  It is the request of Tom Rowan that reimplementations of his algorithm shall
  not use the name `subplex'.
}
\seealso{
  \code{subplex::subplex}
}
\examples{
# Will be filled later
}
