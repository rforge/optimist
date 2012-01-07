\name{Hooke-Jeeves}
\alias{hooke_jeeves}
\title{
  Hooke-Jeeves Method
}
\description{
  An implementation of the Hooke-Jeeves algorithm for derivative-free
  optimization.
}
\usage{
hooke_jeeves(x0, f, h = 1, scale = 1,
             info = FALSE, tol = .Machine$double.eps^(1/2), ...)
}
\arguments{
  \item{x0}{start value.}
  \item{f}{function to be minimized.}
  \item{h}{max. number of iterations.}
  \item{scale}{scale factor, of -1 the maximum will be searched for.}
  \item{info}{logical, whether to print information during the main loop.}
  \item{tol}{relative tolerance, to be used as stopping rule.}
  \item{\ldots}{additional arguments to be passed to the function.}
}
\details{
  This method computes a new point using the values of \code{f} at suitable
  points along the orthogonal coordinate directions around the last point.
}
\value{
  List with following components:
    \item{xmin}{minimum solution found.}
    \item{fmin}{value of \code{f} at minimum.}
    \item{fcalls}{number of function calls.}
    \item{niter}{number of iterations performed.}
}
\references{
  Quarteroni, Sacco, and Saleri (2007), Numerical Mathematics, Springer.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\note{
  For a more elaborate implementation of Hooke-Jeeves see the package
  `dfoptim'.
}
\seealso{
  \code{\link{fletcher_powell}}
}
\examples{
##  Rosenbrock function
hooke_jeeves(c(0, 0), rosenbrock)
# $xmin
# [1] 1 1
# $fmin
# [1] 3.951249e-19
# $niter
# [1] 13
}
\keyword{ optimize }
