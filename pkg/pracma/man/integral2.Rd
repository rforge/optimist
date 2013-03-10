\name{integral2}
\alias{integral2}
\title{
  Numerically Evaluate Double Integral
}
\description{
  Numerically evaluate a double integral.
}
\usage{
integral2(fun, xmin, xmax, ymin, ymax, sector = FALSE,
            reltol = 1e-6, abstol = 0, maxlist = 5000,
            singular = FALSE, vectorized = TRUE, ...)
}
\arguments{
  \item{fun}{function}
  \item{xmin, xmax}{lower and upper limits of x.}
  \item{ymin, ymax}{lower and upper limits of y.}
  \item{sector}{logical.}
  \item{reltol}{relative tolerance.}
  \item{abstol}{absolute tolerance.}
  \item{maxlist}{maximum length of the list of rectangles.}
  \item{singular}{logical.}
  \item{vectorized}{logical.}
  \item{...}{additional parameters to be passed to the function.}
}
\details{
  Implements the `TwoD' algorithm, that is Gauss-Kronrod with (3, 7)-nodes
  on 2D rectangles.
}
\value{
  Returns a list with \code{Q} the integral and \code{error} the error term.
}
\references{
  Shampine, L. F. (2008). MATLAB Program for Quadrature in 2D.
  Proceedings of Applied Mathematics and Computation, 2008, pp. 266--274.
}
\author{
  Copyright (c) 2008 Lawrence F. Shampine for Matlab code and description
  of the program; adapted and converted to R by Hans W Borchers.
}
\note{
  To avoid recursion, a possibly large matrix will be used and passed between
  subprograms. A more efficient implementation may be possible.
}
\seealso{
  \code{\link{integral}}, \code{cubature:adaptIntegrate}
}
\examples{
fun <- function(x, y) cos(x) * cos(y)
integral2(fun, 0, 1, 0, 1, reltol = 1e-10)
# $Q:     0.708073418273571  # 0.70807341827357119350 = sin(1)^2
# $error: 8.618277e-19       # 1.110223e-16
}
\keyword{ math }
