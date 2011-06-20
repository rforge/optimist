\name{quadv}
\alias{quadv}
\title{
  Vectorized Quadrature
}
\description{
  Gaussian quadrature for vector-valued functions.
}
\usage{
quadv(f, a, b, n = 32, ...)
}
\arguments{
  \item{f}{univariate function.}
  \item{a, b}{lower limits of integration; must be finite.}
  \item{n}{number of nodes used.}
  \item{\dots}{additional arguments to be passed to \code{f}.}
}
\details{
  Extends the Gaussian quadrature to vector-valued functions which is
  particularly easy. Function values can even be matrices, but have to
  be of the same size over whole integration interval.

  The function \code{f} needs \emph{not} to be vectorized

  \code{quadv} is not suitable for functions with singularities.
}
\value{
  A numerical value or vector, the computed integral of all the components
  of the function \eqn{f: R \rightarrow R^m}.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\seealso{
\code{\link{quad}}, \code{\link{quad2d}}
}
\examples{
f1 <- function(x) c(x, x^2, x^3)             # row vector
f2 <- function(x) as.matrix(c(x, x^2, x^3))  # column vector
quadv(f1, 0, 2)
quadv(f2, 0, 2)

f3 <- function(x) exp(x[1])                  # not vectorized
quadv(f3, 0, 2)
}
\keyword{ math }
