\name{gradient}
\alias{gradient}
\title{
  Gradient (Matlab Style)
}
\description{
  Discrete numerical derivatives.
}
\usage{
gradient(f, x = rep(1, length(f)))
}
\arguments{
  \item{f}{vector of function values}
  \item{x}{x-coordinates of grid points}
}
\details{
  Returns the numerical gradient of the vector or matrix \code{f}.
}
\value{
  Vector of discrete slopes.
}
\note{
  TODO: two-dimensional version.
}
\seealso{
  \code{\link{fderiv}}
}
\examples{
x <- seq(0, 2*pi, length.out = 100)
y <- sin(x)
f <- gradient(y, x)
max(f - cos(x))      #=> 0.00067086

\dontrun{
plot(x, y, type = "l", col = "blue")
lines(x, cos(x), col = "gray", lwd = 3)
lines(x, f, col = "red")
grid()}
}
\keyword{ math }
