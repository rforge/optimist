\name{newtonpoly}
\alias{newtonpoly}
\title{
  Newtons's Interpolation Method
}
\description{
  Newton's method of polynomial interpolation.
}
\usage{
newtonpoly(x, y, xs = c())
}
\arguments{
  \item{x, y}{x-, y-coordinates of data points defining the polynomial.}
  \item{xs}{either empty, or a vector of points to be interpolated.}
}
\details{
  Straightforward implementation of Newton's method (vectorized in \code{xs}).
}
\value{
  A vector of values at \code{xs} of the polynomial defined by \code{x,y}.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Each textbook on numerical analysis.
}
\seealso{
  \code{\link{neville}}, \code{\link{barylag}}
}
\examples{
p <- poly(c(1, 2, 3))
fp <- function(x) polyval(p, x)

x <- 0:4; y <- fp(x)
xx <- linspace(0, 4, 51)
yy <- newtonpoly(x, y, xx)
\dontrun{
ezplot(fp, 0, 4)
points(xx, yy)}
}
\keyword{ math }
