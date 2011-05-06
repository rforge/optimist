\name{lgammaz}
\alias{lgammaz}
\title{
  Complex Logarithmic Gamma Function
}
\description{
  Logarithmic Gamma function valid in the entire complex plane.
}
\usage{
lgammaz(z)
}
\arguments{
  \item{z}{Complex number or a complex vector of any size.}
}
\details{
  Computes the natural log Gamma function for complex arguments using the
  Lanczos series approximation.
}
\value{
  Returns a complex vector of function values.
}
\references{
  Zhang, Sh., and J. Jin (1996). Computation of Special Functions.
  Wiley-Interscience, New York.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\note{
  Copyright (c) 2001 Paul Godfrey for a Matlab version available on
  Mathwork's Matlab Central under BSD license.
}
\seealso{
  \code{\link{gamma}}, \code{\link{gammaz}}
}
\examples{
log(gamma(1:10)) - lgammaz(1:10)
lgammaz(-1)
lgammaz(c(-2-2i, -1-1i, 0, 1+1i, 2+2i))
}
\keyword{ math }
