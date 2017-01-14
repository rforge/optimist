\name{modlog}
\alias{modlog}
\title{
  Modular (or: Discrete) Logarithm
}
\description{
  Realizes the modular (or discrete) logarithm modulo a prime number
  \code{p}, that is determines the unique exponent \eqn{n} such that
  \eqn{g^n = x \mathrm{mod} p}.
}
\usage{
  modlog(g, x, p)
}
\arguments{
  \item{g}{a primitive root mod p.}
  \item{x}{an integer.}
  \item{p}{prime number.}
}
\details{
  The method is in principle a complete search, cut short by "Shank's 
  trick", see Forster (1996, pp. 65ff). \code{g} has to be a primitive
  root modulo \code{p}, otherwise exponentiation is not bijective.
}
\value{
  Returns an integer.
}
\references{
  Forster, O. (1996). Algorithmische Zahlentheorie. Friedr. Vieweg u. 
  Sohn Verlagsgesellschaft mbH, Wiesbaden.
}
\seealso{
  \code{\link{primroot}}
}
\examples{
modlog(11, 998, 1009)  # 505
}
