\name{Gamma}
\alias{sp.gamma}
\alias{sp.lgamma}
\title{
  Gamma and Logarithmic Gamma Function
}
\description{
  The Gamma function, denoted by \eqn{\Gamma(z)}, is defined by
  \deqn{\Gamma(z) = \int_0^{\infty} t^{z-1} e^{-t} dt ,\qquad Re(z) > 0}
  and can be extended by analytic continuation for \eqn{Re(z) \le 0}.
}
\usage{
sp.gamma(z)
sp.lgamma(z)
}
\arguments{
  \item{z}{real or complex argument.}
}
\details{
  The Gamma function can be calculated by its series expansion, 
  or by an asymptotic expansion, a known as Stirling's formulas.

  \code{sp.gamma} computes the Gamma function combining both approaches, and
  \code{sp.lgamma} computes the (natural) logarithm of the absolute value of
  \eqn{\Gamma(z)}.
}
\value{
  Returns the result of the Gamma function (or its natural logarithm).
}
\examples{
    sp.gamma(1/4); sp.gamma(-1/4)   #  3.625609908221908  -4.90166680986071
    sp.gamma(1/3); sp.gamma(-1/3)   #  2.678938534707748  -4.06235381827920
    sp.gamma(1/2); sp.gamma(-1/2)   #  1.772453850905516  -3.54490770181103
    sp.gamma(2/3); sp.gamma(-2/3)   #  1.354117939426400  -4.01840780206162
    sp.gamma(3/4); sp.gamma(-3/4)   #  1.225416702465178  -4.83414654429588
}
