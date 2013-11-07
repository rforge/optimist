\name{Gamma}
\alias{sp.gamma}
\alias{sp.lgamma}
\alias{sp.gammainc}
\title{
  Gamma Functions
}
\description{
  The Gamma function, denoted by \eqn{\Gamma(z)}, is defined by
  \deqn{\Gamma(z) = \int_0^{\infty} t^{z-1} e^{-t} dt ,\qquad Re(z) > 0}
  and can be extended by analytic continuation for \eqn{Re(z) \le 0}.

  The incomplete lower and upper Gamma functions are defined as
  \deqn{\gamma(s, x) = \int_0^x t^{s-1} e^{-t} dt, \qquad 
        \Gamma(s, x) = \int_x^{\infty} t^{s-1} e^{-t} dt}
  for \eqn{Re(s) > 0} such that \eqn{\Gamma(s) = \gamma(s, x) + \Gamma(s, x)}, 
  and the normalized (lower) incomplete Gamma function is 
  \eqn{P(s, x) = \frac{\gamma(s)}{\Gamma(s)}.}
}
\usage{
sp.gamma(z)
sp.lgamma(z)

sp.gammainc(s, x)
}
\arguments{
  \item{z}{real or complex argument.}
  \item{s,x}{real arguments greater zero.}
}
\details{
  The Gamma function can be calculated by its series expansion, 
  or by an asymptotic expansion, known as Stirling's formula.

  \code{sp.gamma} computes the Gamma function combining both approaches, and
  \code{sp.lgamma} computes the (natural) logarithm of the absolute value of
  \eqn{\Gamma(z)}.

  \code{sp.gammainc} computes the lower and upper incomplete Gamma function 
  as well as the normalized lower incomplete Gamma function. These functions 
  are here only computed for real numbers \code{s, z > 0}.
}
\value{
  Returns the result of the Gamma function (or its natural logarithm), while 
  \code{sp.gammainc} returns a list with \code{gin} the lower, \code{gim} the 
  upper and \code{gip} the normalized lower incomplete gamma function value.
}
\examples{
    sp.gamma(1/4); sp.gamma(-1/4)   #  3.625609908221908  -4.90166680986071
    sp.gamma(1/3); sp.gamma(-1/3)   #  2.678938534707748  -4.06235381827920
    sp.gamma(1/2); sp.gamma(-1/2)   #  1.772453850905516  -3.54490770181103
    sp.gamma(2/3); sp.gamma(-2/3)   #  1.354117939426400  -4.01840780206162
    sp.gamma(3/4); sp.gamma(-3/4)   #  1.225416702465178  -4.83414654429588
}
