\name{Beta}
\alias{sp.beta}
\alias{sp.betainc}
\title{
  Beta Functions
}
\description{
  The Beta function is defined as 
  \deqn{B(p, q) = \int_0^1 t^{p-1} (1 - t)^{q-1} dt , \qquad Re(p), Re(q) > 0}
  where \eqn{p} and \eqn{q} are parameters.

  The incomplete Beta function is defined as
  \deqn{B(x, p, q) = \int_0^x t^{p-1} (1 - t)^{q-1} dt}
  for \eqn{Re(p), Re(q) > 0} and \eqn{0 < x < 1}.
}
\usage{
sp.beta(p, q)

sp.betainc(x, p, q)
}
\arguments{
  \item{p, q}{real arguments greater 0.}
  \item{x}{real number between 0 and 1.}
}
\details{
  The Beta function is computed through its relation to the Gamma function:
  \deqn{B(p, q) = \frac{\Gamma(p) \Gamma(q)}{\Gamma(p + q)}}

  The incomplete Beta function satisfies the symmetry relation
  \deqn{B(x, p, q) = 1 - Beta(1 - x, q, p)}.
}
\value{
  Returns the result of the (incomplete) Beta function.
}
\examples{
sp.beta(1, 1)
}
