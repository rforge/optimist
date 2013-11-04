\name{Beta}
\alias{sp.beta}
\title{
  Beta Function
}
\description{
  The Beta function is defined as 
  \deqn{B(p, q) = \int_0^1 x^{p-1} (1 - x)^{q-1} dx , \qquad Re(p), Re(q) > 0}
  where \eqn{p} and \eqn{q} are parameters.
}
\usage{
sp.beta(p, q)
}
\arguments{
  \item{p, q}{real arguments greater 0.}
}
\details{
  The Beta function is computed through its relation to the Gamma function:
  \deqn{B(p, q) = \frac{\Gamma(p) \Gamma(q)}{\Gamma(p + q)}}
}
\value{
  Returns the result of the Beta function.
}
\examples{
sp.beta(1, 1)
}
