\name{Hypergeometric}
\alias{sp.hypergeom}
\title{
  Hypergeometric Function
}
\description{
  Gauss' hypergeometric function 2F1 is an infinite series defined as
  \deqn{_2F_1(\alpha, \beta; \gamma; z) = \sum_{k=0}^{\infty} \frac{(\alpha)_k (\beta)_k}{(\gamma)_k} \frac{z^k}{k!}}
  where the Pochhammer symbol \eqn{(\lambda)_k} is employed, that is \eqn{(\lambda)_0 = 1} and
  \deqn{(\lambda)_k = \lambda (\lambda + 1) (\lambda + 2) \ldots (\lambda + k - 1)}

  This series converges for \eqn{|z| < 1} and can be analytically continued 
  into the whole complex plane.
}
\usage{
sp.hypergeom(a, b, c, z) 
}
\arguments{
  \item{a,b,c}{real parameters.}
  \item{z}{real or complex variable.}
}
\details{
  The hypergeometric function is difficult to compute in one program because 
  the four arguments need special care in different situations. Basically, the 
  converging series is used together with recurrence relations and handling of 
  special cases.
}
\value{
  Returns the result of the function evaluation.
}
\examples{
hgf <- sp.hypergeom
## Hypergeometric function sometimes assume rational numbers:
hgf(1/3, 2/3, 5/6, 27/32)               # [1] 1.6 = 8/5
hgf(1/4, 1/2, 3/4, 80/81)               # [1] 1.8 = 9/5

## F(1, 1, 1, z) = F(1, 2, 2, z) = 1 /(1-z)
hgf(1, 1, 1, 1i); hgf(1, 2, 2, 1i)      # [1] 0.5+0.5i


## Elliptic integral of the second kind:
k2 <- 1 - 0.5^2
Ek2 <- pi/2 * hgf(-0.5, 0.5, 1, k2)
print(4 * Ek2, digits = 12)             # ellipse arc length, a = 1, b = 0.5
# [1] 4.84422411027

## An example of very bad computation:
hgf(0.77, 0.88, 0.81, 0.3694-0.8346i)   # true value: 0.6890589-0.6718384i
# [1] -48817.81-6297.21i
}