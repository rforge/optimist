\name{integral}
\alias{integral}
\title{
  Adaptive Numerical Integration
}
\description{
  Combines several approaches to adaptive numerical integration of functions
  of one variable.
}
\usage{
integral(fun, xmin, xmax,
            method = c("Kronrod","Richardson","Clenshaw","Simpson","Romberg"),
            tol = 1e-8, ...)
}
\arguments{
  \item{fun}{integrand, univariate (vectorized) function.}
  \item{xmin,xmax}{endpoints of the integration interval.}
  \item{method}{integration procedure, see below.}
  \item{tol}{relative tolerance.}
  \item{...}{additional parameters to be passed to the function.}
}
\details{
  Combines the following methods for adaptive numerical integration:
  \itemize{
    \item Gauss-Kronrod
    \item Gauss-Richardson
    \item Clenshaw-Curtis (not yet made adaptive)
    \item adaptive Simpson
    \item Romberg
  }

  If the interval is infinite, \code{quadinf} will be called that accepts
  the same methods as well.
}
\value{
  Returns the integral, no error terms given.
}
\references{
  Davis, Ph. J., and Ph. Rabinowitz (1984). Methods of Numerical Integration.
  Dover Publications, New York.
}
\note{
  \code{integral} does not provide `new' functionality, everything is already
  contained in the functions called here.
}
\seealso{
  \code{\link{quadgk}}, \code{\link{quadgr}},
  \code{\link{simpadpt}}, \code{\link{romberg}}
}
\examples{
##  Very smmoth function
fun <- function(x) 1/(x^4+x^2+0.9)
val <- 1.582232963729353
for (m in c("Kron", "Rich", "Clen", "Simp", "Romb")) {
    Q <- integral(fun, -1, 1, tol = 1e-12, method = m)
    cat(m, Q, abs(Q-val), "\n")}
# Kron 1.582233 3.197442e-13 
# Rich 1.582233 3.197442e-13 
# Clen 1.582233 3.199663e-13 
# Simp 1.582233 3.241851e-13 
# Romb 1.582233 2.555733e-13 

##  Highly oscillating function
fun <- function(x) sin(100*pi*x)/(pi*x)
val <- 0.4989868086930458
for (m in c("Kron", "Rich", "Clen", "Simp", "Romb")) {
    Q <- integral(fun, 0, 1, tol = 1e-12, method = m)
    cat(m, Q, abs(Q-val), "\n")}
# Kron 0.4989868 2.775558e-16 
# Rich 0.4989868 4.440892e-16 
# Clen 0.4989868 2.231548e-14
# Simp 0.4989868 6.328271e-15 
# Romb 0.4989868 1.508793e-13

## Evaluate improper integral
fun <- function(x) log(x)^2 * exp(-x^2)
val <- 1.9475221803007815976
for (m in c("Kron", "Rich", "Clen", "Simp", "Romb")) {
    Q <- integral(fun, 0, Inf, tol = 1e-12, method = m)
    cat(m, Q, abs(Q-val), "\n")}
# Kron 1.947522 1.101341e-13 
# Rich 1.947522 2.928655e-10 
# Clen 1.948016 0.0004935834 
# Simp 1.947522 9.416912e-12 
# Romb 1.952683 0.00516102 
}

\keyword{ math }
