\name{lsqnonlin}
\alias{lsqnonlin}
\title{
  Nonlinear Least-Squares Fitting
}
\description{
  \code{lsqnonlin} solves nonlinear least-squares problems, including
  nonlinear data-fitting problems, through the Levenberg-Marquardt approach.
}
\usage{
lsqnonlin(fun, x0, options = list(), ...)
}
\arguments{
  \item{fun}{User-defined, vector-valued function.}
  \item{x0}{starting point.}
  \item{...}{additional parameters passed to the function.}
  \item{options}{list of options, for details see below.}
}
\details{
  \code{lsqnonlin} computes the sum-of-squares of the vector-valued function
  \code{fun}, that is if \eqn{f(x) = (f_1(x), \ldots ,f_n(x))} then
  \deqn{min || f(x) ||_2^2 = min(f_1(x)^2 + \ldots + f_n(x)^2)}
  will be minimized.

  \code{x=lsqnonlin(fun,x0)} starts at point \code{x0} and finds a minimum
  of the sum of squares of the functions described in fun. \code{fun} shall
  return a vector of values and not the sum of squares of the values.
  (The algorithm implicitly sums and squares fun(x).)

  \code{options} is a list with the following components and defaults:
  \itemize{
    \item \code{tau}: used in starting value for Marquardt parameter.
    \item \code{tolx}: stopping parameter for step length.
    \item \code{tolg}: stopping parameter for gradient.
    \item \code{maxeval} the maximum number of function evaluations.
  }
  Typical values for \code{tau} are from \code{1e-6...1e-3...1} with small
  values for good starting points and larger values for not so good or known
  bad starting points.
}
\value{
  Returns a list with the following elements:
  \itemize{
    \item \code{x}: the point with least sum of squares value.
    \item \code{ssq}: the sum of squares.
    \item \code{ng}: norm of last gradient.
    \item \code{nh}: norm of last step used.
    \item \code{mu}: damping parameter of Levenberg-Marquardt.
    \item \code{neval}: number of function evaluations.
    \item \code{errno}: error number, corresponds to error message.
    \item \code{errmess}: error message, i.e. reason for stopping.
  }
}
\note{
  The refined approach, Fletcher's version of the Levenberg-Marquardt
  algorithm, may be added at a later time; see the references.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Madsen, K., and H. B.Nielsen (2010).Introduction to Optimization and
  Data Fitting. Technical University of Denmark, Intitute of Computer
  Science and Mathematical Modelling.

  Fletcher, R., (1971). A Modified Marquardt Subroutine for Nonlinear Least
  Squares. Report AERE-R 6799, Harwell.
}
\seealso{
  \code{\link{nls}}
}
\examples{
##  Rosenberg function as least-squares problem
x0  <- c(0, 0)
fun <- function(x) c(10*(x[2]-x[1]^2), 1-x[1])
lsqnonlin(fun, x0)

\dontrun{
#   Lanczos1 data (artificial data)
#   f(x) = 0.0951*exp(-x) + 0.8607*exp(-3*x) + 1.5576*exp(-5*x)
x <- linspace(0, 1.15, 24)
y <- c(2.51340000, 2.04433337, 1.66840444, 1.36641802, 1.12323249, 0.92688972,
       0.76793386, 0.63887755, 0.53378353, 0.44793636, 0.37758479, 0.31973932,
       0.27201308, 0.23249655, 0.19965895, 0.17227041, 0.14934057, 0.13007002,
       0.11381193, 0.10004156, 0.08833209, 0.07833544, 0.06976694, 0.06239313)
f2 <- function(b) y - b[1]*exp(-b[2]*x) + b[3]*exp(-b[4]*x) + b[5]*exp(-b[6]*x)
lsqnonlin(f2, c(1.2, 0.3, 5.6, 5.5, 6.5, 7.6))
}
}
\keyword{ fitting }
