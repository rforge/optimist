\name{Nelder-Mead}
\alias{nelder_mead}
\title{
  Hooke-Jeeves Minimization Method
}
\description{
  An implementation of the Nelder-Mead algorithm for derivative-free
  optimization.
}
\usage{
nelder_mead(x0, f, maxiter = 1000, scale = 1,
            tol = .Machine$double.eps^(2/3), show = FALSE, ...)
}
\arguments{
  \item{x0}{starting vector.}
  \item{f}{nonlinear function to be minimized.}
  \item{maxiter}{max. number of iterations.}
  \item{scale}{scale factor, of -1 the maximum will be searched for.}
  \item{tol}{relative tolerance, to be used as stopping rule.}
  \item{show}{logical, whether to store intermediate information.}
  \item{\ldots}{additional arguments to be passed to the function.}
}
\details{
  Also called a `simplex' method for finding the local minimum of a function
  of several variables. The method is a pattern search that compares function
  values at the vertices of the simplex. The process generates a sequence of
  simplices with ever reducing sizes.
}
\value{
  List with following components:
    \item{xmin}{minimum solution found.}
    \item{fmin}{value of \code{f} at minimum.}
    \item{niter}{number of iterations performed.}
}
\references{
  J. H. Mathews and K. D. Fink (2006). Numerical Methods Using Matlab.
  Fourth Edition, Pearson Education.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\note{
  Used Matlab code as provided in the textbook by Mathews and Fink.

  For a more elaborate implementation of Nelder-Mead see the package
  `dfoptim'.
}
\seealso{
  \code{hooke_jeeves}
}
\examples{
##  Rosenbrock function
# nelder_mead(c(0, 0), rosenbrock)
# $xmin
# [1] 1 1
# $fmin
# [1] 3.951249e-19
# $niter
# [1] 13
}
\keyword{ optimize }
