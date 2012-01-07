\name{Wood Function}
\alias{wood}
\title{
  Wood Test Function
}
\description{
  Wood test function
}
\usage{
wood(x)
}
\arguments{
  \item{x}{numeric vector of length 4.}
}
\details{
  The Wood test function is defined as
    \deqn{wood(x_0,x_1,x_2,x_3) = 100 (x_1-x_0^2) + (1-x_0^2) + 
                                  90 (x_3-x_2^2) + (1-x_2^2) +
                                  10 (x_1+x_3-2)^2 + (x1-x3)^2 / 10}
  and has a global minimum at \code{f(1,1,1,1)=0}.

  \code{wood} is known to be difficult to handle by some optimization
  routines, especially with a starting point \code{x0=[-3,-1,-3,-1]}.
  Try, for example, `Nelder-Mead' in \code{optim}.
}
\value{
  Returns the value of the Wood function at the specified point.
}
\references{
  More, J. J., B. S. Grabow, and K. E. Hillstrom. Test unconstrained
  optimization software. ACM Trans. on Mathematical Software, Vol. 7,
  pp. 17--24, 1981.
  (TOMS algorithm 566, \url{http://www.netlib.org/toms/566})
}
\examples{
wood(c(0, 0, 0, 0))
wood(c(-3, -1, -3, -1))
}
