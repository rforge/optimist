\name{eulerHeun}
\alias{eulerHeun}
\title{
  Euler-Heun ODE Solver
}
\description{
  Euler and Euler-Heun ODE solver.
}
\usage{
eulerHeun(f, a, b, y0, n = 100, improved = TRUE, ...)
}
\arguments{
  \item{f}{function in the differential equation \eqn{y' = f(x, y)}.}
  \item{a, b}{start and end points of the interval.}
  \item{y0}{starting value at a.}
  \item{n}{number of grid points.}
  \item{improved}{logical; shall the Heun method be used; default TRUE.}
  \item{...}{additional parameters to be passed to the function.}
}
\details{
  \code{eulerHeun} is an integration method for ordinary differential
  equations using the simple Euler resp. the improved Euler-Heun Method.
}
\value{
  List with components \code{t} for grid (or `time') points, and \code{y} 
  the vector of predicted values at those grid points.
}
\references{
  Quarteroni, A., and F. Saleri (). Scientific Computing with MATLAB
  and Octave. Second Edition, Springer-Verlag, Berlin Heidelberg, 2006.
}
\seealso{
  \code{\link{ode23}}
}
\examples{
## See Crank-Nicolson
}
\keyword{ ode }
