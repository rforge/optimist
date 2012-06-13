\name{Test functions}
\alias{rastrigin}
\alias{rosenbrock}
\alias{wood}
\title{
  Optimization Test Function
}
\description{
  Simple and often used test function defined in higher dimensions.
  The dimension is determined by the length of the input vector.
}
\usage{
rastrigin(x)
rosenbrock(x)
wood(x)
}
\arguments{
  \item{x}{numeric vector of some length \code{n >= 2}.}
}
\details{
  The Rosenbrock test function is defined as
    \deqn{\sum_{i=1}^{n-1} (1-x_i)^2 + 100 (x_{i+1}-x_i^2)^2}
  and has a global minimum at \code{f(1,\ldots,1)=0} and local minima at
  \code{f(-1,1,\ldots,1)} for \code{n>=4}.
}
\value{
  Returns the value of the function at the specified point.
}
\references{
  See the Wikipedia.
}
\examples{
rosenbrock(c(1, 1, 1, 1))
}
