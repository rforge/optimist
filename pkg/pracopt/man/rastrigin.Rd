\name{Rastrigin Function}
\alias{rastrigin}
\title{
  Rastrigin Test Function
}
\description{
  Simple and often used test function defined in higher dimensions.
  The dimension is determined by the length of the input vector.
}
\usage{
rastrigin(x)
}
\arguments{
  \item{x}{numeric vector of some length \code{n >= 2}.}
}
\details{
  The Rastrigin test function is defined as
    \deqn{10^n + \sum_{i=1}^n x_i^2 - 10 \cos(2 pi x_i)}
  and has a global minimum at \code{f(0,\ldots,1)=0} and many local minima.
}
\value{
  Returns the value of the function at the specified point.
}
\references{
  See the Wikipedia.
}
\examples{
rastrigin(c(0, 0, 0, 0))
}
