\name{wood}
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
}
\value{
  Returns the value of the Wood function at the specified point.
}
\references{
  Unknown
}
\examples{
wood(c(0, 0, 0, 0))
}
