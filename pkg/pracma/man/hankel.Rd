\name{Hankel}
\alias{Hankel}
\title{Hankel Matrix}
\description{
Generate Hankel matrix from column and row vector
}
\usage{
  Hankel(a, b)
}
\arguments{
  \item{a}{vector that will be the first column}
  \item{b}{vector that if present will form the last row.}
}
\details{
\code{Hankel(a)} returns the square Hankel matrix whose first column is
\code{a} and whose elements are zero below the secondary diagonal. (I.e.:
\code{b} may be missing.)

\code{Hankel(a, b)} returns a Hankel matrix whose first column is \code{a}
and whose last row is \code{b}. If the first element of \code{b} differs
from the last element of \code{a} it is overwritten by this one.
}
\value{
  matrix of size \code{(length(a), length(b))}
}
\seealso{
  \code{\link{Toeplitz}}, \code{\link{Hadamard}}
}
\examples{
Hankel(1:5, 5:1)
}
\keyword{ array }
