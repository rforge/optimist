\name{Rosser}
\alias{Rosser}
\title{Rosser Matrix}
\description{
Generate the Rosser matrix.
}
\usage{
Rosser()
}
\details{
This is a classic symmetric eigenvalue test problem.
It has a double eigenvalue, three nearly equal eigenvalues, dominant
eigenvalues of opposite sign, a zero eigenvalue, and a small, nonzero
eigenvalue.
}
\value{
matrix of size 8 x 8
}
\seealso{
  \code{\link{Wilkinson}}
}
\examples{
Rosser()
}
\keyword{ array }
