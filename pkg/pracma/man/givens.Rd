\name{givens}
\alias{givens}
\alias{qrGivens}
\title{Givens Rotation}
\description{
  Givens Rotations and QR decomposition
}
\usage{
givens(xk, xl)
qrGivens(A)
}
\arguments{
  \item{xk, xl}{real numbers.}
  \item{A}{numeric square matrix.}
}
\details{
  \code{givens(xk,xl)} generates a unitary 2-by-2 matrix \code{U} such that
  \code{U * [xk;xl] = [x,0]} where \code{x=sqrt(xk^2+xl^2)}.

  \code{qrGivens(A)} returns a QR decomposition (or factorization) of the
  square matrix \code{A}.
}
\value{
  List with two matrices \code{Q} and \code{R}, \code{Q} orthonormal and
  \code{R} upper triangular, such that \code{A=Q\%*\%R}.
}
\author{
  hwb \email{hwborchers@googlemail.com}
}
\references{
  Golub, G. H., and Ch. F. van Loan (1996). Matrix Computations.
  Third edition, John Hopkins University Press, Baltimore.
}
\seealso{
  \code{\link{qrHouseholder}}
}
\examples{
##  QR decomposition
A <- matrix(c(0,-4,2, 6,-3,-2, 8,1,-1), 3, 3, byrow=TRUE)
gv <- qrGivens(A)
(Q <- gv$Q); (R <- gv$R)
zapsmall(Q \%*\% R)

qrGivens(magic(5))

##  Solve a linear system of equations
b <- c(-2, -6, 7)
m <- nrow(A); n <- ncol(A)
b <- t(Q[, 1:n]) \%*\% b
x <- numeric(n)
x[n] <- b[n] / R[n, n]
for (k in (n-1):1)
    x[k] <- (b[k] - R[k, (k+1):n] \%*\% x[(k+1):n]) / R[k, k]
A \%*\% x
}
\keyword{ array }
