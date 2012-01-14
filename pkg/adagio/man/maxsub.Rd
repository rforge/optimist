\name{maxsub}
\alias{maxsub}

\title{
Maximal Sum Subarray
}
\description{
Find a subarray with maximal positive sum.
}
\usage{
maxsub(x, inds = TRUE, compiled = TRUE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{inds}{logical; shall the indices be returned?}
  \item{compiled}{logical; shall the compiled version be used?}
}
\details{
Find a contiguous subarray whose sum is maximally positive.

\code{maxsub} can use a compiled and very fast version with a running time of
\code{O(n)} where \code{n} is the length of the input vector \code{x}.
}
\value{
Either just a maximal sum, or a list this sum as component \code{sum} plus
the start and end indices as a vector \code{inds}.
}
\references{
Bentley, Jon (1986). ``Programming Pearls'', Column 7.
Addison-Wesley Publ. Co., Reading, MA.
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\examples{
\dontrun{
    set.seed(8237)
    x <- rnorm(1e6)
    system.time(res <- maxsub(x, inds = TRUE, compiled = TRUE))
    res}
}
\keyword{ optimization }
