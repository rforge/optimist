\name{maxsub}
\alias{maxsub}
\alias{maxsub_r}

\title{
Maximal Subarray
}
\description{
Find a subarray with maximal positive sum.
}
\usage{
maxsub(x, inds = TRUE)
maxsub_r(x, inds = FALSE)
}
\arguments{
  \item{x}{numeric vector.}
  \item{inds}{logical; shall the indices be returned?}
}
\details{
Find a contiguous subarray whose sum is maximally positive.

\code{maxsub_r} is an R script kept here for documentation and tests, while
\code{maxsub} is a compiled and very fast version with a running time of
\code{O(n)} where \code{n} is the length of the input vector \code{x}.
}
\value{
Either just a maximal sum, or a list this sum as component \code{sum} plus
the start and end indices as a vector \code{inds}.
}
\references{
Jon Bentley, ``Programming Pearls''
}
\author{
  HwB  <hwborchers@googlemail.com>
}
\examples{
\dontrun{
    set.seed(8237)
    x <- rnorm(1e6)
    system.time(res <- maxsub(x, inds = TRUE))
    res}
}
\keyword{ optimization }
