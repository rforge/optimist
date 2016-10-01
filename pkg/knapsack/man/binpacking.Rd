\name{binpacking}
\alias{binpacking}
\title{
  The Bin Packing Problem
}
\description{
  Solves the bin packing problem for integers
}
\usage{
binpacking(weights, cap)
}
\arguments{
  \item{weights}{integer vector of weights of items}
  \item{cap}{capacity of bins}
}
\details{
Solves the bin packing problem for integer weights. It implements the
algorithm described in section 8.5 of the book "Knapsack Problems" by
S. Martello and P. Toth.
}
\value{
  Returns the minimum number of bins found as a solution.
}
\references{
Martello, S., and P. Toth (1990). "Knapsack Problems: Algorithms and
Computer Implementations". John Wiley & Sons, Ltd.
}
\author{
The Fortran routines used are copyright of S. Martello and P. Toth and are
distributed under a free license only for personal research and academic use.
}
\note{
This routine will be slow even for medium-sized problems.
}
\seealso{
adagio::binpacking
}
\examples{
wts <- 20:1
cap <- 21
binpacking(wts, cap)

# First example from the following Web page:
# www2.wiwi.uni-jena.de/Entscheidung/binpp/bin1dat.htm
ws <- c(99,99,96,96,92,92,91,88,87,86,
        85,76,74,72,69,67,67,62,61,56,
        52,51,49,46,44,42,40,40,33,33,
        30,30,29,28,28,27,25,24,23,22,
        21,20,17,14,13,11,10, 7, 7, 3)
cp <- 100
binpacking(ws, cp)
# nbins
# [1] 25
# xbins
#  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 20 21 18 22
# [22] 19 19 22 18 23 23 24 16 17 15 20 21 14 24 24 13 12 25 25 25 25
# [43] 23 10  9  8 11  5  6  3
}
