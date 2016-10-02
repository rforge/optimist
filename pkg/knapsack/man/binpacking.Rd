\name{binpacking}
\alias{binpacking}
\title{
  The Bin Packing Problem
}
\description{
  Solves the bin packing problem for integers
}
\usage{
binpacking(weights, cap, back = -1, jck = 0, lb = 0)
}
\arguments{
  \item{weights}{integer vector of weights of items}
  \item{cap}{capacity of all bins}
  \item{back}{number of backtrackings; default -1, ie. exact solution.}
  \item{jck}{input checking wanted; default 0 for no checks.}
  \item{lb}{lower bound; default 0 means: will be calculated.}
}
\details{
  Solves the bin packing problem for integer weights. It implements the
  algorithm described in section 8.5 of the book "Knapsack Problems" by
  S. Martello and P. Toth.

  \code{back} is -1 for exact solutions, 0 for approximate solutions, and
  positive integer restricts the number of backtrackings. Input checking
  is not necessary, as the R wrapper does this. \code{lb} will be set to
  the obvious lower bound \code{ceil(sum(weights)/cap)}.
}
\value{
  Returns a list with components \code{nbins} the minimum number of bins
  found for a solution, and \code{xbins} the number of the bin each item
  is assigned to.
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

# Medium-sized example provided by Eduardo Mesa
ws <- c(2083, 1836, 1736, 1536, 1236, 1136, 1086, 1036,  636,  536,  436)
qs <- c(126,    94,    6,  126,    4,   42,   82,    2,  282,    4,  132)
weights <- c()
for (i in 1:length(ws))
    weights <- c(weights, rep(ws[i], qs[i]))
capacity <- 6150

# optimal solution is nbins=169, lb>=167
res <- binpacking(weights, capacity, back = 0)
res$nbin                        # 171
}
