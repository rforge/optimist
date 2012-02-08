\name{subsetsum}
\alias{subsetsum}
\title{
  The Subset Sum Problem
}
\description{
  Solves the subset sum problem for integers.
}
\usage{
subsetsum(weights, target, check = TRUE)
}
\arguments{
  \item{weights}{integer vector of weights.}
  \item{target}{target sum to be reached, integer.}
  \item{check}{logical; Fortran check routine enabled; cannot be turned off.}
}
\details{
  Solves the subset sum problem for integer weights. It implements the mixed
  algorithm described in section 4.2.3 of the book``Knapsack Problems'' by
  S. Martello and P. Toth.

  The subset sum problem can be described as follows: Given integer weights
  \code{w[j], j=1,...,n} and a target value \code{W}, find a subset of weights,
  defined by a 0-1 vector \code{x[j]}, such that

  Maximize \code{Z = w[1]*x[1] + ... + w[n]*x[n]}

  subject to: \code{w[1]*x[1] + ... + w[n]*x[n] <= W}

  The input problem description must satisfy the following conditions:
  \itemize{
  \item maximum weight is smaller than the target
  \item sum of all weights is greater or equal to the target
  }
}
\value{
  A list with \code{ssum} the maximal summation reached, and \code{inds}
  the indices of weights.
}
\author{
  HwB  email: <hwborchers@googlemail.com>

  The Fortran routines used are copyright of S. Martello and P. Toth and are
  distributed under a free license only for personal research and academic use.
}
\references{
  Martello, S., and P. Toth (1990). ``Knapsack Problems: Algorithms and
  Computer Implementations''. John Wiley & Sons, Ltd.
}
\seealso{
  \code{adagio::subsetsum}
}
\examples{
##  Example from R-help
amount <- 4748652
products <- c(
    1198000, 1110000, 1031655, 829350, 800000, 533600, 456000, 123289,
      90900,   83800,   80000,  77382,  77000,  76950,  71190,  71040,
      70740,   70740,   63051,  59700,  59500,  55000,  53730,  51186,
      49600,   48048,   42000,  42000,  42000,  42000,  42000,  42000,
      42000,   42000,   42000,  42000,  40000,  35100,  35000,  33000,
      33000,   30500,   30500,  30500,  30500,  28000,  28000,  28000,
      28000,   28000,   27465,  27000,  27000,  27000,  27000,  27000,
      27000,   26140,   25908,  21716,  15025,  15025,  15025,  15025)
# find a subset that sums up to amount:
R <- subsetsum(products, amount)
R$ssum == amount            # TRUE
R$inds                      # 1  2  3  4  8  9 10 12 24 25 40 46 51 64

##  Some weights are greater than the target
n <- length(products)
amount <- 1000000
imax <- which(products > amount)
prods <- products[products <= amount]

# reproduce the indices
R <- subsetsum(prods, amount)
inds <- setdiff(1:n, imax)[R$inds]              # 4 25 39 47 49 61 64
sum(products[inds]) == amount                   # TRUE
}

