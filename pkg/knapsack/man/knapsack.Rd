\name{knapsack}
\alias{knapsack}
\title{
  The Single Knapsack Problem
}
\description{
  Solves the 0-1 (single), bounded, and unbounded single knapsack problem.
}
\usage{
knapsack(profits, weights, capacity, bounds = NULL, check = TRUE)
}
\arguments{
  \item{profits}{integer vector of profits, with length greater or equal 2.}
  \item{weights}{integer vector of weights, of the same length.}
  \item{capacity}{total capacity of the knapsack, integer.}
  \item{bounds}{bounds on the multiplicity of items; can be \code{NULL},
                \code{Inf}, or an integer vector.}
  \item{check}{logical; Fortran check routine enabled; cannot be turned off.}
}
\details{
  Solves the 0-1 single knapsack problem for integer profits and weights
  when \code{is.null(bounds)}, the default.
  It implements the branch-and-bound algorithm described in section 2.5.2
  of Martello and Toth's book ``Knapsack Problems''.

  With \code{bounds} an integer vector it solves the bounded knapsack problem 
  for integer profits and weights, i.e. \eqn{x_j \le b_j} for all components.
  If \code{bounds} is a single integer, it will be expanded to the length of
  \code{profits} (and \code{weights}).
  It implements the transformation method (that is, transformed into an 
  equivalent 0-1 knapsack problem) described in section 3.2 of Martello
  and Toth's book.

  With \code{bounds==Inf} it solves the unbounded knapsack problem, that is
  the components \code{x_j} can get as large as possible.
  It implements the enumerative algorithm described in section 3.6.3 of the
  book ``Knapsack Problems''.

  A mixed bounded/unbounded knapsack problem can be formulated by mixing
  \code{Inf} with integers in one vector as \code{bounds}.

  Problem formulation: Given a set of \eqn{n} items and a knapsack, with

  \eqn{p_j = profit} of item \eqn{j},

  \eqn{w_j = weight} of item \eqn{j},

  \eqn{c = capacity} of the knapsack,

  select a subset of the items, described by the 0/1 vector \eqn{x}, so as to

  \bold{maximize} \eqn{z = \sum_{j=1}^n p_j x_j}

  \strong{subject to} \eqn{\sum_{j=1}^n w_j x_j \le c}

  where \eqn{x_j = 0} or 1 for \eqn{j = 1,\ldots,n}.

  The input problem description must satisfy the following conditions:
  \itemize{
  \item maximum weight is smaller than the capacity
  \item sum of all weights is greater or equal to the capacity
  \item \code{profits[j]/weights[j], j=1..(n-1)} must be monotonically 
  decreasing.
  }
}
\value{
  Vector of indices, with components only 0 or 1.
}
\note{
  This routines do not work correctly if too many or all of the elements
  in the sequence \code{profit/weight} are equal instead of decreasing.
  Thus, it cannot be usedfor, e. g., for the subset sum problem. This will
  be corrected in the future.
}
\author{
  HwB  email: <hwborchers@googlemail.com>
}
\references{
  Martello, S., and P. Toth (1990). ``Knapsack Problems: Algorithms and
  Computer Implementations''. John Wiley & Sons, Ltd.
}
\seealso{
  Other packages implementing knapsack routines.
}
\examples{
## Example 1
p <- c(15, 100, 90, 60, 40, 15, 10,  1)
w <- c( 2,  20, 20, 30, 40, 30, 60, 10)
cap <- 102
(is <- knapsack(p, w, cap))
# [1] 1 2 3 4 6 , capacity 102 and total profit 280

## Example 2
p <- c(70, 20, 39, 37, 7, 5, 10)
w <- c(31, 10, 20, 19, 4, 3,  6)
cap <- 50
(is <- knapsack(p, w, cap))
# [1] 1 4 , capacity 50 and total profit 107

##  Example 3: 2012 as sum of a minimal number of squares
w <- (1:44)^2; p <- w-1; p[1] <- 1
knapsack(p, w, 2012)
# 1^2 + 7^2 + 21^2 + 39^2 == 2012 # as capacity is fully eaten up.

##  Example 4: <rosettacode.org/wiki/Knapsack_problem/0-1>
K <- read.table(textConnection(
"item                   weight  value   pieces
map                       9     150     1
compass                  13      35     1
water                   153     200     2
sandwich                 50     160     2
glucose                  15      60     2
tin                      68      45     3
banana                   27      60     3
apple                    39      40     3
cheese                   23      30     1
beer                     52      10     3
suntan_cream             11      70     1
camera                   32      30     1
T-shirt                  24      15     2
trousers                 48      10     2
umbrella                 73      40     1
waterproof_trousers      42      70     1
waterproof_overclothes   43      75     1
note-case                22      80     1
sunglasses                7      20     1
towel                    18      12     2
socks                     4      50     1
book                     30      10     2"), header=TRUE, row.names="item")

ks <- knapsack(K$value, K$weight, 400)    # profit: 1030, capacity: 396
is <- sort(ks$indices)
row.names(K)[is]
#  [1] "map"                "compass"                "water"                 
#  [4] "sandwich"           "glucose"                "banana"                
#  [7] "suntan_cream"       "waterproof_trousers"    "waterproof_overclothes"
# [10] "note-case"          "sunglasses"             "socks"                 

##  Example 5: Bounded knapsack problem
ks <- knapsack(K$value, K$weight, 400, K$pieces)
                                          # profit: 1200, capacity: 385
for (i in 1:length(ks$indices))
    cat(row.names(K)[ks$indices[i]], "\t", ks$nitems[i], "\n")
# map 	 1 
# socks 	 1 
# suntan_cream 	 1 
# glucose 	 2 
# note-case 	 1 
# sandwich 	 2 
# sunglasses 	 1 
# compass 	 1 
# banana 	 3 
# waterproof_overclothes 	 1 
# waterproof_trousers 	 1 
# cheese 	 1 

##  Example 6: Unbounded knapsack
p <- c(20, 39, 52, 58, 31, 4, 5)
w <- c(15, 30, 41, 46, 25, 4, 5)
knapsack(p, w, 101, bounds = Inf)
# indices: 1 3 ; nitems: 4 1 ; profit: 132 ; capacity: 101
}

