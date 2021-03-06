\name{Primes}
\alias{Primes}
\alias{atkin_sieve}
\title{Prime Numbers}
\description{
  Eratosthenes resp. Atkin sieve methods to generate a list of prime numbers
  less or equal \code{n}, resp. between \code{n1} and \code{n2}.
}
\usage{
  Primes(n1, n2 = NULL)

  atkin_sieve(n)
}
\arguments{
  \item{n, n1, n2}{natural numbers with \code{n1 <= n2}.}
}
\details{
  The list of prime numbers up to \code{n} is generated using the "sieve of
  Eratosthenes". This approach is reasonably fast, but may require a lot of
  main memory when \code{n} is large.

  \code{Primes} computes first all primes up to \code{sqrt(n2)} and then
  applies a refined sieve on the numbers from \code{n1} to \code{n2}, thereby
  drastically reducing the need for storing long arrays of numbers.

  The sieve of Atkins is a modified version of the ancient prime number sieve 
  of Eratosthenes. It applies a modulo-sixty arithmetic and requires less
  memory, but in R is not faster because of a double for-loop.

  In double precision arithmetic integers are represented exactly only up to
  2^53 - 1, therefore this is the maximal allowed value.
}
\value{
  vector of integers representing prime numbers
}
\references{
  A. Atkin and D. Bernstein (2004), Prime sieves using quadratic forms.
  Mathematics of Computation, Vol. 73, pp. 1023-1030.
}
\seealso{
  \code{\link{isPrime}}, \code{gmp::factorize}, \code{pracma::expint1}
}
\examples{
Primes(1000)
Primes(1949, 2019)

atkin_sieve(1000)

\dontrun{
##  Appendix:  Logarithmic Integrals and Prime Numbers (C.F.Gauss, 1846)

library('gsl')
# 'European' form of the logarithmic integral
Li <- function(x) expint_Ei(log(x)) - expint_Ei(log(2))

# No. of primes and logarithmic integral for 10^i, i=1..12
i <- 1:12;  N <- 10^i
# piN <- numeric(12)
# for (i in 1:12) piN[i] <- length(primes(10^i))
piN <- c(4, 25, 168, 1229, 9592, 78498, 664579,
         5761455, 50847534, 455052511, 4118054813, 37607912018)
cbind(i, piN, round(Li(N)), round((Li(N)-piN)/piN, 6))

#  i     pi(10^i)      Li(10^i)  rel.err  
# --------------------------------------      
#  1            4            5  0.280109
#  2           25           29  0.163239
#  3          168          177  0.050979
#  4         1229         1245  0.013094
#  5         9592         9629  0.003833
#  6        78498        78627  0.001637
#  7       664579       664917  0.000509
#  8      5761455      5762208  0.000131
#  9     50847534     50849234  0.000033
# 10    455052511    455055614  0.000007
# 11   4118054813   4118066400  0.000003
# 12  37607912018  37607950280  0.000001
# --------------------------------------}
}
