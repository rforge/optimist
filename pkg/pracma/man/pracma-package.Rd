\name{pracma-package}
\alias{pracma-package}
\alias{pracma}
\docType{package}
\title{
Practical Numerical Math Routines
}
\description{
  This package provides R implementations of more advanced
  math functions from Matlab and Octave (and the Euler Math Toolbox)
  with a special view on optimization and time series routines.

  Some of these realizations are the result of courses on Scientific
  Computing (``Wissenschaftliches Rechnen'') and are mostly intended to
  demonstrate how to implement certain algorithms in R/S.
}
\details{
\tabular{ll}{
Package: \tab pracma\cr
Type: \tab Package\cr
Version: \tab 0.3-1\cr
Date: \tab 2011-04-17\cr
License: \tab GPLv3\cr
LazyLoad: \tab yes\cr
}
  Several implementations of root finding routines and numerical sorting
  are provided as examples.
}
\author{
Hans W Borchers

Maintainer: HwB  <hwborchers@googlemail.com>
}
\references{
  Abramowitz, M., and I. A. Stegun (1965). Handbook of Mathematical Functions
  (with Formulas, Graphs, and Mathematical Tables). Dover, New York.

  Cormen, Th. H., Ch. E. Leiserson, and R. L. Rivest (2009). Introduction
  to Algorithms. Third Edition, The MIT Press, Cambridge, MA.

  Gentle, J. E. (2009). Computational Statistics.
  Springer Science+Business Media LCC, New York.

  Hazewinkel, M., Editor (2002). Encyclopaedia of Mathematics.
  Springer-Verlag, Berlin Heidelberg New York,
  Internet: \url{http://eom.springer.de/}.

  Press, W. H., S. A. Teukolsky, W. T Vetterling, and B. P. Flannery (2007).
  Numerical Recipes: The Art of Numerical Computing. Third Edition, incl.
  Numerical Recipes Software, Cambridge University Press, New York.

  Quarteroni, A., and F. Saleri (2006). Scientific Computing with Matlab
  and Octave. Second Edition, Springer-Verlag, Berlin Heidelberg.

  Quarteroni, A., R. Sacco, and F. Saleri (2007). Numerical Mathematics.
  Second Edition, Springer-Verlag, Berlin Heidelberg.

  Skiena, St. S. (2008). The Algorithm Design Manual. Second Edition,
  Springer-Verlag, London.

  Stoer, J., and R. Bulirsch (2002). Introduction to Numerical Analysis.
  Third Edition, Springer-Verlag, New York.

  Strang, G. (2007). Computational Science and Engineering.
  Wellesley-Cambridge Press.
}
\keyword{ package }
\seealso{
  The R package `matlab' contains some of the basic routines from Matlab,
  but unfortunately not any of the higher math routines.
}
\examples{
\dontrun{
tail(primes(1e7))  # the last 5 prime numbers below 10,000,000 .}
}
