\name{Psi}
\alias{sp.psi}
\title{
  Psi Function
}
\description{
  The Psi function is defined as 
  \deqn{\psi(z) = \frac{d}{dz} \log \Gamma(z) = \frac{\Gamma'(z)}{\Gamma(z)}}
  which is also called the digamma function.
}
\usage{
sp.psi(z)
}
\arguments{
  \item{z}{real or complex argument.}
}
\details{
  The computation utilizes an asymptotic expansion and a recurrence relation.
}
\value{
  Returns the value of the Psi function.
}
\examples{
Psi <- Vectorize(sp.psi)
print(Psi(1), digits = 16)
# [1] -0.5772156649015329               # Euler's Gamma constant
print(Psi(1/2), digits = 16)
# [1] -1.963510026021424                # -Gamma - 2*log(2)

\dontrun{
x <- seq(-5, 5, length.out = 501)
y <- Psi(x)
plot(x, y, type = 'l', col = "blue",
        main = "Psi function")
grid()}
}
