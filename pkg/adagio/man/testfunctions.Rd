\name{Testfunctions}
\alias{fnRosenbrock}
\alias{grRosenbrock}
\alias{fnNesterov}
\alias{grNesterov}
\alias{fnRastrigin}
\alias{grRastrigin}
\alias{fnHald}
\alias{grHald}
\title{
  Optimization Test Function
}
\description{
  Simple and often used test function defined in higher dimensions and with
  analytical gradients, especially suited for performance tests. Analytical 
  gradients are provided with the \code{gr} prefix.
  The dimension is determined by the length of the input vector.
}
\usage{
fnRosenbrock(x)
grRosenbrock(x)
fnNesterov(x)
grNesterov(x)

fnRastrigin(x)
grRastrigin(x)
fnHald(x)
grHald(x)
}
\arguments{
  \item{x}{numeric vector of some length \code{n >= 2}.}
}
\details{
\code{fnRosenbrock} is Rosenbrock's famous valley function from 1960. It can also be seen as a least-squares problem:
\deqn{\sum_{i=1}^{n-1} (1-x_i)^2 + 100 (x_{i+1}-x_i^2)^2}
\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: \tab at f(-1, 1, \ldots, 1) for n >= 4 \cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 1, i = 1:n \cr
}

\code{fnNesterov} is Nesterov's smooth adaptation of Rosenbrock, based on the idea of Chebyshev polynomials. This function is even more difficult to optimize than Rosenbrock's:
\deqn{(x_1 - 1)^2 / 4 + \sum_{i=1}^{n-1} (1 + x_{i+1} - 2 x_i^2)}
\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: ?\cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 1, i = 1:n \cr
}

Rastrigin's function is a famous, non-convex example from 1989 for global optimization. It is a typical example of a multimodal function with many local minima:
  \deqn{10 n + \sum_1^n (x_i^2 - 10 \cos(2 \pi x_i))}
\tabular{ll}{
  No. of Vars.: \tab n >= 2 \cr
  Bounds: \tab -5.12 <= xi <= 5.12 \cr
  Local minima: \tab many \cr
  Minimum: \tab 0.0 \cr
  Solution: \tab xi = 0, i = 1:n \cr
}

Hald's function is a typical example of a non-smooth test function, from Hald and Madsen in 1981.
\deqn{\max_{1 \le i \le n} \frac{x_1 + x_2 t_i}{1 + x_3 t_i + x_4 t_i^2 + x_5 t_i^3} - \exp(t_i)}
where \eqn{t_i = -1 + (i - 1)/10} for \eqn{1 \le i \le 21}.
\tabular{ll}{
  No. of Vars.: \tab n =5 \cr
  Bounds: \tab -1 <= xi <= 1 \cr
  Local minima: \tab ? \cr
  Minimum: \tab 0.0001223713 \cr
  Solution: \tab (0.99987763,  0.25358844, -0.74660757,  0.24520150, -0.03749029) \cr
}
}
\value{
  Returns the value of the function at the specified point.
}
\references{
  See the Wikipedia.
}
\examples{
x <- runif(5)
fnRosenbrock(x); grRosenbrock(x)
fnNesterov(x);   grNesterov(x)
fnRastrigin(x);  grRastrigin(x)
fnHald(x);       grHald(x)
}
