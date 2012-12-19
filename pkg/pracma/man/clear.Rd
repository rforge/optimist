\name{clear, what, who}
\alias{clear}
\alias{ver}
\alias{what}
\alias{who}
\alias{whos}
\title{
Clear function (Matlab style)
}
\description{
  List or remove items from workspace, or display system information.
}
\usage{
clear(lst)
ver()

what(dname = getwd())
who()
whos()
}
\arguments{
  \item{lst}{Character vector of names of variables in the global environment.}
  \item{dname}{(relative or absolute) directory path.}
}
\details{
  Remove these or all items from the workspace, i.e. the global environment,
  and freeing up system memory.

  \code{what()} lists all files in a directory.\cr
  \code{who()}  lists all items on the workspace.\cr
  \code{whos()} lists all items and their class and size.

  \code{ver()} displays version and license information for R and all the
  loaded packages.
}
\value{
  Invisibly NULL.
}
\seealso{
  \code{\link{ls}}, \code{\link{rm}}, \code{\link{sessionInfo}}
}
\examples{
# clear()  # DON'T
# what()
# who()
# whos()
# ver()
}
