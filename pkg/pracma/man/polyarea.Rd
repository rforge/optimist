\name{polyarea}
\alias{polyarea}
\alias{polycenter}
\alias{polylength}
\title{Area of a Polygon}
\description{
  Calculates the area and length of a polygon given by the vertices in the 
  vectors \code{x} and \code{y}.
}
\usage{
  polyarea(x, y)

  polylength(x, y)
  polycenter(x, y)
}
\arguments{
  \item{x}{x-coordinates of the vertices defining the polygon}
  \item{y}{y-coordinates of the vertices}
}
\details{
  \code{polyarea} calculates the area of a polygon defined by the vertices
  with coordinates \code{x} and \code{y}. Areas to the left of the vertices
  are positive, those to the right are counted negative.

  The computation is based on the Gauss polygon area formula. The polygon
  automatically be closed, that is the last point need not be / should not
  be the same as the first.

  If some points of self-intersection of the polygon line are not in the
  vertex set, the calculation will be inexact. The sum of all areas will be
  returned, parts that are circulated in the mathematically negative sense
  will be counted as negative in this sum.

  If \code{x}, \code{y} are matrices of the same size, the areas of all
  polygons defined by corresponding columns are computed.

  \code{polycenter} calculates the center (of mass) of the figure defined by
  the polygon. Self-intersections should be avoided in this case.
  The mathematical orientation of the polygon does not have influence on the
  center coordinates.

  \code{polylength} calculates the length of the polygon
}
\value{
  Area or length of the polygon resp. sum of the enclosed areas; or the
  coordinates of the center of gravity.
}
\seealso{
  \code{\link{trapz}}, \code{\link{arclength}}
}
\examples{
  # Zu Chongzhi's calculation of pi (China, about 480 A.D.),
  # approximating the circle from inside by a regular 12288-polygon(!):
  phi <- seq(0, 2*pi, len=3*2^12+1)
  x <- cos(phi)
  y <- sin(phi)
  pi_approx <- polyarea(x, y)
  print(pi_approx, digits=8)    #=> 3.1415925 or 355/113

  polylength(x, y)              #=> 6.2831852 where 2*pi is 6.2831853

  x1 <- x + 0.5; y1 <- y + 0.5
  x2 <- rev(x1); y2 <- rev(y1)
  polycenter(x1, y1)            #=> 0.5 0.5
  polycenter(x2, y2)            #=> 0.5 0.5
}
\keyword{ math }
