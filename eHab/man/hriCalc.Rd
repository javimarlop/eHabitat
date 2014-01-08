\name{hriCalc}
\alias{hriCalc}
\title{
Computes hri from spatial or raster objects
}
\description{
Functions that calculates similarities between habitat of interest and surroundings
and/or future.
}
\usage{
hriCalc(object, habitat, pvals, populationDensity) 
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{object}{
  A \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} or a 
  \code{\link[raster]{raster}} object with the similarities computed 
  by \code{\link{pHabitat}} or \code{\link{hri}}}
\item{habitat}{
  A \code{\link[sp]{SpatialPolygons}}* object with the boundaries of the 
  protected area or habitat of interest}
\item{pvals}{
  The probability limit(s) used to decide if a pixel
  could be a suitable replacement for the habitat or not}
\item{populationDensity}{
  A \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} or a 
  \code{\link{data.frame}} with the population density of the same 
  locations as used for the indicators. If not available, the population density
  is set equal to 1 within the protectedArea and 0 outside.}
}
\details{


\code{hriCalc} is usually called directly from \code{\link{hri}} or
\code{\link{mhri}}, but can also be called separately. It computes the relative
size of the area with a similarity above the threshold-value(s) defined by
\code{pvals}.




}
\value{
The result of a call to this function is a data.frame with one row for each
threshold, and the hri inside and outside the habitat for each of the similarity
layers. The relative size of the habitat is computed from the number of 
non-NA cells within the habitat.

}
\author{
Jon Olav Skoien
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{pHabitat}}
}
\examples{
library(eHab)
options(error = recover)
data(eHab)
hrep = mhri(indicators,protectedArea, pval = seq(0.1,0.9,0.1), nsim = 4, nboot = 4)
attr(hrep, "hri")
hriCalc(hrep, protectedArea, pval = seq(0.05,0.95,0.05))

}
\keyword{ spatial }
