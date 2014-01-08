\name{pHabitatUnc}
\alias{pHabitatUnc}
\title{
Estimate habitat replaceability with uncertainty using an analytical solution
}
\description{
Function for estimating similarities between the indicator values at different pixels  
and a region of interest (protected area, important bird area or presence only
observations), using the Mahalanobis distance.
}


\usage{
pHabitatUnc(indicators = NULL, habitat = NULL, lmc)
}

\arguments{
  \item{indicators}{
      A \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} with 
      the indicators used to estimate the similarity to the habitat}
  \item{habitat}{
      A \code{\link[sp]{SpatialPolygons}}*-object with the boundaries of the 
      protected area or habitat of interest, or a \code{\link[sp]{SpatialPoints}}
      object with precence locations. See also details below}
  \item{lmc}{Linear model of coregionalization describing the uncertainty of the 
      indicators, typcailly a result from a call 
      to \code{\link[gstat:fit.lmc]{fit.lmc}}.}
}

\details{
The function computes similarities
to a region of interest from the surroundings.

The uncertainty of the indicators is propagated using an analytical solution. 
This is still in a test phase, and it is quite common that the solution is
not stable. 

}
\value{
The function returns a \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} 
equal to \code{indicators} but extended with two columns:
\item{mDist}{Mahalanobis distance between the pixel and the average of the 
habitat/observation locations}
 \item{pHab}{The probability of a pixel being a suitable habitat }
 
}


\references{
Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of female 
black bear habitat use for a geographic information system. Journal of Wildlife 
Management, 57, 519-526. }
\author{
Jon Olav Skoien}


\seealso{
\code{\link{hri}}
}
\examples{
data(eHab)
unc = indicators
unc@data = unc@data/10
lmc = eHab:::makeLmc(unc)
pp = pHabitatUnc(indicators,protectedArea, lmc)
}
\keyword{spatial}
