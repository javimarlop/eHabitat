\name{pHabitat}
\alias{pHabitat}
\title{
Estimate habitat replaceability 
}
\description{
Function for estimating similarities between the indicator values at different pixels  
and a region of interest (protected area, important bird area or presence only
observations).
}


\usage{
pHabitat(indicators = NULL, habitat = NULL, meanPoly = NULL, covPoly = NULL,
               polyData = NULL, meanCovPoly = NULL, na.rm = FALSE, indicators2 = NULL, 
               forecast = FALSE, nonNum = NULL, method = "mahalanobis", bootsample = FALSE)
}

\arguments{
  \item{indicators}{
      A \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} with 
      the indicators used to estimate the similarity to the habitat}
  \item{habitat}{
      A \code{\link[sp]{SpatialPolygons}}*-object with the boundaries of the 
      protected area or habitat of interest, or a \code{\link[sp]{SpatialPoints}}
      object with precence locations. See also details below}
  \item{meanPoly}{Mean of the indicators at the locations of \code{habitat}, 
      should only be given if already computed or if it can not be found from
      habitat and indicators}
  \item{covPoly}{Covariance matrix of indicators at the locations of \code{habitat}
        should only be given if already computed or if it can not be found from
      habitat and indicators. See details below for usage.}
  \item{polyData}{matrix or data.frame with the values of the indicators valid 
       for the habitat. See details below for usage.}
  \item{meanCovPoly}{list of covariance matrix, mean, polyData etc from a call
       to the internal function eHab:::habMeanCovPoly. See details below for usage.}
  \item{na.rm}{
      Logical; defining whether NA-values in the indicator data set should be removed or
      not, or whether they should be interpolated (\code{na.rm = "int")}.}
  \item{indicators2}{A second \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} 
      with forecasted indicators}
  \item{forecast}{Logical; if TRUE and neither \code{indicators2} nor 
      \code{covPoly} and \code{meanPoly} are not submitted,
      the data of the indicators will be split in two columnwise, assuming that the
      first part of the data frame represents the current state and the second part 
      the forecast. If covPoly and meanPoly are submitted. If \code{forecast} is TRUE,
      also the area within the habitat will be included for computation of the similarity.}
  \item{nonNum}{list of non-numerical inputs, usually from a call to ehab:::habMeanCovPoly.}
  \item{method}{Which method to use for computing similarites. The possibilities are
      at the moment "mahalanobis", "maxent", "bioclim" and "domain". The three last
      methods are mainly available through the \code{\link{dismo}}-package.}
  \item{bootsample}{Whether bootsampling should be done.}
}

\details{
The function computes similarities
to a region of interest from the surroundings and/or future. Different methods
are available, such as: "mahalanobis", "maxent", "bioclim" and "domain". The last
three are available through the package \code{\link[dismo]{dismo}}, maxent requires
also installation of the java package maxent.jar in the java directory of dismo.


The habitat will normally refer to the habitat or the presence only observations
of one species, or the \code{\link[sp]{SpatialPolygons}}*-object defining a region 
of interest, typically a protected area.

There are four different ways of submitting the habitat characteristics.
\itemize{
 \item{indicators} {Let the function find the habitat characteristics from
         the set of indicators. A habitat must be submitted.}
 \item{meanCovPoly} {The result of an earlier call to \code{pHabitat}, stored in
         found as \code{attr(pHabImg, "meanCovPoly")} if \code{pHabImg} is the  
         result of a call to \code{link{pHabitat}} or \code{link{hri}}.}
 \item{meanPoly and covPoly} {The means and the covariance matrix that describes
         the characteristics of the habitat. This is only useful if the method
         is "mahalanobis".}
  \item{polyData} {A matrix or a data.frame of the habitat characteristics}
 }

If any of the last three options are chosen, the indicators are not used. 
It is only possible to choose one of the last three options, combinations 
will cause an error.

If bootsampling is enabled, the similarity will instead be computed from a
bootstrap sample. 

}
\value{
The function returns a \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} 
equal to \code{indicators} but extended with two columns:
\item{mDist}{Mahalanobis distance between the pixel and the average of the 
habitat/observation locations}
 \item{pHab}{The probability of a pixel being a suitable habitat }
 
The object will also have a range of attributes:

\code{attr(spdf,"mean") = meanPoly} The mean values of indicators used in 
computation of similarity.

\code{attr(spdf,"cov") = covPoly} The covariance matrix of indicators used
in computation of similarity.

\code{attr(spdf,"inHA") = inHA} An array indicating which rows from the 
      \code{data.frame} of \code{indicators} are within the habitat

\code{attr(spdf,"polyData") = polyData} A data.frame with the values of the 
       indicators for the selected habitat.

\code{attr(spdf, "hist") = hist} A histogram of the values of the indicators 
       for the selected habitat

\code{attr(spdf, "meanCovPoly")} A list with information about the habitat
      characteristics, and, if necessary, information about categorical data,
      constant values etc.
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
pp = pHabitat(indicators,protectedArea)
}
\keyword{spatial}
