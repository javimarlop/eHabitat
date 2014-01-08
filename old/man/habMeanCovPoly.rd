\name{habMeanCovPoly}
\alias{habMeanCovPoly}
\title{
Extract training data from the region of interest
}
\description{
Function for extracting and creating the training data from the region of interest, according to method.
}


\usage{
habMeanCovPoly(indicators, habitat, inHA = NULL, polyData = NULL, 
    method = "mahalanobis", mAPoints = 1000, na.rm = FALSE, bootsample = FALSE)
}

\arguments{
  \item{indicators}{
      A \code{\link[sp:SpatialPixels]{SpatialGridDataFrame}} with 
      the indicators used to estimate the similarity to the habitat}
  \item{habitat}{
      A \code{\link[sp]{SpatialPolygons}}*-object with the boundaries of the 
      protected area or habitat of interest, or a \code{\link[sp]{SpatialPoints}}
      object with precence locations. See also details below}
  \item{inHA}{array with indices of indicators that define the region of interest}
  \item{polyData}{matrix or data.frame with already extracted values of the indicators valid 
       for the habitat}
  \item{method}{Which method to use for computing similarites. The possibilities are
      at the moment "mahalanobis", "maxent", "bioclim" and "domain". The three last
      methods are mainly available through the \code{\link{dismo}}-package.}
  \item{mAPoints}{Number of pseudo-absence points for the Maxent method}
  \item{na.rm}{
      logical; defining whether NA-values in the indicator data set should be removed or
      not}
  \item{bootsample}{logical; whether a single bootstrap sample should be used instead.}
}

\details{

First of all, this  function extracts the training data for the region of interest. It will fail if it find 2 or less points, as no similarity method can be used for such small training data sets. 

The second part is separation of the categorical variables from the numerical variables. 

The third step is method dependent. For Mahalanobis distance, this means creating the mean vector and the covariance matrix. For the Maxent method, it means creation of a set of pseudo-absence locations. 

These are sampled randomly from all locations that are not included in the habitat.

For bootstrapping, only a single bootstrap sample is returned.


}
\value{
The result is a list including some or all of the following:

\code{polyData}: A data frame with the observations of habitat. If maxent, also the observations from the absence locations

\code{meanPoly}: The mean vector of the observations

\code{covPoly}: The covariance matrix of the observations

\code{inHA}: The indices of the habitat observations

\code{outHA}: The indices of the absence locations for Maxent, NULL for other methods

\code{pDat}: An array indicating presence and absence locations for Maxent, NULL for other methods

\code{nonNum}: A vector of column numbers that have been identified as categorical data


}


\author{
Jon Olav Skoien}


\seealso{
\code{\link{pHabitat}}
}
\examples{
data(eHab)
pDat = habMeanCovPoly(indicators,protectedArea)
}
\keyword{spatial}
