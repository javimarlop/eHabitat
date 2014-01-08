\name{ecohri}
\alias{ecohri}
\title{
Computation of hri for multiple parks
}
\description{
Function for computing the habitat replacement index for a set of protected
areas.
}
\usage{
ecohri(ecoregions, parks, indicators, ecoID = names(ecoregions)[1],
  pvals = seq(0.05, 1, 0.05), tiffdir = "tiffs", pngdir = "nopng",
  hriRes = NULL, hriRes2 = NULL, hriInRes = NULL, hriInRes2 = NULL, 
  minVar = NULL, istart = 1, tparks = 0, ecoBuffer = 1,
  wdpaid = "wdpaid", keepErrors = FALSE, pamax = 60)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{ecoregions}{\code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}
   with one or more subregions (typically ecoregions) within which to compute the hri.}
\item{parks}{\code{\link[sp:SpatialPolygons]{SpatialPolygonsDataFrame}}
   with one or more protected areas for which to compute the hri}
\item{indicators}{
  A \code{\link[raster:stack]{stack}} of the raster-package with 
  the indicators used to estimate the similarity of the habitat.}
\item{ecoID}{the column name of the id of the ecoregions. This is necessary as 
  some ecoregions can be multiple polygons}
\item{pvals}{
  The probability limit for the Mahalanobis distance used to decide if a pixel
  could be a suitable replacement for the protectedArea or not}
\item{tiffdir}{the directory where tiff-files of the similarity will be saved}
\item{pngdir}{the directory where png-files of the similarity will be saved}
\item{hriRes}{a data frame with existing results for the same parks and ecoregions, 
  to be used for a possible restart, see details below}
\item{hriRes2}{a data frame with existing results for the same parks and ecoregions, 
  to be used for a possible restart, see details below}
\item{hriInRes}{a data frame with existing results for the same parks and ecoregions, 
  to be used for a possible restart, see details below}
\item{hriInRes2}{a data frame with existing results for the same parks and ecoregions, 
  to be used for a possible restart, see details below}
\item{minVar}{a function defining the minimum variance for a set of indicators}
\item{istart}{which number in the list of ecoregions to start with - to be used
   for restarts}
\item{tparks}{total number of parks so far - to be used for restarts}
\item{ecoBuffer}{The buffer to use around ecoregion, in the metric of ecoregions}
\item{wdpaid}{the column name of the id of the parks}
\item{keepErrors}{logical; whether to keep records with errors in the resulting data.frames}
\item{pamax}{The maximum number of parks in one iteration, could be necessary to reduce for parallel computation as the maximum number of connections could be reached.}
}


\details{

The \code{ecohri} function is a method for computing the hri for one or
more parks within one or more ecoregions. The side effect is that a similarity
map is printed to a tiff-file.

}
\value{
The function will return two or four (if \code{minVar} is different from \code{NULL})
\code{data.frames} in a list. \code{hriRes} give the estimated habitat 
replaceability index (hri) from the function \code{\link{hri}} for the protected area within the
ecoregion of interest, \code{hriRes2} the similar value computed with a minimum 
variability. \code{hriInRes} and \code{hriInRes2} are the similar
hri-values inside the boundaries of the park. This is mainly for explaining purposes,
not so much for analysis.

In addition, the function will return an array with the error messages for those 
ecoregions/PAs which caused errors during computation.

}
\author{
Jon Olav Skoien
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{hri}}
}
\examples{
\dontrun{
library(eHab)
options(error = recover)
data(ecohri)
system.time(hriImage<-ecohri(ecoregions=ecoregions, parks=parks, indicators=indicators,
    tiffdir="e:/tmp/", pngdir="nopng", minVar=minVar, ecoBuffer=10, pval = 0.5))
}
}
\keyword{ spatial }
