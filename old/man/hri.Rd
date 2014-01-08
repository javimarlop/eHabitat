\name{hri}
\alias{hri}
\alias{mhri}
\alias{whri}
\title{
Habitat replacement index
}
\description{
Functions that calculates similarities between habitat of interest and surroundings
and/or future.
}
\usage{
hri(indicators = NULL, habitat = NULL, geometry = NULL, populationDensity = NULL, 
  pvals = seq(0,0.95,0.05), forecast = FALSE, hriCalc = TRUE, ...)
mhri(indicators, habitat, indicators2 = NULL, stdind = NULL,
              forecast = FALSE, nsim = 1, nboot = 0, range = NULL, 
              exc = 0.5, nugget = NULL,  lmc, debug.level = -1, 
              method = "mahalanobis",
              pvals = seq(0,0.95,0.05), nclus = 1, ...) 
whri(habitat, geometry, indicators, forecast = FALSE, indicators2 = NULL, ...)



}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{indicators}{
  A \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} with 
  the indicators used to estimate the similarity of the habitat. 
  The \code{data.frame} can include local variances when calling \code{mhri} if
  \code{forecast = FALSE}. }
\item{habitat}{
  A \code{\link[sp]{SpatialPolygons}}* object with the boundaries of the 
  protected area or habitat of interest, or a \code{\link[sp]{SpatialPoints}}
  object with precence locations}
\item{geometry}{Either a \code{\link[sp]{SpatialPolygons}}* object with the 
     ranges of all species of interest from a protected area (PA) or important
     bird area (IBA), or a \code{\link[sp:SpatialPoints]{SpatialPointsDataFrame}}
     with observations of all the species of interest from the PA or IBA}     
\item{populationDensity}{
  A \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} or a 
  \code{\link{data.frame}} with the population density of the same 
  locations as used for the indicators. If not available, the population density
  is set equal to 1 within the protectedArea and 0 outside.}
\item{pvals}{
  The probability limit for the Mahalanobis distance used to decide if a pixel
  could be a suitable replacement for the protectedArea or not}
\item{indicators2}{A second \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridPixelsFrame}} 
  with forecasted indicators. The \code{data.frame} can also contain local variances 
  of the forecast data when calling \code{mhri} }
\item{forecast}{Logical; if TRUE and neither \code{indicators2} nor 
    \code{covPoly} and \code{meanPoly} are not submitted,
    the data of the indicators will be split in two columnwise, assuming that the
    first part of the data frame represents the current state and the second part 
    the forecast. If covPoly and meanPoly are submitted. If \code{forecast} is TRUE,
    also the area within the habitat will be included for computation of hri.}
\item{stdind}{A \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridPixelsFrame}} with 
    spatial distribution of the standard deviations of the observations.}
\item{nsim}{Number of simulations to use when doing Monte Carlo simulations
    of uncertain input. Function will call \code{hri} directly when \code{nsim = 1}}
\item{nboot}{The number of bootstrap samples for estimates of the model uncertainty}
\item{range}{Correlation length to be used in conditional simulations, 
    \code{range = NULL} gives automatic fitting of range, see details below.}
\item{exc}{The similarity level for computation of exceedance probabilities}
\item{nugget}{Nugget effect to be used in conditional simulations}
\item{lmc}{Linear model of coregionalization, to be used in conditional simulations, 
     see also \code{\link[gstat]{fit.lmc}} for more details. The model should
     have \code{sill = 1}, as the local variance is given separately.}
\item{debug.level}{Passed to \code{\link[gstat]{predict.gstat}} when creating
     realizations of uncertain input}
\item{method}{Which method to use for computing similarites. The possibilities are
      at the moment "mahalanobis", "maxent", "bioclim" and "domain". The three last
      methods are mainly available through the \code{\link{dismo}}-package.}
\item{nclus}{The number of clusters to use in case of spatial simulations with 
     paralell processing of the simulations}
\item{hriCalc}{logical; whether the hri should be computed, mainly used by higher level functions}
\item{...}{Variables to be passed to \code{hri} from \code{mhri} or \code{whri}, 
     or from \code{whri} to \code{mhri}}
}
\details{


\code{hri} is called by the two other functions
after some preprocessing of their input. The base of all these functions are
\code{\link{pHabitat}}, which computes similarities
to a region of interest from the surroundings and/or future, based on the 
Mahalanobis distance to the mean and variance of a set of indicators at the 
region of interest. See the reference below for more details.

If the function is used for forecasting, the indicators for the future can either 
be given in variable \code{indicators2}, or as additional columns of \code{indicators}
with \code{forecast = TRUE}.

\code{mhri} is a function for estimating habitat similarity, assuming uncertainty 
of the indicators. The correlation of the uncertainty can either be supplied as a linear model 
of coregionalization (see also \code{\link[gstat]{fit.lmc}}), or through 
nugget and range. If none are supplied, the function will estimate range from
the average of ranges fitted to the sample variograms of the individual 
indicators with \code{\link[automap]{autofitVariogram}}. 

The local variability can be given as standard deviation, using the same
column names with an additional \code{u} first, so that original column names
of \code{bio, epratio and prec} would give standard deviation column names
of \code{ubio, uepratio and uprec}. These columns can be added to 
\code{indicators2}, or to \code{indicators} if \code{indicators2} is not given
and \code{forecast = FALSE}.

The function will for simplicity assume no uncertainty for the indicators 
used to estimate the mean and covariance at the region of interest. Instead,
the Linear model of coregionalization is used to create sets of uncoditional simulations
of errors, that are multiplied with the local standard deviations and added
to the original set of indicators to be used for estimation of similarity. 

\code{whri} is called when the mean and variance to be used for computation
of habitat similarity should be based on a weighted set of indicators. The variable
\code{geometry} will usually give a set of range maps or occurence observations
of trigger species 
observed in a protected area or important bird area. Mean and variance for the
region of interest is then computed as a weighted mean of the values of the indicators 
for these ranges/occurence locations, rather than for the area of the PA/IBA
itself.


}
\value{
\code{hri} calls \code{\link{pHabitat}}, which adds the columns 
\code{mDist} (Mahalanobis distance) and \code{pHab} (similarity between the 
indicators of pixels and their values at the locations of interest) to the 
indicators object. In addition, there is a column \code{replace} with values 0/1, with 1 if the
pixel has a value larger than \code{pval}.
 
The object will also have three attributes:

\code{attr(spdf, "hri"} A data frame with the total size of the area with 
   \code{pHab > pval}, compared to the original protected area, with one
   row for each value of \code{pval}. The hri
   is split between the hri inside and outside the park boundaries. 
   For \code{mhri}, this attribute gives the hri for each individual simulation, 
   in addition to the mean over simulations and the standard deviation.

\code{attr(spdf,"mean") = meanPoly} The mean values of indicators used in 
computation of similarity.

\code{attr(spdf,"cov") = covPoly} The covariance matrix of indicators used
in computation of similarity.


\code{mhri} will also return the results of the individual realisations and some
summary statistics.

\code{whri} wirll return the same results as \code{hri} if called with
\code{nsim = 1} and as \code{whri} if called with \code{nsim > 1}.


}
\references{
Clark, J.D., Dunn, J.E. and Smith, K.G. (1993) A multivariate model of female 
black bear habitat use for a geographic information system. Journal of Wildlife 
Management, 57, 519-526. 
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
hrep = hri(indicators,protectedArea, pval = 0.5)
spplot(hrep,"pHab",col.regions = bpy.colors(), 
      panel = function(x,y, ...){
        panel.gridplot(x,y, ...)
        sp.polygons(protectedArea, col=2,fill=0,lw = 2)
      }
)

hrep = mhri(indicators,protectedArea, pval = 0.5, nsim = 4)

}
\keyword{ spatial }
