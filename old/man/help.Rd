\name{eHab-help}
\alias{scatterplot3dm}
\alias{recenterPoints}
\alias{size}
\alias{makePolygon}
\alias{bbArea}
\title{ eHab help functions}

\description{Some useful functions in the eHab package}
\usage{
recenterPoints(spdf, action = NA)
makePolygon(xcor = c(0,1), ycor = c(0,1), coor = NULL, pstring = NA, id = "1",
              resol = NULL, data = NULL )
size(envir = .GlobalEnv)
bbArea{bb}
scatterplot3dm(x, y=NULL, z=NULL, color=par("col"), pch=NULL,
    main=NULL, sub=NULL, xlim=NULL, ylim=NULL, zlim=NULL,
    xlab=NULL, ylab=NULL, zlab=NULL, scale.y=1, angle=40,
    axis=TRUE, tick.marks=TRUE, label.tick.marks=TRUE,
    x.ticklabs=NULL, y.ticklabs=NULL, z.ticklabs=NULL,
    y.margin.add=0, grid=TRUE, box=TRUE, lab=par("lab"),
    lab.z=mean(lab[1:2]), type="p", highlight.3d=FALSE,
    mar=c(5,3,4,3)+0.1, col.axis=par("col.axis"),
    col.grid="grey", col.lab=par("col.lab"),
    cex.symbols=par("cex"), cex.axis=0.8 * par("cex.axis"),
    cex.lab=par("cex.lab"), font.axis=par("font.axis"),
    font.lab=par("font.lab"), lty.axis=par("lty"),
    lty.grid=par("lty"), lty.hide=NULL, lty.hplot=par("lty"),
    log="", ...)

}
\arguments{
  \item{spdf}{SpatialPointsDataFrame with geographical coordinates}
  \item{action}{"sub" or "add", for change from Pacific to Atlantic view, or the opposite,
              will be guessed from the data if not given.}

 \item{xcor}{limits of x-coordinates of polygon to be created or x-coordinate of centre of gravity}
  \item{ycor}{limits of y-coordinates of polygon to be created or y-coordinate of centre of gravity}
  \item{coor}{alternatively a matrix/data.frame of the limits, or xy-coordinates of centre of gravity}
  \item{pstring}{project4string}
  \item{id}{id for polygon}
  \item{resol}{size of polygon if centre of gravity is given instead of boundaries}
  \item{data}{data for the spatial polygons}
   \item{envir}{environment}
 \item{bb}{either a \code{Spatial}-object or a bounding box from a call to \code{\link[sp]{bbox}}}
 \item{x}{the coordinates of points in the plot.}
 \item{y}{the y coordinates of points in the plot, optional if \code{x} is an appropriate structure.}
 \item{z}{the z coordinates of points in the plot, optional if \code{x} is an appropriate structure.}
 \item{color}{colors of points in the plot, optional if \code{x} is an appropriate structure.
   Will be ignored if \code{highlight.3d = TRUE}.}
 \item{pch}{plotting "character", i.e. symbol to use.}
 \item{main}{an overall title for the plot.}
 \item{sub}{sub-title.}
 \item{xlim, ylim, zlim}{the x, y and z limits (min, max) of the plot. Note that setting enlarged limits
    may not work as exactly as expected (a known but unfixed bug).}
 \item{xlab, ylab, zlab}{titles for the x, y and z axis.}
 \item{scale.y}{scale of y axis related to x- and z axis.}
 \item{angle}{angle between x and y axis (Attention: result depends on scaling).}
 \item{axis}{a logical value indicating whether axes should be drawn on the plot.}
 \item{tick.marks}{a logical value indicating whether tick marks should
   be drawn on the plot (only if \code{axis = TRUE}).}
 \item{label.tick.marks}{a logical value indicating whether tick marks should be labeled on the plot
   (only if \code{axis = TRUE} and \code{tick.marks = TRUE}).}
 \item{x.ticklabs, y.ticklabs, z.ticklabs}{vector of tick mark labels.}
 \item{y.margin.add}{add additional space between tick mark labels and
   axis label of the y axis}
 \item{grid}{a logical value indicating whether a grid should be drawn on the plot.}
 \item{box}{a logical value indicating whether a box should be drawn around the plot.}
 \item{lab}{a numerical vector of the form c(x, y, len).  The values of
   x and y give the (approximate) number of tickmarks on the x and y axes.}
 \item{lab.z}{the same as \code{lab}, but for z axis.}
 \item{type}{character indicating the type of plot: "p" for points, "l"
   for lines, "h" for vertical lines to x-y-plane, etc.}
 \item{highlight.3d}{points will be drawn in different colors related to y coordinates
   (only if \code{type = "p"} or \code{type = "h"}, else \code{color} will be used).\cr
   On some devices not all colors can be displayed. In this case try the
   postscript device or use \code{highlight.3d = FALSE}.}
 \item{mar}{A numerical vector of the form c(bottom, left, top, right)
   which gives the lines of margin to be specified on the four sides of the plot.}
 \item{col.axis, col.grid, col.lab}{the color to be used for axis / grid / axis labels.}
 \item{cex.symbols, cex.axis, cex.lab}{the magnification to be used for
   point symbols, axis annotation, labels relative to the current.}
 \item{font.axis, font.lab}{the font to be used for axis annotation / labels.}
 \item{lty.axis, lty.grid}{the line type to be used for axis / grid.}
 \item{lty.hide}{line style used to plot \sQuote{non-visible} edges (defaults of the \code{lty.axis} style)}
 \item{lty.hplot}{the line type to be used for vertical segments with \code{type = "h"}.}
 \item{log}{Not yet implemented!  A character string which contains "x"
   (if the x axis is to be logarithmic), "y", "z", "xy", "xz", "yz", "xyz".}
 \item{\dots}{more graphical parameters can be given as arguments,
   \code{pch = 16} or \code{pch = 20} may be nice.}
}




\note{
scatterplot3dm is modified from the
scatterplot3d function of the package with same name, using bpy.colors instead of rgb for colors,
and fixing the y-colors for repeated calls to the function with different data.

size() gives the size of all objects in the global environment

recenterPoints can recenter a SpatialPointsDataFrame from Pacific to Atlantic centred, or the opposite

makePolygon makes a square polygon based on the x and y coordinates of the extremes
}

\author{ Jon Olav Skoien }
\keyword{ spatial }
