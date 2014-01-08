\name{findPatches}
\alias{findPatches}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Identify continuous patches }
\description{
Function that identifies continuous patches of possible replacement areas for 
a habitat}
\usage{
findPatches(hrep, plot = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{hrep}{
     an \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} object from the function \code{\link{hri}} with pixels that could
     possible replace the pixels in the original habitat
}
  \item{plot}{logical; defines whether the function should plot the patches}
}


\value{
  The function returns the submitted \code{\link[sp:SpatialPixelsDataFrame]{SpatialGridDataFrame}} 
  including a new column with IDs of patches
}
\author{
Jon Olav Skoien
}

\keyword{spatial}
