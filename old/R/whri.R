


whri = function(habitat, geometry = NULL, indicators, forecast = FALSE, indicators2 = NULL, 
    ...) {
  
  indispdf = indicators
  if (inherits(geometry, "SpatialPolygons")) {
    w2 = over(indicators, geometry, returnList = TRUE)
    indispdf$weight = sapply(w2, length)
  } else if (inherits(geometry, "SpatialPoints")) {
    print("points not implemented yet")
  } else if (is.null(geometry)) {
    return(mhri(indicators, habitat, meanPoly = meanPoly, covPoly = covPoly, forecast, ...))  
  } else stop(paste("weighted eHab not possible, class of geometry:", class(geometry)))
  
  ndim = dim(indicators)[2]
  nvar = ifelse(forecast && is.null(indicators2), ndim/2,ndim)
  meanPoly = sapply(indispdf@data[,1:nvar], weight = indispdf$weight, 
            FUN = function(x, weight) weighted.mean(x, weight))
  covPoly = cov.wt(indispdf@data[,1:nvar], indispdf$weight)$cov
  mhri(indicators, habitat, meanPoly = meanPoly, covPoly = covPoly, forecast = forecast, ...)   
}





#whri2(pa, mammr, indispdf1)
