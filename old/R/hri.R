#hri = function(indicators = NULL, habitat = NULL, populationDensity, pval = 0.5, inHA = NULL, 
#               meanPoly = NULL, covPoly = NULL, polyData = NULL, na.rm = FALSE, forecast = FALSE,
#               indicators2 = NULL) {
#  pHab = pHabitat(indicators, indicators2 = indicators2, habitat, inHA, meanPoly, 
#        covPoly, polyData, na.rm = na.rm, forecast = forecast)



hri = function(indicators = NULL, habitat = NULL, geometry = NULL, populationDensity = NULL, 
               pvals = seq(0,0.95,0.05),  forecast = FALSE, hriCalc = TRUE, ...) {
  if (!is.null(geometry)) {
    pa = habitat
    habitat = geometry
  }
  if (!is.null(habitat) && dim(coordinates(habitat))[2] == 3 && 
          dim(coordinates(indicators))[2] == 2 &&  
          length(unique(coordinates(habitat)[,3])) == 1) {
     habitat@coords = habitat@coords[,1:2, drop = FALSE]
     habitat@bbox = habitat@bbox[1:2,]
  }

  pHab <- pHabitat(indicators = indicators, habitat = habitat, forecast = forecast, ...) # changed by @javier (26.07.13)
  if (!is.null(geometry)) habitat = pa
  if (!is(pHab,"try-error") && inherits(habitat, "SpatialPolygons")) {
    polyData = attr(pHab, "polyData")
    if (hriCalc && !is.null(habitat)) 
         attr(pHab, "hri") = hriCalc(pHab, habitat, pvals, populationDensity)
  } else {attr(pHab,"hri") = NA}
  pHab
}





