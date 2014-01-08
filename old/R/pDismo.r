pDismo = function(habitat, indicators, method = "mahalanobis", ...) {

#, inHA = NULL, meanPoly = NULL, covPoly = NULL,
#               polyData = NULL, na.rm = FALSE, indicators2 = NULL, forecast = FALSE, nonNum = NULL,
#               method = "mahalanobis"){

indicators = stack(indicators)
habitat = SpatialPoints(habitat)

if (method == "mahalanobis") {
  object = mahal(indicators, habitat)
} else if (method == "maxent") {
  object = maxent(indicators, habitat)
} else if (method == "bioclim") {
  v = extract(indicators, habitat)
  object <- bioclim(v)
} else if (method == "domain") {
  object = domain(indicators, habitat)
}

predict(object, indicators)
}

 
 
