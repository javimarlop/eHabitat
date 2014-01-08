# There are two potential calls to this function in the pHabitat script: L. 78 and 183 @javier
# There is one potential call to this function in the ecohri script: L. 152 @javier


habMeanCovPoly = function(indicators, habitat, inHA = NULL, polyData = NULL, 
           method = "mahalanobis", mAPoints = 1000, na.rm = FALSE, bootsample = FALSE) {
  weight = 1
  if (is.null(inHA) & is.null(polyData)) {
    if (inherits(habitat, "SpatialPoints")) {
      inHA = over(habitat, geometry(indicators))
      habitat$inHA = inHA
      inHA = inHA[!is.na(inHA)]
      if ("weight" %in% names(habitat)) {
        indicators$weight = 0
        indicators$weight[inHA] = habitat$weight[!is.na(habitat$inHA)]
      }
    } else if (inherits(habitat, "list") && inherits(habitat[[1]], "SpatialPoints")) {
      weight = c(1:dim(indicators)[1])
      inHA = NULL
      for (ilist in 1:length(habitat)) {
        linHA = overlay(indicators,habitat[[ilist]])
        linHA = linHA[!is.na(linHA)]
        weight[linHA] = weight[linHA] +1
        inHA = c(inHA, linHA)
      }
      inHA = unique(inHA)
    } else if (inherits(indicators, "Spatial") && inherits(habitat,"SpatialPolygons")) {
#
#
      rr <- over(indicators, SpatialPolygons(habitat@polygons, pO = habitat@plotOrder,
               CRS(proj4string(habitat))), 
            returnList = FALSE)
      inHA = which(!is.na(rr))
      if (length(habitat) > 1) {
#        r2 = sp:::.invert(rr, length(habitat), length(indicators))
#        indicators$weight = unlist(lapply(r2, FUN = function(x) length(x))) 
        rr[is.na(rr)] = 0
        indicators$weight = rr
      }
    } else if (!is.null(indicators) && !is.null(habitat)) {
      stop(paste("Overlay for class ", class(habitat), "not yet defined"))
    }
  }
  if (is.null(polyData)) {
    pl = length(inHA)
    if (pl <= 2) {
      stop("Area of interest covers 2 or fewer pixels, not possible to calculate statistics")
    }

    if ("weight" %in% names(indicators)) {
      weight = indicators$weight[inHA]
      indicators = indicators[,-which(names(indicators) == "weight")]  
    }
    polyData = indicators@data[inHA,, drop = FALSE]
  } else weight = 1

  naPoly = (rowSums(is.na(polyData)) > 0)
  if (sum(naPoly) >=1) {
    polyData = polyData[!naPoly,,drop = FALSE]
    if (length(weight) > 1) weight = weight[!naPoly]
  }
  if (bootsample) {
    bootsample = min(dim(polyData)[1], 1000)
    polyData = polyData[sample(1:dim(polyData)[1],bootsample, replace = TRUE),, drop = FALSE]
  }
  dimP = dim(polyData)[2]
  numData = which(unlist(lapply(as.data.frame(polyData), FUN = function(x) is.numeric(x))))
  nonNum = which(!(1:dimP %in% numData))
  if (length(nonNum) == 0) nonNum = NULL
  if (all(weight == 1)) {
    meanPoly = colMeans(as.matrix(polyData[,numData]), na.rm = na.rm)
    covPoly = var(polyData[,numData],,na.rm = na.rm)
  } else {
#    txt = paste("length habitat:", dim(habitat)[1], "length polyData", dim(polyData)[1],
#                "length naPoly", length(naPoly), naPoly)
#    stop(txt)
    wcov = cov.wt(polyData[,numData],weight)
    meanPoly = wcov$center
    covPoly = wcov$cov      
  }
  if (method == "maxent") {
    outHA = sample((1:dim(indicators)[1])[-inHA],min(dim(indicators)[1]-length(inHA),max(length(inHA), mAPoints)))
    outData = indicators@data[outHA,]
    polyData = rbind(polyData, outData)
    pDat = c(rep(1,length(inHA)), rep(0, length(outHA)))
  } else {
    outHA = NULL
    pDat = NULL
  }
  list(meanPoly = meanPoly, covPoly = covPoly, inHA = inHA, polyData = polyData,
       nonNum = nonNum, outHA = outHA, pDat = pDat)



}
   

#if (FALSE) {  
#library(eHab)
#load("e:/temp/debug2.rda")
#options(error = recover)
##debug(mhri)
#hriim = mhri(imageSpatial, spatialGeometry, nboot = 5, nsim = 3)
#
#nn = names(ind1)[3]
#spplot(ind1, nn,col.regions = bpy.colors(), 
# panel = function(x,y, ...){
#    panel.gridplot(x,y, ...)
#    sp.points(habitat, col = "red", lwd = .5)})
#    
#}




identicalCRS = function(x,y) {
  if ("rgdal" %in% .packages()) {
    identical(CRS(proj4string(x)), CRS(proj4string(x)))
  } else {
    x = gsub("^ .", "", proj4string(x))
    y = gsub("^ .", "", proj4string(y))
    identical(x,y)
  }
} 
