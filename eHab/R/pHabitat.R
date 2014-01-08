pHabitat = function(indicators = NULL, habitat = NULL, meanPoly = NULL, covPoly = NULL,
               polyData = NULL, meanCovPoly = NULL, na.rm = FALSE, indicators2 = NULL, 
               forecast = FALSE, nonNum = NULL, method = "mahalanobis", bootsample = FALSE){

                               
  preComp = (!is.null(meanPoly) & !is.null(covPoly) | !is.null(polyData) | !is.null(meanCovPoly))


# Either indicators or indicators2 has to be present
  if ((is.null(dim(indicators))) & is.null(dim(indicators2)))
     stop(paste("indicators submitted without content, class(indicators)",
                 class(indicators), class(indicators2)))
#
  if (forecast && is.null(indicators2) && dim(indicators)[2] %%2) {
    stop("indicators2 not supplied, and number of variables in indicators is odd. ", 
          "Not able to determine which part of data frame referst to current and future conditions")
  }
  full = FALSE
  if (!is.null(indicators) && fullgrid(indicators)) {
    full = TRUE
    fullgrid(indicators) = FALSE
  } 
  if (!is.null(indicators2) && fullgrid(indicators2)) {
    full = full + 1
    fullgrid(indicators2) = FALSE
  }
  if (full == 1 && !is.null(indicators) & !is.null(indicators2)) 
    stop("One set of indicators is fullgrid and the other is not")
#
# Setting NA values for extremely small values and setting forecast = TRUE
# if both indicators and indicators2 are given 
  if (!is.null(indicators) & !is.null(indicators2)) forecast = TRUE
  if (!is.null(indicators)) indicators@data[indicators@data < -1e10] = NA
  if (!is.null(indicators2)) indicators2@data[indicators2@data < -1e10] = NA
# Checking that it will be possible to compute statistics for the habitat
  if (inherits(habitat,"SpatialPolygons")) {
    if (!preComp && gridparameters(indicators)$cellsize[1]^2 > bbArea(habitat)) {
      stop("Area of interest smaller than a pixel")
    }
  }
#  
# Setting indicators2 as the grid to compute similarity for
  if (is.null(indicators2)) {
    if (preComp | !forecast) {
      indicators2 = indicators
    } else if (!preComp & forecast) {
# Forecast, but there is no precomputed data and indicators2 is missing
      ncols = dim(indicators@data)[2]
      indicators1 = indicators[,1:(ncols/2)]
      indicators2 = indicators[,(ncols/2 + 1):ncols]
    }
  }
# If indicators1 is needed and was not created above
  if (!("indicators1" %in% objects())) {
    if (!preComp) {
      indicators1 = indicators
    } else indicators1 = NULL
  }
#
# Changing habitat from 3D to 2D if necessary  
  if (!is.null(habitat) && dim(coordinates(habitat))[2] == 3 && 
          dim(coordinates(indicators))[2] == 2 &&  
          length(unique(coordinates(habitat)[,3])) == 1) {
     habitat@coords = habitat@coords[,1:2]
     habitat@bbox = habitat@bbox[1:2,]
  }
#
  rm(indicators)
  gc()
# interpolate NA values if necessary
  if (!preComp && na.rm == "int" && sum(is.na(indicators1@data)) > 0) 
    indicators1 = interpolateInd(indicators1)
  if (na.rm == "int" && sum(is.na(indicators2@data)) > 0) 
    indicators2 = interpolateInd(indicators2)
  if (na.rm == "int") na.rm = FALSE
#
  if (is.null(covPoly) | is.null(meanPoly)) {
    if (is.null(meanCovPoly)) meanCovPoly = habMeanCovPoly(indicators1, habitat, polyData = polyData, # call to habmeanCovPoly = meanCovPoly @javier
        na.rm = na.rm, method = method, bootsample = bootsample)
    meanPoly = meanCovPoly$meanPoly
    covPoly = meanCovPoly$covPoly
    inHA = meanCovPoly$inHA
    polyData = meanCovPoly$polyData
    nonNum = meanCovPoly$nonNum
    if (method == "maxent") {
      pDat = meanCovPoly$pDat
      outHA = meanCovPoly$outHA
    }
  } else  {
    if (!is.null(meanCovPoly)) stop("You cannot give both meanPoly and covPoly together with meanCovPoly")
    if (!all(!is.null(meanPoly),!is.null(covPoly))) 
       stop("You have to give either both meanPoly and covPoly, or neither")
    if (!is.null(polyData)){
      dimP = dim(polyData)[2]
      numData = which(unlist(lapply(polyData, FUN = function(x) is.numeric(x))))
      if (is.null(nonNum)) nonNum = which(!(1:dimP %in% numData))
    }
  }
  covPoly = as.matrix(covPoly)
  keep = 1:dim(indicators2)[2]
  if (!is.null(nonNum) & length(nonNum) > 0) {
    if (is.null(nnTable)) {
      nnTable = list()
      np = dim(polyData)[1]
      for (ivar in nonNum) {
        nTab  = table(polyData[,ivar])        
        nTab = nTab[nTab > 0.1*np]
        nnTable[[ivar]] = nTab
      }
    } 
    keep = keep[-nonNum]
  }    
  if (!"dimP" %in% objects()) dimP = length(meanPoly)
  noVar = which(diag(covPoly) == 0)
  pCov = covPoly
  dimnames(pCov) = list(1:dimP, 1:dimP)
  if (length(noVar) > 0) pCov = pCov[-noVar, -noVar] 
  
  dcov = abs(pCov/pCov[,1])
  idup = which(duplicated(signif(dcov,4)))
  dupVar = as.integer(dimnames(pCov)[[1]][idup])
  pKeep = 1:dimP
  if (length(noVar) > 0 | length(dupVar) > 0) {
    if (length(noVar) > 0) 
#        warning(paste("Removing variable(s) without variance inside area",
#                 names(noVar)))
    if (length(dupVar) > 0) {
      dnames = names(polyData)[dupVar]
      duplic = duplicated(abs(covPoly)) | duplicated(abs(covPoly),fromLast = TRUE)
      onames = names(polyData)[duplic][!(names(polyData)[duplic] %in% dnames)]
#      warning(paste("Removing variable(s) 1:1 correlated with other variable(s)",
#          dnames, "correlated with", onames ))
    }
    remVar = c(noVar,dupVar)
    keep = keep[-remVar]
    if (length(keep) == 0) 
      stop("All indicators are flagged as constant within habitat. Not possible to compute covariances.")
    if (length(noVar) > 0) noVarMean = meanPoly[noVar]
    meanPoly = meanPoly[-remVar]
    covPoly = covPoly[-remVar,-remVar]
    pKeep = pKeep[-remVar]
  }
  if (method == "mahalanobis") {
    res = try(mahalanobis(indicators2@data[,keep, drop = FALSE],meanPoly,covPoly ))
    if (!is(res, "try-error")) {
      indicators2$mDist = res
      indicators2$pHab = 1-pchisq(indicators2$mDist,length(keep))
    } 
  } else { 
    ind2 = stack(indicators2)
    if (method == "maxent") {
      object = maxent(polyData[,pKeep], pDat)
    } else if (method == "bioclim") {
      object = bioclim(polyData[,pKeep])
    } else if (method == "domain") {
      object = domain(polyData[,pKeep])
    }
    res = predict(object, ind2)
    fullgrid(indicators2) = TRUE
    if (!is(res,"try-error")) indicators2$pHab = extract(res, extent(res))
    fullgrid(indicators2) = FALSE    
  }
  if (!is(res,"try-error")) {
    if (length(noVar) > 0) 
      for (inam in 1:length(noVar)) 
         indicators2$pHab[!(indicators2@data[,noVar[inam]] == noVarMean[inam])] <- 0
    if (!is.null(nonNum)) {
      for (ivar in nonNum) {
        indicators2$pHab[!(indicators2@data[,ivar] %in% names(nnTable[[ivar]]))] = 0
      }
    }
    gridded(indicators2) = TRUE
    if (full) fullgrid(indicators2) = TRUE
    if (!is.null(meanCovPoly)) attr(indicators2, "meanCovPoly") = meanCovPoly
    attr(indicators2,"cov") = covPoly
    attr(indicators2,"mean") = meanPoly
    if (exists("polyData")) {
      attr(indicators2,"polyData") = polyData
      attr(indicators2, "hist") = lapply(polyData[,pKeep], 
                   FUN = function(x) hist(x, plot = FALSE))
    }
    if (!is.null(habitat)) {
      polyData2 = try(habMeanCovPoly(indicators2, habitat, # call to habmeanCovPoly = polyData2 @javier
        na.rm = na.rm, method = method)$polyData)
      attr(indicators2, "polyDataForecast") = polyData2
      if (!is(polyData2, "try-error")) attr(indicators2, "histForecast") = 
                   lapply(polyData2[,pKeep], 
                   FUN = function(x) hist(x, plot = FALSE))
    }
    indicators2
  } else {
    warning(paste("not able to compute Mahalanobis distances \n",
       paste(dim(covPoly), collapse = " "), length(meanPoly), "\n", 
       paste(meanPoly, collapse = " ")))
    attr(res,"cov") = covPoly
    attr(res,"mean") = meanPoly
    if (!is.null(meanCovPoly)) attr(res, "meanCovPoly") = meanCovPoly
    if (exists("polyData") && !is.null(polyData)) {
      attr(res,"polyData") = polyData
      attr(res, "hist") = lapply(polyData[,!(1:dim(polyData)[2] %in% nonNum)], 
                   FUN = function(x) hist(x, plot = FALSE))
    }
    res
  }
}

   
   

interpolateInd = function(indicators) {
  for (i in 1:dim(indicators@data)[2]) {
    if (sum(is.na(indicators@data[,i])) > 0) {
      isna = which(is.na(indicators@data[,i]))
      obs = try(SpatialPointsDataFrame(coordinates(indicators)[-isna,],
                  data = data.frame(obs = indicators@data[-isna,i])))
      if (is(obs,"try-error")) {
        stop(paste("Cannot identify indicators with data for interpolation.\n",
            "dim indicators:",paste(dim(indicators), collapse = " "), "\n",
            "length(isna):", length(isna)))
      }
      ploc = SpatialPoints(coordinates(indicators)[isna,])
      require(automap)
      if (dim(obs@data)[1] > 1000) {
        data_variogram = obs[sample(c(1:dim(obs@data)[1]),1000),]
      } else data_variogram = obs
      plocD = autoKrige(obs~1, obs, data_variogram = data_variogram, 
                        ploc, nmax = 8, maxdist = 0.01)$krige_output
      indicators@data[isna,i] = plocD$var1.pred
    }
  }
  indicators
}


bbArea = function (bb) {
    if (inherits(bb, "Spatial")) bb = bbox(bb)
    xd = bb[[3]] - bb[[1]]
    yd = bb[[4]] - bb[[2]]
    abs(xd) * abs(yd)
}

