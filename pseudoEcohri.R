

ecohri = function(ecoregions = NULL, parks, indicators, ecoID = names(ecoregions)[1],
                  pvals = seq(0.05, 1, 0.05), tiffdir = "tiffs", pngdir = "nopng",
                  hriRes = NULL, hriRes2 = NULL, hriInRes = NULL, hriInRes2 = NULL, 
                  minVar = NULL, istart = 1, tparks = 0, ecoBuffer = 1,
                  wdpaid = "wdpaid", keepErrors = FALSE, pamax = 60) {
  
  parks$ecoID = find EcoID
  create empty HRI matrix if not existing
  ecoregs = unique(parks$ecoID) , and removing those not in ecoregions
  for (iie in istart:length(ecoregs)) {
    ecorr = ecoregions[iie,] + buffer
    indicators2 = mask(crop(indicators, ecorr))
    cpm = canProcessInMemory(indicators2, n = 15) # @javier
    pai = how many groups of 60 parks in ecoregion
    for (ipai in 1:pai) {
      if (!cpm) {
        Create temp files for similarity rasters
      }
      for (it in 1:numberOfTiles) {
        icells = cellFromRowColCombine(indicators2,tr$row[it]: (tr$row[it]+tr$nrows[it]-1), 1:cols)
        ind <- getValuesBlock(indicators2, row=tr$row[it], nrows=tr$nrows[it])
        indispdf = SpatialPointsDataFrame(xyFromCell(indicators2,icells), data = data.frame(ind))
        for (jpark in 1:min(length(parksInEcoregion), 60)){
          if (it == 1) {
            pdat = extract(indicators2, pa, df = TRUE, small = TRUE)[,-1]
            hmc = try(habMeanCovPoly(polyData = pdat)) # Find mean and covariance
            habMeanCov[[jpark]] = hmc
          }
          hriLoc = hri(indispdf, habMeanCov[[jpark]])
          if (!is.null(minVar)) hriLoc2 = hri(indispdf, minVar(habMeanCov[[jpark]]))
        }
        if (cpm) {
          writePngs(hriLoc and hriLoc2)
          writeRasters(hriLoc and hriLoc2)
          hriRes[jpark,] = hriCalc(hriLoc and hriLoc2)
        } else {
          writeValues(tempfiles, hriLoc and hriLoc2)
        }
      }
      if (!cpm) {
        for (jpark in 1:min(length(parksInEcoregion), 60)){
          writeRasters(tempfiles[[jpark]])
          hriRes[jpark,] = hriCalc(tempfiles[[jpark]])
        }
      }
    }
  }
  return(list(hriRes, hriInRes, hriRes2, hriInRes2, errors))
}
