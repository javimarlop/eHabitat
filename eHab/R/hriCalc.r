hriCalc = function(object, habitat, pvals, populationDensity) {
  sims = grep("sim", names(object))
  if (length(sims) > 0) sims = names(object)[sims]
  lnames = c("pHab", sims)
  if (!inherits(object, "Raster")) {
    if (fullgrid(object)) fullgrid(object) = FALSE
    object = brick(object[, lnames]) # @javier 
# The raster:::subset construction should not be necessary, but sometimes
# the method did not dispatch correctly.
  } else object = raster:::subset(object, lnames) # @javier 
  inHA = try(mask(object, habitat)) # @javier 
  if (is(inHA, "try-error")) {
    inHA= NULL
    outHA = object
  } else outHA = mask(object, habitat, inverse = TRUE) # @javier 
# Necessary until update of raster package
  if (!is.null(inHA) && !is(inHA, "RasterBrick")) inHA = brick(inHA) # @javier 
  if (!is(outHA, "RasterBrick")) outHA = brick(outHA) # @javier 

# The method below is more correct for computing the excact number of pixels
# but it does not subtract NA pixels
#  paarea = sum(sapply(pa@polygons, FUN = function(x) x@area))
#  ores = xres(object) * yres(object)
#  inPixels = paarea/ores
  inPixels = ifelse(is.null(inHA), 0, sum(!is.na(getValues(raster(inHA,1)))))  
  np = length(pvals)+1  
  ret = data.frame(pvals = pvals, matrix(NA, ncol = 2*length(sims), nrow = length(pvals)))
  no = nlayers(object)
  for (i in 1:nlayers(object)) {
    if (!is.null(inHA)) {
      inStats = hist(raster(inHA,i), breaks = c(0,pvals, 1.1), plot = FALSE)$counts
      inStats = cumsum(inStats[np:1])[np:1]
    } else inStats = rep(0, length(pvals)+1)
    ret[,(i*2)] = inStats[2:length(inStats)]/inPixels
    outStats = hist(raster(outHA,i), breaks = c(0,pvals, 1.1), maxpixels = Inf, plot = FALSE)$counts
    outStats = cumsum(outStats[np:1])[np:1]
    ret[,i*2+1] = outStats[2:length(outStats)]/inPixels
  }
#  if (nlayers == 1) {
#    names(ret) = c("inHRI", "outHA")
  names(ret)[2:(nlayers(object)*2+1)] = paste(rep(c("in", "out"), length(lnames)), 
                 rep(lnames,each = 2), sep = "_")
  rnames = names(ret)
  if (length(grep("sim",rnames)) > 0) {
    snames = rnames[grep("sim", rnames)]
    ret$inHAm = rowMeans(ret[,snames[grep("in",snames)]])
    ret$inHAstd = apply(ret[,snames[grep("in",snames)]], MARGIN = 1, FUN = function(x) sd(x))
    ret$outHAm = rowMeans(ret[,snames[grep("out",snames)]])
    ret$outHAstd = apply(ret[,snames[grep("out",snames)]], MARGIN = 1, FUN = function(x) sd(x))
  }
  ret
} 
 
  

hriCalcOld = function(object, pa, pvals, populationDensity) {
  if (fullgrid(object)) fullgrid(object) = FALSE
  if (!inherits(object, "Raster")) object = raster(object[, "pHab"])
  inHA = try(mask(object, pa))
  if (is(inHA, "try-error")) {
    inHA= NULL
    outHA = object
  } else outHA = mask(object, pa, inverse = TRUE)
#  ores = res(object)
#  inPixels = polyArea(pa)/ores[1]/ores[2]
  inPixels = ifelse(is.null(inHA), 0, sum(!is.na(getValues(inHA))))  
  np = length(pvals)+1
  ret = data.frame(pvals = pvals, inHA = NA, outHA = NA )
  inStats = hist(inHA, breaks = c(0,pvals, 1.1), plot = FALSE)$counts
  inStats = cumsum(inStats[np:1])[np:1]
  ret[,2] = inStats[2:length(inStats)]/inPixels
  outStats = hist(outHA, breaks = c(0,pvals, 1.1), maxpixels = Inf, plot = FALSE)$counts
  outStats = cumsum(outStats[np:1])[np:1]
  ret[,3] = outStats[2:length(outStats)]/inPixels
  ret
}
  
 
maskEco = function (x, mask, ...) 
{
    .local <- function (x, mask, filename = "", inverse = FALSE, 
        ...) 
    {
        if (inverse) {
            mask <- rasterize(mask, x, -1)
            mask(x, mask, filename = filename, inverse = TRUE, 
                ...)
        }
        else {
            if (nlayers(x) > 1) {
                mask <- rasterize(mask, x, -1)
                mask(x, mask, filename = filename, ...)
            }
            else {
                rasterize(mask, x, filename = filename, mask = TRUE, 
                  ...)
            }
        }
    }
    .local(x, mask, ...)
}

