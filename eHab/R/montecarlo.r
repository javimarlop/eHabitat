
makeUnc = function(indispdf){
  for (ii in 1:dim(indispdf)[2]) {
    indispdf@data[,ii] = min(indispdf@data[indispdf@data[,ii]>0,ii], na.rm = TRUE)/10+
                   0.05*indispdf@data[,ii]
  }
  indispdf
}



mvarSim = function(lmc, spdf, stdind, nsim = 1, nmax = Inf, maxdist = Inf, debug.level = debug.level, nclus = 1) {
  if (length(spdf) < 1e5) {
    lmcm = lmc$model
    spdf@data[spdf@data < -1e10] = NA
    vars = names(lmcm)[-grep("\\.", names(lmcm))]
    nvar = length(vars)
#    tvar = nvar*2
    range = lmcm[[1]]$range
#    dia = sqrt(bbArea(bbox(spdf)))
#    ninit = as.integer(dia/range*4)
    gg.dummy = NULL
    for (i in 1:nvar) {
      gg.dummy <- gstat(gg.dummy, vars[i], formula = as.formula(paste(vars[i], "~1")), 
         locations = ~x+y, dummy = TRUE,  beta = 0, nmax = nmax, maxdist = maxdist)
    }
    gg.dummy$model = lmcm
    for (itest in 1:5) {
      yy <- try(predict.gstat(gg.dummy, newdata = spdf, nsim = nsim, debug.level = -1))
        if (!is(yy, "try-error")) break
    }
  } else yy = mvSim(lmc, newdata = geometry(spdf), nmax = 10, nsim = 1, debug.level = -1)
 
  spdf@data = spdf@data + yy@data * stdind@data
  spdf
}



mhri = function(indicators, habitat, indicators2 = NULL, stdind = NULL,
              forecast = FALSE, nsim = 1, nboot = 0, range = NULL, 
              exc = 0.5, nugget = NULL,  lmc, debug.level = -1, 
              method = "mahalanobis", 
              pvals = seq(0,0.95,0.05), nclus = 1, ...) {

if (!identical(proj4string(indicators),proj4string(habitat)))
   stop(paste("mhri: projections of indicators and habitat are not equal:",
               proj4string(indicators),proj4string(habitat)))
if (!is.null(indicators2) && !identical(proj4string(indicators2),proj4string(habitat)))
   stop(paste("mhri: projections of indicators2 is different:",
               proj4string(indicators),proj4string(indicators2)))
if (!is.null(stdind) && !identical(proj4string(stdind),proj4string(habitat)))
   stop(paste("mhri: projections of stdind is different:",
               proj4string(indicators),proj4string(stdind)))
                    
proj4string(indicators) = CRS(proj4string(indicators))
proj4string(habitat) = CRS(proj4string(indicators))
if (!is.null(indicators2)) proj4string(indicators2) = CRS(proj4string(indicators))
if (!is.null(stdind)) proj4string(stdind) = CRS(proj4string(indicators))

  if (!is.null(habitat) && dim(coordinates(habitat))[2] == 3 && 
          dim(coordinates(indicators))[2] == 2 &&  
          length(unique(coordinates(habitat)[,3])) == 1) {
     habitat@coords = habitat@coords[,1:2]
     habitat@bbox = habitat@bbox[1:2,]
  }

indicators@data[indicators@data < -1e10] = NA
if (!is.null(indicators2)) indicators2@data[indicators2@data < -1e10] = NA
if (nsim == 1) {
  return(hri(indicators, habitat = habitat, indicators2 = indicators2, forecast = forecast, 
             pvals = pvals, hriCalc = TRUE, method = method,  ...))
} else {
  if (forecast) {
    if (is.null(indicators2)) {
      nvar = dim(indicators)[2]/2
      indicators2 = indicators
      indicators2@data = indicators2@data[,(nvar+1):(nvar*2)]
      indicators@data = indicators@data[,1:nvar]
    } else nvar = dim(indicators2)[2]
  } else {
    indicators2 = indicators
    nvar1 = dim(indicators)[2]
    nv1 = nvar1/2
    if (!(nvar1 %% 2) && sum(names(indicators)[1:nv1] == 
            sub("u","",names(indicators)[(nv1+1):nvar1])) == nv1) {
#     Removing uncertainty from explanatory indicators, in case indicators 2 was created directly 
#     from these
      indicators@data = indicators@data[,1:nv1]
    }
    nvar2 = dim(indicators2)[2]
    nv2 = nvar2/2
    if (!(nvar2 %% 2) && sum(names(indicators2)[1:nv2] == 
            sub("u","",names(indicators2)[(nv2+1):nvar2])) == nv2) {
      stdind = indicators2[,(nv2+1):nvar2]
      indicators2@data = indicators2@data[,1:nv2]
    }
  }
  names(indicators) = gsub("\\.", "", names(indicators))
  names(indicators2) = gsub("\\.", "", names(indicators2))

  if (fullgrid(indicators)) {
    full = TRUE
    fullgrid(indicators) = FALSE
    fullgrid(indicators2) = FALSE
    sftry = try(fullgrid(stdind), silent = TRUE)
    if (!is(sftry, "try-error") && sftry) fullgrid(stdind) = FALSE
  } else full = FALSE   
  hrepp = hri(indicators, habitat = habitat, indicators2 = indicators2, forecast = forecast, 
               pvals = pvals, hriCalc = FALSE, method = method,  ...)
  if (!is(hrepp,"try-error")) {
    gridded(indicators2) = FALSE
    if (is.null(stdind)) {
#      nvar not divisible with 2 or names of indicators dont match
      stdind = makeUnc(indicators2)
    }      
    gridded(indicators2) = FALSE
    btrans = FALSE
    if (length(grep("4326",proj4string(indicators2))) == 1) {
      indicators2 = spTransform(indicators2, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"))
      stdind = spTransform(stdind, CRS("+proj=moll +lon_0=0 +x_0=0 +y_0=0 +ellps=WGS84"))
      btrans = TRUE
    }
    dia = sqrt(bbArea(bbox(indicators2)))
    maxdist = ifelse(is.null(range), dia, range)
    if (missing(lmc) || is.null(lmc)) {
      lmc = makeLmc(indicators2, model = "Exp", range = range, nugget = nugget)   
      if (is.null(range)) range = lmc[[1]]$range[2]
    }
    jsim = 0
    for (isim in 1:nsim) {
      ers = mvarSim(lmc, indicators2, stdind, maxdist = maxdist, nmax = 10, 
          debug.level = debug.level, nclus = nclus,  ...)     
      for (ivar in 1:dim(ers)[2]) {
        ers@data[,ivar][ers@data[,ivar] < min(indicators2@data[,ivar])] = min(indicators2@data[,ivar])
        ers@data[,ivar][ers@data[,ivar] > max(indicators2@data[,ivar])] = max(indicators2@data[,ivar])
      }
#      ers$epratio[ers$epratio > 100] = 100
#      ers$epratio[ers$epratio < 0.001] = min(indicators2$epratio)
#      ers$bio[ers$bio < 0.1] = 0.1
#      ers$prec[ers$prec < 0.1] = 0.1     
      if (btrans) ers = spTransform(ers, CRS(proj4string(indicators)))
      gridded(ers) = TRUE

      for (iboot in 1:max(1,nboot)) {
        if (forecast) {
          hrepl = hri(indicators, habitat = habitat, indicators2 = ers, forecast = TRUE,
                     pvals = pvals,  hriCalc = FALSE, bootsample = nboot, method = method, ...)
        } else hrepl = hri(indicators = ers, habitat = habitat, forecast = FALSE, 
                     pvals = pvals, hriCalc = FALSE, bootsample = nboot, method = method,  ...)
        if (!is(hrepl,"try-error")) {
          jsim = jsim + 1
          if (jsim == 1) {
             hrepp$sq = 0
             hrepp$mean = 0
           } 
           hrepp$sq = hrepp$sq + hrepl$pHab^2
           hrepp$mean = hrepp$mean + hrepl$pHab
           hrepp@data = cbind(hrepp@data,hrepl$pHab)
           names(hrepp)[dim(hrepp@data)[2]] = paste("sim",isim,sep="")
        } else {
          print(paste("hri-error, jsim, iboot:", isim, iboot))
        }
      }
    }
    
    print(paste("done with simulation", isim, "of", nsim))
    hrepp$pmean = hrepp$mean/jsim
    hrepp$var = (hrepp$sq - jsim*(hrepp$pmean^2))/(jsim-1)
    hrepp$stdev = sqrt(hrepp$var)
    hrepp$cv = hrepp$stdev/hrepp$pmean
    sims = hrepp@data[,grep("sim",names(hrepp))]
    sims2 = sims > exc 
    hrepp$exc = rowSums(sims2)/jsim
  }
  if (full) fullgrid(hrepp) = TRUE
  attr(hrepp, "hri") = hriCalc(hrepp, habitat, pvals)
  hrepp
}
}
               
# hte = mhri(indispdforg, pa, indicators2 = indispdf, forecast = TRUE, nsim = 5)
               
