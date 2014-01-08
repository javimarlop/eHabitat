makeLmc = function(spdf, model = "Sph", range = NULL, nugget = NULL) {
#  rho = cor(spdf@data, use = "complete.obs")
#  rho[1:3,1:3] = c(1.0000000, 0.9774792, 0.9790330,
#                0.9774792, 1.0000000, 0.9552906,
#                 0.9790330, 0.9552906, 1.0000000)
  nd = dim(spdf)[1]
  if (nd > 5000) spdf = spdf[sample(nd,5000),]
  gridded(spdf) = FALSE
  narem = which(rowSums(is.na(spdf@data)) > 0)
  if (length(narem) > 0) spdf = spdf[-narem,]
  g.dummy = NULL
  for (ivar in 1:dim(spdf@data)[2]) {
    vname = names(spdf)[ivar]
    g.dummy = gstat(g.dummy, vname, as.formula(paste(vname, "~ 1")), data = spdf)
  }
  dia = sqrt(bbArea(bbox(spdf)))
  if (is.null(range)) range = dia/5
  g.dummy = gstat(g.dummy, model = vgm(1,"Sph", range), fill.all = TRUE)
  lmc = fit.lmc(variogram(g.dummy), g.dummy)
  smax = max(unlist(sapply(lmc$model, FUN = function(x) x$psill)))
  for (im in 1:length(lmc$model)) lmc$model[[im]]$psill = lmc$model[[im]]$psill/smax
  lmc
}



sumSim = function(isims, xbox, ybox, yy, boxes, yxw0, yxw1, isub, jsub, newdata,
                        nmax, nsim, md, lmcm, xad, debug.level) { 
  pred = NULL
  for (isim in isims) {
    if (debug.level > 1) print(isim)
    jbox = (isim-1) %% ybox + 1
    ibox = (isim-1) %/% ybox + 1
    jsim = (isub-1)*2 +jsub
    if (debug.level > 1) print(paste(isim, ibox, jbox))
    t0 = proc.time()[1]
    lbox = boxes[[ibox]][[jbox]]
    bb = bbox(lbox)
    xmin = bb[1,1]
    xmax = bb[1,2]
    ymin = bb[2,1]
    ymax = bb[2,2]
    xpmin = xmin - ifelse(isub ==1,yxw0,yxw1)
    xpmax = xmax + ifelse(isub ==1,yxw0,yxw1)
    ypmin = ymin - ifelse(jsub ==1,yxw0,yxw1)
    ypmax = ymax + ifelse(jsub ==1,yxw0,yxw1)
    if (is.null(md))  md = yxw0/5
    lnew = newdata[newdata$id2 == isim*10+jsim,]
    
    if (debug.level > 1) {print(paste("isim", length(lnew))); print(bbox(lnew))}
    if (length(lnew) > 0) {
      lids = lnew$id
      SpP2 = makePolygon(xcor = c(xpmin,xpmax), ycor = c(ypmin,ypmax), pstring = proj4string(newdata))
      
      #ids = unlist(over(gBuffer(SpP, width = yx), geometry(yy), returnList = TRUE))
      ids = unlist(over(SpP2, geometry(yy), returnList = TRUE))
      obs = yy[ids,-dim(yy)[2]]
      dup = which(duplicated(rbind(coordinates(lnew), coordinates(obs))))
      if (length(dup) > 0) obs = obs[-(dup-length(lnew)),]
      if (length(lnew) == 0) break
      nvar = dim(obs)[2]
      vars = names(obs)
      if (debug.level > 1) print(paste("vars", paste(vars, collapse = " ")))
      if (debug.level > 1) print(paste("lmcm names",  paste(names(lmcm)[-grep("\\.", names(lmcm))], collapse = " ")))  
      g.dummy = NULL
      for (i in 1:nvar) {
        g.dummy <- gstat(g.dummy, vars[i], formula = as.formula(paste(vars[i], "~1")), 
            beta = mean(obs@data[,vars[i]]), data = obs, nmax = nmax, maxdist = md)
      }
      g.dummy = gstat(g.dummy, model = lmcm, fill.all = TRUE)
      g.dummy$model = lmcm
      if (debug.level > 1) print(g.dummy)
      s1 = system.time(zz <- try(predict.gstat(g.dummy, newdata = geometry(lnew), 
             nsim = nsim, debug.level = debug.level, nmax = nmax)))[1]
      if (debug.level > 1) print(paste("names(zz)",names(zz)))
      if (debug.level > 1) print(paste("vars",names(vars)) )
      t1 = proc.time()[1]
      if (debug.level >= 1) print(paste("time", isub,jsub, ibox, jbox,  dim(obs)[1],length(lnew),
         round(s1,2), round(t1-t0,2)))
      if (!is(zz, "try-error")) {
        zz@data = cbind(zz@data, id= lids)
      } else zz = NULL
      pred = rbind(pred, as.data.frame(zz))
    }
    if (debug.level > 1) print(paste("sumSim:", isim,length(lnew), dim(pred)[1], dim(pred)[1]- length(lnew)))
  }
  pred
}


mvSim = function(lmc, newdata, nsim = 1, debug.level = -1, nmax = 10, xbox = 4, 
    ybox = 4, md = NULL, nclus = 1) {
#
  if (is(newdata, "SpatialPoints")) newdata = SpatialPointsDataFrame(newdata, data = data.frame(id = rep(NA, length(newdata))))
  newdata$id = 1:length(newdata)
  lmcm = lmc$model
  vars = names(lmcm)[-grep("\\.", names(lmcm))]
  nvar = length(vars)
  ispdf = spsample(newdata, 10000, "regular")
  dia = sqrt(bbArea(bbox(ispdf)))
  xad = dia/sqrt(length(newdata))
  gg.dummy = NULL
  for (i in 1:nvar) {
    gg.dummy <- gstat(gg.dummy, vars[i], formula = as.formula(paste(vars[i], "~1")), 
     locations = ~x+y, dummy = TRUE,  beta = 0, nmax = nmax, maxdist = Inf)
  }
  gg.dummy$model = lmcm
  yy <- try(predict.gstat(gg.dummy, newdata = ispdf, nsim = nsim, debug.level = -1))
  names(yy) = vars
  yy$id = NA
  ydim = dim(yy)[1]
  
#
  xd = diff(bbox(newdata)[1,])/xbox
  xlms = bbox(newdata)[1,1]+c(0:xbox)*xd
  yd = diff(bbox(newdata)[2,])/ybox
  ylms = bbox(newdata)[2,1]+c(0:ybox)*yd
  dia2 = sqrt((xd/xbox)^2+(yd/ybox)^2)
  yxw0 =  dia2/2 
  nnd = sqrt(bbArea(bbox(newdata))/length(newdata))
  yxw1 = nnd*2
  boxes = list()
  il = 0
  res = as.data.frame(yy)
  newdata$id2 = 0
  for (ibox in 1:xbox) {
    lbox = list()
    for (jbox in 1:ybox) {
      il = il + 1
      xmin = xlms[ibox]
      xmax = xlms[ibox+1]
      ymin = ylms[jbox]
      ymax = ylms[jbox+1]
      SpP = makePolygon(xcor = c(xmin, xmax), ycor = c(ymin, ymax), pstring = proj4string(newdata))
      lbox[[jbox]] = SpP
      lids = unlist(over(SpP, geometry(newdata), returnList = TRUE))
      idbas = il*10
      xs = c(xmin,mean(c(xmin,xmax)), xmax)
      ys = c(ymin,mean(c(ymin,ymax)), ymax)
      xadd = ifelse(ibox == xbox, xad/1e6,0)
      xsub = ifelse(ibox == 1, xad/1e6,0)
      yadd = ifelse(jbox == ybox, xad/1e6,0)
      ysub = ifelse(jbox == 1, xad/1e6,0)
      jl = 0
      for (isub in 1:2) {
        for (jsub in 1:2) {
          jl = jl + 1
          xmin = xs[isub]-xsub
          xmax = xs[isub+1]+xadd
          ymin = ys[jsub]-ysub
          ymax = ys[jsub+1]+yadd
          xpmin = xmin - ifelse(isub ==1,yxw0,yxw1)
          xpmax = xmax + ifelse(isub ==1,yxw0,yxw1)
          ypmin = ymin - ifelse(jsub ==1,yxw0,yxw1)
          ypmax = ymax + ifelse(jsub ==1,yxw0,yxw1)
          if (is.null(md))  md = yxw0/5
          SpP = makePolygon(xcor = c(xmin, xmax), ycor = c(ymin, ymax), pstring = proj4string(newdata))
          lids = unlist(over(SpP, geometry(newdata), returnList = TRUE))
          newdata$id2[lids] = idbas + jl
        }
      }
    }
    boxes[[ibox]] = lbox
  }
  for (isub in 1:2) {
    for (jsub in 1:2) {
### Possibility for paralellization 
      nPred = xbox*ybox

      if (nclus > 1) {
        if (!suppressMessages(suppressWarnings(require(doParallel)))) 
           stop("nclus is > 1, but package doParalell is not available")
        if (isub == 1 & jsub == 1) {
          clus <- c(rep("localhost", nclus))
          cl <- makeCluster(clus, type = "SOCK")
          registerDoParallel(cl, nclus)
          clusterEvalQ(cl, library(eHab))
          splt = rep(1:nclus, length.out = nPred)
          clist = split(1:nPred, splt)
          iclus = 1  # to avoid check error
        }
        pred <- foreach(iclus = 1:nclus, .combine = rbind) %dopar% {
         sink("log.txt", append=TRUE)
          sumSim(clist[[iclus]], xbox, ybox, yy, boxes, yxw0, yxw1, isub, jsub, newdata,
                        nmax, nsim, md, lmcm, xad, debug.level) 
        }
      } else {
        pred = NULL
        for (isim in 1:nPred) {
          ploc = sumSim(isim, xbox, ybox, yy, boxes,yxw0, yxw1, isub, jsub, newdata,
                        nmax, nsim, md, lmcm, xad, debug.level) 
          pred = rbind(pred, ploc)
        }
      }
      names(pred) = names(res)
      pred = SpatialPointsDataFrame(pred[,1:2],pred[,3:(3+nvar)], proj4string = CRS(proj4string(yy)))
      yy = rbind(yy, pred)
    }
  }

  if (nclus > 1) stopCluster(cl)

  yy = yy[-(1:ydim),]
  if (gridded(newdata)) {
    gridded(yy) = TRUE
  } else yy = yy[order(yy$id),]
  yy[,-dim(yy)[2]]
}
