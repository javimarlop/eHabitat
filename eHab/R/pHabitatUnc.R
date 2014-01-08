pHabitatUnc = function(indicators = NULL, habitat = NULL, lmc){

if (inherits(habitat,"SpatialPolygons")) {
  habitat = rasterToPoints(rasterize(habitat, raster(indicators)), spatial = TRUE)  
}
inHA = over(habitat, geometry(indicators))
inHA = inHA[!is.na(inHA)]

polyData = indicators@data[inHA,]
cordm = coordinates(indicators[inHA,])
hm = mean(dist(cordm))


indicators$mDist = NA
indicators$mDistVar = NA

delta = solve(var(polyData, na.rm = TRUE))
mum = colMeans(polyData, na.rm = TRUE)
ndim = dim(indicators)[1]
dd = apply(coordinates(indicators), 1, FUN = function(x) mean(spDistsN1(cordm, x))) 
vars = names(lmc$data)
nvars = length(vars)

sigma = matrix(NA, ncol =  nvars, nrow = nvars)

variograms = matrix(NA, nrow = ndim, ncol = sum(1:nvars))
icol = 0
gammahm = 1:sum(1:nvars)
for (ivar in 1:nvars) {
  for (jvar in ivar:nvars) {
    id = ifelse(ivar == jvar, vars[ivar], paste(vars[ivar],vars[jvar], sep = "."))
    icol = icol + 1
    model = lmc$model[[id]]
    variograms[,icol] = variogramLine(model, dist_vector = dd)$gamma
    gammahm[icol] = variogramLine(model, dist_vector = hm)$gamma
  }
}

mDist = 1:ndim
mDistVar = 1:ndim
indata = as.matrix(indicators@data[,1:nvars])

mm = matrix(NA, nrow = ndim, ncol = nvars+13)
stest = exp(seq(-2,5,.3))
ssup = diff(c(0,stest,100))
sdists = (c(0,stest)+c(stest,100))/2
sims = 1-pchisq(sdists, nvars)
quants = c(0.25,0.75)
simi = mm
for (i in 1:ndim) {
  icol = 0
  mu = indata[i,] - mum
  t(mu) %*% delta %*% mu
  for (ivar in 1:nvars) {
    for (jvar in ivar:nvars) {
      icol = icol + 1
      if (ivar == jvar) {
        sigma[ivar, ivar] = variograms[i,icol]
      } else {
        sigma[ivar,jvar] = 2*variograms[i,icol]-gammahm[icol]
        sigma[jvar,ivar] = sigma[ivar,jvar]
      }
    }
  }
  delsig = delta %*% sigma
  mumd = t(mu) %*% delta %*% mu
  mumsig = t(mu) %*% delsig %*% delta %*% mu
  nfree = sum(diag(delsig))
  mm[i,1:nvars] = mu
  mm[i,nvars+1] = nfree
  mm[i,nvars+2] = mumd
  mm[i,nvars+3] = 2*sum(diag(delsig %*% delsig))
  mm[i,nvars+4] = 4 * mumsig
  mDist[i] = nfree + mumd
  mDistVar[i] = 2*sum(diag(delsig %*% delsig)) + 4 * mumsig
  if (i %% 100 == 0) print(i)
  sfreq = dchisq(sdists, nfree, mumd/2) 
  smean = sum(sfreq*sims*ssup)
  tdists = qchisq(quants, nfree, mumd/2)
  tsims = 1-pchisq(tdists, nvars, 0)
  mm[i,nvars+5] = 1-pchisq(mDist[i],nvars)
  mm[i,nvars+6] = smean
  mm[i,nvars+7] = mean(tsims)
  mm[i,nvars+8] = tsims[1]
  mm[i,nvars+9] = tsims[2]
  mm[i,nvars+(10:13)] = chiCont(delta, sigma, mu)
}
mm = as.data.frame(mm)
names(mm) = c(paste("mu",vars,sep = "_"), "nfree", "muest", "mvar1", "mvar2",
    "eSim", "smean", "tmean", "sim75", "sim25", "vavdiff", "mudiff", "vamudiff", "avadiff")

indicators$mDist = mDist
indicators$mDistVar = mDistVar
indicators$pHab = 1 - pchisq(indicators$mDist, nvars) 
indicators$pHabL = 1 - pchisq(indicators$mDist+sqrt(indicators$mDistVar), nvars) 
indicators$pHabH = 1 - pchisq(indicators$mDist-sqrt(indicators$mDistVar), nvars) 
indicators@data = cbind(indicators@data, mm)
indicators
}




makeLmcN = function (spdf, model = "Sph", range = NULL, nugget = NULL, maxvar = .2) 
{
    nd = dim(spdf)[1]
    if (nd > 5000) 
        spdf = spdf[sample(nd, 5000), ]
    gridded(spdf) = FALSE
    g.dummy = NULL
    for (ivar in 1:dim(spdf@data)[2]) {
        vname = names(spdf)[ivar]
        g.dummy = gstat(g.dummy, vname, as.formula(paste(vname, 
            "~ 1")), data = spdf)
    }
    dia = sqrt(eHab:::bbArea(bbox(spdf)))
    if (is.null(range)) 
        range = dia/5
    g.dummy = gstat(g.dummy, model = vgm(1, "Sph", range), fill.all = TRUE)
    lmc = fit.lmc(variogram(g.dummy), g.dummy, fit.ranges = TRUE)
    smax = max(unlist(sapply(lmc$model, FUN = function(x) x$psill)))
    for (im in 1:length(lmc$model)) lmc$model[[im]]$psill = lmc$model[[im]]$psill*maxvar^2
    lmc
}


chiCont = function(a, v, mu, debug.level = 0) {
  ava = a %*% v %*% a
  va = v %*% a
  vavav = va %*% va %*% v
  vav = va %*% v
  muavamu = t(mu) %*% ava %*% mu
  muamu = t(mu) %*% a %*% mu
  vavamu = va %*% va %*% mu
  vamu = va%*% mu
  if (debug.level >= 1) {
    print("vavav")
    print(vavav)
    print("vav")
    print(vav)
    print(" ")
    print("muavamu")
    print(muavamu)
    print("muamu")
    print(muamu)
    print(" ")
    print("vavamu")
    print(vavamu)
    print("vamu")
    print(vamu)
    print(" ")
    print("ava")
    print(ava)
    print("a")
    print(a)
  }
  vavdiff = sum(abs(vavav/vav - 1))
  mudiff = (muavamu-muamu)/(muavamu+muamu)
  vamudiff = sum(abs(vavamu/vamu-1))
  avadiff = sum(abs(ava/a-1))
  c(vavdiff, mudiff, vamudiff, avadiff)
}
#chiCont(delta, sigma, mu)
