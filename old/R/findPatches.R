findPatches  = function(hrep, plot = TRUE) {
  hrep$reID = 0
  hrep$reID[hrep$replace] = c(1:sum(hrep$replace))
  hrepRaster = raster(hrep[,"reID"])
  fromCells = which(getValues(hrepRaster) > 0)
  hrepAdj = adjacency(hrepRaster,fromCells,toCells = c(1:ncell(hrepRaster)), directions=8)

  hrepAdj = cbind(hrepAdj,hrepRaster[hrepAdj[,1]])
  hrepAdj = cbind(hrepAdj,hrepRaster[hrepAdj[,2]])
  hrepAdj = hrepAdj[!is.na(hrepAdj[,4]),]

  hrepAdj = hrepAdj[order(hrepAdj[,1]),]
  cells = unique(hrepAdj[,1])
  icell = list()
  for (i in 1:length(cells)) icell[[i]] = which(hrepAdj[,1] == cells[i])
  jcell = list()
  for (i in 1:length(cells)) jcell[[i]] = which(hrepAdj[,2] == cells[i])
   
  ichange = 1
  while (ichange > 0) {
    ichange = 0
    for (i in 1:length(cells)) {
      if ((i %% 1000) == 0) print(paste(i,length(cells)))
      cmax = max(hrepAdj[,4][icell[[i]]])
      if (hrepAdj[icell[[i]][1],3] < cmax) {
        hrepAdj[icell[[i]],3] = cmax
        hrepAdj[jcell[[i]],4] = cmax
        ichange = ichange + 1
      }
    }
    print(ichange)
  }
  
  pID = unique(hrepAdj[,3])
  for (i in 1:length(pID)) hrepAdj[,3][hrepAdj[,3] == pID[i]] = i 
  hrepRaster[hrepAdj[,1]] = hrepAdj[,3]
  if (plot) plot(hrepRaster)
  print(unique(hrepRaster[cells]))
  print(length(unique(hrepRaster[cells])))
  print(length(cells))
  hrepRaster[cells]
}

