ecohri = function(ecoregions = NULL, parks, indicators, ecoID = names(ecoregions)[1],
  pvals = seq(0.05, 1, 0.05), tiffdir = "tiffs", pngdir = "nopng",
  hriRes = NULL, hriRes2 = NULL, hriInRes = NULL, hriInRes2 = NULL, 
  minVar = NULL, istart = 1, tparks = 0, ecoBuffer = 1,
  wdpaid = "wdpaid", keepErrors = FALSE, pamax = 60) {

  s1 = s2 = s3 = s4 = s5 = s6 = 0
  exlim = function(x,y) {
    if (x@xmin < y@xmin) x@xmin = y@xmin
    if (x@xmax > y@xmax) x@xmax = y@xmax
    if (x@ymin < y@ymin) x@ymin = y@ymin
    if (x@ymax > y@ymax) x@ymax = y@ymax
    x
  }
  nparks = dim(parks)[1]
  if (is.null(ecoregions)) {
    parks$ecoID = 1
  }
  # if (is.null(parks$ecoID)) {
    # s01 = system.time(paeco <- over(SpatialPoints(coordinates(parks),
        # proj4string = CRS(proj4string(ecoregions))), ecoregions))
    # leids = which(is.na(paeco[ecoID]))
    # print(s01[[1]])
    # if (length(leids > 0)) {
      # print(paste("will do second attempt to find ecoregions of ", length(leids),"parks"))
      # pps = matrix(NA, nrow = length(leids), ncol = 10) 
      # for (il in 1:length(leids)) {
        # ptt = spsample(parks[leids[il],], 100, "regular")
        # ptt = SpatialPointsDataFrame(ptt, data = data.frame(id = rep(il, length(ptt))))
        # if (il == 1) {
          # pts = ptt
        # } else pts = rbind(pts, ptt)
      # }
      # s1 = system.time(pp <- over(pts, ecoregions)[,ecoID])[[1]]
      # pts$ecoID =  pp
      # for (il in 1:length(leids)) {
        # pp = pts$ecoID[pts$id == il]      
        # pp = pp[!is.na(pp)]
        # ppd = c(table(pp))
        # pp = unique(pp)
        # print(paste(il, paste(pp, collapse = " ")))
        # pps[il,1:(length(pp)+1)] = c(leids[il], pp)
        # if (length(ppd) > 0) pps[il,6:(length(pp)+5)] = ppd
        # if (sum(pps[il,6:10], na.rm = TRUE) > 10) paeco[leids[il], ecoID] = pps[il, which.max(pps[il,6:10])+1]
      # }
      # pps = cbind(pps, rowSums(pps[,6:10],na.rm = TRUE))
    # }
 # # Assign ecoregion
# #
    # parks$ecoID = paeco[,ecoID]
  # }
bb = bbox(indicators)
#
if (is.null(hriRes)) {
  hriRes = data.frame(matrix(NA, ncol = length(pvals)+2, nrow = dim(parks)[1]))
  names(hriRes) = c("paID", "ecoID", paste("hri",pvals*100,sep = "_"))
  hriInRes = hriRes
  if (is.null(hriRes2) & !is.null(minVar)) {
    hriRes2 = hriRes
    hriInRes2 = hriRes
  }
}
#
#################
ecoregsp = unique(parks$ecoID)
ecoregse = unique(ecoregions@data[,ecoID])
ecoregs = ecoregse[ecoregse %in% ecoregsp]
if (length(ecoregs) == 0) return(NULL)
errors = vector("list", length(parks))
print(ecoregs)

gc()

for (iie in istart:length(ecoregs)) {
  ieco <- ecoregs[iie] # IMPORTANT TO CREATE IDENTIFIABLE TMP FILES! @javier
  print(paste("looping over ecoregs", iie, "of", length(ecoregs)))
  if (is.na(ieco) || ieco < -1) next
  ecor = ecoregions[ecoregions@data[,ecoID] == ieco,] # creation of 'ecor' @javier
  if (length(ecor@polygons[[1]]@Polygons) > 100 | length(ecor@polygons) > 100) ecor = gConvexHull(ecor)
  ecor = createSPComment(ecor) # library 'rgeos' @javier (just to clean up the polygons @jon)
  ecorr = gBuffer(ecor, width = ecoBuffer) # for buffers @javier
  if (is.null(intersect(extent(indicators), extent(ecorr)))) next # CONTINUE HERE! @javier
                                                        # Just to be clear: Continue with next iteration of iie
  execor = exlim(extent(ecorr),extent(indicators)) 
  				  	# fn7<-rasterTmpFile()
					# print(fn7)
  fnn = indicators[[1]]@file@name
  print(fnn)
  print(GDALinfo(fnn))
  s1 = round(system.time(indicators2 <- crop(indicators, execor))[[1]],2) # creation of target raster? @javier
                  # This is cropping the global rasterStack of indicators to a rectangular one that is just 
                  # big enough to include the ecoregion + buffer @jskoien
				  	 #fn5<-rasterTmpFile()
					#print(fn5)
  s2 = round(system.time(indicators2 <- mask(indicators2, ecorr))[[1]],2) # MASK, 
                  # i.e., blanking (setting to NA) all pixels in the cropped rasterStack which 
                  # are outside the ecoregion + buffer @jskoien
					 #fn6<-rasterTmpFile()
					 #print(fn6)
  gc()

  cpm = TRUE # canProcessInMemory(indicators2, n = 15) @javier
  s22 = proc.time()
  paidss = which(parks$ecoID == ieco)
  pai = ceiling(length(paidss)/pamax)

####################################################################################################

  for (ipai in 1:pai) { # starts looping over subsets of parks up to line 330
    paids = paidss[((ipai-1)*pamax+1): (ipai*pamax)]
    paids = paids[!is.na(paids)]
    print(paste("looping over subsets of parks", ipai, "of", pai))
    tmpfiles = list()
    tmpfiles2 = list()

    if (!cpm) { # if cpm is FALSE, so the processing cannot take place in the memory @javier
	print('cpm FALSE javier')
      fls <- NULL # @javier
      for (i in 1:length(paids)) {
	    
        tmpfiles[[i]] = raster(indicators2)
        tmpfiles[[i]] = writeStart(tmpfiles[[i]], filena=rasterTmpFile()) # @javier
		print('raster tmp created')
        f1 = tmpfiles[[i]]@file@name
        fls <- c(fls, f1)
        f2 =  gsub(".grd", ".gri", f1)
        fls <- c(fls, f2)
		  	#fn9<-rasterTmpFile()
			#print(fn9)
        if (!is.null(minVar)) {
          tmpfiles2[[i]] = raster(indicators2)
          tmpfiles2[[i]] = writeStart(tmpfiles2[[i]], filena=rasterTmpFile()) # @javier
		  print('raster tmp mst created')
          f1 = tmpfiles2[[i]]@file@name
          fls <- c(fls, f1)
          f2 =  gsub(".grd", ".gri", f1)
          fls <- c(fls, f2)
		  print(fls)
		   	#fn10<-rasterTmpFile()
			#print(fn10)
        }
      }
      tr <- blockSize(indicators2, chunksize = 5e8) # do the tiles! @javier
      # @jskoien This sets the block size (number or chunks, where it starts etc)
    } else {
	print('cpm TRUE javier')
      tr = list(n = 1, row = 1, nrows = dim(indicators2)[1]) # don't do the tiles! @javier
      # @jskoien This creates a tr-object when cpm = TRUE, i.e., when the processing can take place in the memory. 
	    			#fn11<-rasterTmpFile()
					#print(fn11)
    }

    cols = dim(indicators2)[2]
    hriRes[paids,] = NA
    hriInRes[paids,] = NA
    if (!is.null(minVar)) {hriRes2[paids,] = NA; hriInRes2[paids,] = NA}
    habMeanCov = list()
    inHA = list()
    paErr = NULL
    for (it in 1:tr$n) { # what is n if not 1? (after lines 116 and 118)
                         # @jskoien n is different from 1 when the ecoregion cannot be processed in memory and 
                         # is accessed in chunks/tiles.
      # n is the 
    	print(paste("looping over blocks", it, "of", tr$n))
      t1 = proc.time()[[1]]
      icells = cellFromRowColCombine(indicators2,tr$row[it]: (tr$row[it]+tr$nrows[it]-1), 1:cols)
      ind <- getValuesBlock(indicators2, row=tr$row[it], nrows=tr$nrows[it])
      gc()
      indispdf = SpatialPointsDataFrame(xyFromCell(indicators2,icells), data = data.frame(ind))
      proj4string(indispdf) = CRS(proj4string(parks))
      t2 = proc.time()[[1]]
      rm(icells)
      rm(ind)
#    if (it == tr$n & it > 1) rm(indicators2)

      gc()

      gridded(indispdf) = TRUE
      for (jpark in 1:length(paids)) {
        print(paste("looping over parks", jpark, "of", length(paids)))
        ipark = paids[jpark]
        if (ipark %in% paErr) break
        t3 = proc.time()[[1]]
        pa = parks[ipark,]
        paid = pa@data[,wdpaid]
        print(paste("ecoregion, park", ieco, paid))

# Extract the indicators for the park

        if (it == 1) {       
          pdat = extract(indicators2, pa, df = TRUE, small = TRUE)[,-1]
		    		#fn8<-rasterTmpFile()
					#print(fn8)
          hmc = try(habMeanCovPoly(polyData = pdat))
          if (is(hmc, "try-error")) hmc$polyData = matrix(NA) 
          habMeanCov[[jpark]] = hmc
          tparks = tparks + 1
        }
        gc()
        pDat = habMeanCov[[jpark]]$polyData
        if (dim(pDat)[2] > 1 && sum(rowSums(is.na(pDat)) == 0) > dim(pDat)[2])  {
          polyMean = colMeans(habMeanCov[[jpark]]$polyData,na.rm = TRUE)
          polyCov = var(habMeanCov[[jpark]]$polyData, na.rm = TRUE)
          if (!is.null(minVar)) {
            sdmin = minVar(polyMean)
            pCov = polyCov
            if (any(diag(pCov) == 0)) pCov = pCov +  min(abs(polyCov[abs(polyCov) > 0]))/1000
            for (ivar in 1:length(polyMean)) {
              if (sqrt(diag(pCov)[ivar]) < sdmin[ivar]) {
                ovar = diag(pCov)[ivar]
                pCov[ivar,] = pCov[ivar,]*(sdmin[ivar]/sqrt(ovar))
                pCov[,ivar] = pCov[,ivar]*(sdmin[ivar]/sqrt(ovar))
              }
            }
          }
          s3 = round(system.time(hriLoc <- hri(indispdf, 
               meanCovPoly = habMeanCov[[jpark]], pvals = pvals, na.rm = TRUE))[[1]],2)
#          
          if (is.factor(paid)) paid = levels(paid)[paid]
          if (is(hriLoc, "try-error")) {
            paErr = c(paErr, paid)
            hriRes[ipark,] = -999; hriInRes[ipark,] = -999
            errors[[ipark]] = c(errors[[ipark]] ,hriLoc)
            if (!is.null(minVar)) {hriRes2[ipark,] = -999; hriInRes2[ipark,] = -999}
          } else {
            pHab1 = hriLoc[,"pHab"]
            pHab1$pHab[pHab1$pHab < 1e-9] = NA
            if (!is.null(minVar)) {
              hmc = habMeanCov[[jpark]]
              hmc$covPoly = pCov
              hmc$nonNum = NULL
              s4 = round(system.time(hriLoc2 <- hri(indispdf,  # modified by Javier (26.07.13)
                       meanCovPoly = hmc, pvals = pvals, na.rm = TRUE))[[1]],2)
              pHab2 = hriLoc2[,"pHab"]
              pHab2$pHab[pHab2$pHab < 1e-9] = NA
            }
            if (cpm) { # process taking place in memory! @javier
              if (length(pa@polygons[[1]]@Polygons) -
                sum(sapply(pa@polygons[[1]]@Polygons,FUN = function(x) x@hole)) > 100 & FALSE) {
                pa = createSPComment(pa)
                pa = gSimplify(pa, gridparameters(indispdf)$cellsize[1]/2)
                pa = SpatialPolygonsDataFrame(pa, parks[ipark,]@data)
              }
              if (pngdir != "nopng") {
                png(paste(pngdir,paid, ".png", sep = ""))
                print(spplot(pHab1, "pHab", col.regions = bpy.colors(),   panel = function(x,y, ...){
                  panel.gridplot(x,y, ...)
                  sp.polygons(pa, col="red")}))
                dev.off()
              }
              s5 = round(system.time(writeGDAL(pHab1, 
                    fname = paste(tiffdir,paid, "_", ieco, ".tif", sep = "")))[[1]],2)
					print('done tiff')
					 #fn3<-rasterTmpFile()
					 #print(fn3)
              hric = attr(hriLoc, "hri")
              if (is.null(hric) | is.na(hric)) hric = hriCalc(hriLoc, pa, pvals)
              hriRes[ipark,] = c(list(paid, ieco), as.list(hric$out_pHab))
              hriInRes[ipark,] = c(list(paid, ieco), as.list(hric$in_pHab))               
              if (!is.null(minVar)) {
                if (pngdir != "nopng") {
                  png(paste(pngdir, "mst", paid,  "_", ieco, ".png", sep = ""))
                  print(spplot(pHab2, "pHab", col.regions = bpy.colors(),
                    panel = function(x,y, ...){
                    panel.gridplot(x,y, ...)
                    sp.polygons(pa, col="red")}))
                  dev.off()
                }
                writeGDAL(pHab2, fname = paste(tiffdir, "mst", paid,  "_", ieco, ".tif", sep = ""))
				 #fn4<-rasterTmpFile()
				 #print(fn4)
				 print('done mst tiff')
                hric = attr(hriLoc2, "hri")
                hric = NULL
                if (is.null(hric)) hric = hriCalc(hriLoc2, pa, pvals)
                hriRes2[ipark,] = c(list(paid, ieco), as.list(hric$out_pHab))
                hriInRes2[ipark,] = c(list(paid, ieco), as.list(hric$in_pHab))
              }
            } else { # process NOT taking place in memory! @javier
              rloc1 = raster(pHab1)
              tmpfiles[[jpark]] = writeValues(tmpfiles[[jpark]], getValues(rloc1),tr$row[it]) # creation of tmp files do not work! @javier
			  print('writting tmp rasters')
              if (!is.null(minVar)) {
                rloc2 = raster(pHab2)
                tmpfiles2[[jpark]] = writeValues(tmpfiles2[[jpark]], getValues(rloc2),tr$row[it]) # creation of tmp files do not work! @javier
				print('writting tmp mst rasters')
              }
              s5 = 0
            }
            t4 = proc.time()[[1]]
            tts = c(cpm, s1,s2,s3,s4,s5,round(t2-t1, 2),round(t4-t3,2))
            print(paste("eco", iie, "of", length(ecoregs), "it", it, "of", tr$n,
                "pa", jpark, "of", length(paids), "tpark", tparks, "of", nparks,
            paste(sprintf("%3.3f",ifelse(cpm,hric[floor(mean(1:length(pvals))),2],0)), collapse = " ")))      
            print(tts)
          }
        } else {
          hriRes[ipark,] = -999
          hriInRes[ipark,] = -999
          if (!is.null(minVar)) {hriRes2[ipark,] = -999 ; hriInRes2[ipark,] = -999}
          errors[[ipark]] = c(errors[[ipark]] ,paste("dim(pDat)[2]", dim(pDat)[2], 
              "dim(pDat)[2] > 1", dim(pDat)[2] > 1, 
              "sum(rowSums(is.na(pDat)) == 0)", sum(rowSums(is.na(pDat)) == 0)))
          paErr = c(paErr, paid)
        }
      }
    }

    if (!cpm) { #bracket1_init	# if cpm is FALSE, so the processing cannot take place in the memory @javier

      for (jpark in 1:length(paids)) { #bracket2_init
        ipark = paids[jpark]
        paid = parks@data[ipark,wdpaid]
        tmpfiles[[jpark]] = writeStop(tmpfiles[[jpark]] )

        if (!paid %in% paErr) { #bracket3_init
          projection(tmpfiles[[jpark]]) = projection(indicators2)
          names(tmpfiles[[jpark]]) = "pHab"
          hric = hriCalc(tmpfiles[[jpark]], parks[paids[jpark],], pvals)
          if (is.factor(paid)) paid = levels(paid)[paid]
          hriRes[ipark,] = c(list(paid, ieco), as.list(hric$out_pHab))
          hriInRes[ipark,] = c(list(paid, ieco), as.list(hric$in_pHab))
        } #bracket3_end

        if (!is.null(minVar)) { #bracket4_init
          names(tmpfiles2[[jpark]]) = "pHab"
          tmpfiles2[[jpark]] = writeStop(tmpfiles2[[jpark]] )
          if (!paid %in% paErr) { #bracket5_init
            projection(tmpfiles2[[jpark]]) = projection(indicators2)
            hric = hriCalc(tmpfiles2[[jpark]], parks[paids[jpark],], pvals)
            hriRes2[ipark,] = c(list(paid, ieco), as.list(hric$out_pHab))
            hriInRes2[ipark,] = c(list(paid, ieco), as.list(hric$in_pHab))
          } #bracket5_end
        } #bracket4_end

        if (!paid %in% paErr) { #bracket6_init

#          if (FALSE) {                              # is it OK to comment this? @javier
#            if (pngdir != "nopng") {
#              png(paste(pngdir, paid,  "_", ieco, ".png", sep = ""))
#              print(spplot(tmpfiles[[jpark]], col.regions = bpy.colors(),   panel = function(x,y, ...){
#                panel.gridplot(x,y, ...)
#                sp.polygons(pa, col="red")}))
#              dev.off()
#            }
#            if (!is.null(minVar)) {
#              if (pngdir != "nopng") {
#                png(paste(pngdir, "mst", paid, "_", ieco, ".png", sep = ""))
#                print(spplot(tmpfiles2[[jpark]], col.regions = bpy.colors(),   panel = function(x,y, ...){
#                  panel.gridplot(x,y, ...)
#                  sp.polygons(pa, col="red")}))
#                dev.off()
#              }
#            }
#          } # can I comment up to here? @javier
          # Sure, go ahead :-)
          # There were some problems making pngs from tiled rasters, and as I never needed them I 
          # outcommented them. @jskoien

          writeRaster(tmpfiles[[jpark]], filename = paste(tiffdir, paid, "_", ieco, ".tif", sep = "", SPARSE_OK=TRUE),
                     format = "GTiff", overwrite = TRUE)
					 #fn1<-rasterTmpFile()
					 print('done tmp tiff')
          if (!is.null(minVar)) writeRaster(tmpfiles2[[jpark]],
                     filename = paste(tiffdir, "mst", paid, "_", ieco, ".tif", sep = ""),
                     format = "GTiff", overwrite = TRUE)
					 print('done tmp mst tiff')
					 #fn2<-rasterTmpFile()
					 #print(fn2)
# deleting was here

        } #bracket6_end

	#fls<<-list.files('/tmp/R_raster_tmp/majavie/',pattern=paste(ieco)) # @javier /local1/majavie/tmp/
	#flss<<-paste('/tmp/R_raster_tmp/majavie/',fls,sep='') # @javier /local1/majavie/tmp/
	#file.remove(flss)-> kk # added by @javier (30.07.13)
	#print(paste('file ', fls,' removed:',kk,sep='')) # added by @javier (30.07.13)

      } #bracket2_end
  file.remove(fls) -> kk
  print(paste('file ', fls,' removed:',kk,sep=''))
    } #bracket1_end
  
  } # This bracket is the end of the for-loop starting at line 96 (looping over 60 parks at a time) @jskoien

####################################################################################################

  s33 = proc.time()
  s44 = s33-s22
  print(s44)
  if (exists("trr")) trr = c(trr, tr) else trr = tr
} #@javier
if (!keepErrors) {
  if (sum(is.na(hriRes[,1])) > 0) {
    hriRes = hriRes[!is.na(hriRes[,1]),]
    hriInRes = hriInRes[!is.na(hriInRes[,1]),]
    if (!is.null(minVar)) {
      hriRes2 = hriRes2[!is.na(hriInRes2[,1]),]
      hriInRes2 = hriInRes2[!is.na(hriInRes2[,1]),]
    }
    errorsall = errors
    errors = list()
    j = 0
    for (i in 1:length(errorsall)) if (!is.null(errorsall[[i]])) {j = j+1; errors[[j]] = c(i,errorsall[[i]])} 
  }
  if (sum(hriRes[,1] < 0) > 0) {
    hriRes = hriRes[hriRes[,1] >0,]    
    hriInRes = hriInRes[hriRes[,1] >0,]    
    if (!is.null(minVar)) {
      hriRes2 = hriRes2[hriRes[,1] >0,]    
      hriInRes2 = hriInRes2[hriRes[,1] >0,]    
    }
  }
}
res1 <- list(hriRes = hriRes,  hriInRes = hriInRes, errors = errors) # name changed to res1 by @javier
if (!is.null(minVar)) {
  res1$hriInRes2 = hriInRes2
  res1$hriRes2 = hriRes2
}
if (exists("trr")) attr(res1, "tr") = trr
res1
#file.remove(c(fn7,fn8,fn9,fn10,fn11))->kk2
#print(paste('temp files removed: ',kk2,sep=''))
#removeTmpFiles(h=6)
}


