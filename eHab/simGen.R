#if (file.exists("hriRes.csv")) stop("Please, remove old csv files first") # Next step could be read from the csv files in order to restart the process when it finished (= reading processed ecoregions and filtering them out from the script).

# Example:
# Rscript --no-init-file simGen.R ecoID=ECO_ID inddir=/local0/skoiejo/hri/Variables 
#   minVar=function(x){x = x+1/2;c(0.1/x[1],0.1/x[2],5/x[3],5/x[4],0.1/x[5],200/x[6]+0.1,2/x[7],0.3/x[8]+0.1,10/x[9]+0.1)*x}

# setwd("/local0/skoiejo/hri")

# Required libraries
library(eHab)
library(rgeos)
library(rgdal)
#options(error = recover)
##rasterOptions(tmpdir='/local1/majavie/tmp/') # Not working in parallel! @javier

args = commandArgs(trailing = TRUE)
argss = strsplit(args, "=") # split arguments by '='

#wind = .Platform$OS.type == "windows"
sysinf = Sys.info()
comp = sysinf[["nodename"]]
login = sysinf[["login"]]

# This is to clean up the raster temp directory, which tends to be full rather fast on hanks,
# as the hard drive is small. 
removeTmpFiles(h=6) # files older than 6 hours are removed

# The default assumption is that the script started in the directory where the results
# should be stored, and that all input data and result data should be in subfolders.
wdir = "."
tiffdir = "tiffs/"
pngdir = "nopng"      # nopng gives no png output, otherwise they are dumped in this directory
ecoBuffer = 250000    # The size of the ecobuffer outside the ecoregion of the park - must be relative to the projection of the data
parkdir = "parks/"
ecodir = "ecoregions/"
wdpaid = "WDPA_ID"     # The column name giving the id of the parks/polygons
ecoID = "eco_id"      # The column name giving the ecoregion id in the ecoregions shapefile
inddir = "Variables"  # The directory of the tiff files with the indicator variables
minVar = 1
#nclus = 10 # 8 max


# The default parameters are the ones above. 
# These can be overridden by settings in the file settings.R,
# which will again be overridden by command line arguments.
if (file.exists("settings.R")) source("settings.R")
if (length(argss) > 0) {
  for (i in 1:length(argss)) {
    lis = argss[[i]]
    assign(lis[1], lis[2])
  }
}


# inputData.R is a file that can contain a script for reading the input data
# If not, it is assumed that the data is in the directories specified above
if (file.exists("inputData.R")) source("inputData.R")
if (!exists("ecoregions")) {
  shpfiles = list.files(ecodir, pattern = "shp")
  shpfiles = strsplit(shpfiles, ".shp")
  if (length(shpfiles) > 1) {
    warning(paste("several shapefiles found in directory", ecodir, "using", shpfiles[[1]]))
  }
  shpfile = shpfiles[[1]]
  ecoregions = readOGR(ecodir,shpfile) 
}
if (!exists("parks")) {
  shpfiles = list.files(parkdir, pattern = "shp")
  shpfiles = strsplit(shpfiles, ".shp")
  if (length(shpfiles) > 1) {
    warning(paste("several shapefiles found in directory", parkdir, "using", shpfiles[[1]]))
  }
  shpfile = shpfiles[[1]]
  parks = readOGR(parkdir,shpfile) 
}  

if (!exists("indicators")) {
  tiffs = list.files(inddir,pattern = ".tif", full.names = TRUE)
  indicators = stack(tiffs)
}

setwd(wdir)
# checking that all proj4strings are equal, transforming if not. 
# Also modifying the first one to make sure to have a valid string.
pstring = CRSargs(CRS(projection(indicators)))
if (!pstring == proj4string(indicators)) proj4string(indicators) = CRS(pstring)
if (!pstring == proj4string(parks)) parks = spTransform(parks, CRS(pstring))
if (!pstring == proj4string(ecoregions)) ecoregions = spTransform(ecoregions, CRS(pstring))

options(rasterMaxMemory = 3e100, chunksize=1e+20) # changed from 3e7 and 1e+7 by @javier


# Building a function for the minimum variance
if (minVar == 1) {                # minVarFuncDict can be added as an argument to the script or through 
                                  # inputData.R, but is usually defined on the lines below.
  if (!exists("minVarFuncDict")) minVarFuncDict = list(herb = "0.1/x[idx]", # where does 'minVarFuncDict' come from? @javier
                      tree = "0.1/x[idx]",
                      ndvi = "5/x[idx]",
                      ndwi = "5/x[idx]",
                      srtm_ramp2_slope_world = "0.1/x[idx]",
                      srtm_ramp2_world = "200/x[idx]+0.1",
                      epratio = "0.3/x[idx]+0.1",
                      bio = "2/x[idx]",
                      prec = "10/x[idx] + 0.1")
  flist = NULL
  for (li in 1:length(names(indicators))) {
    ln = names(indicators)[li]
    for (il in 1:length(minVarFuncDict)) {
      im = grep( names(minVarFuncDict)[il], ln) # pattern, input
      if (length(im) > 0) break
    }
    if (length(im) == 0) {
      warning(paste("could not match layername", ln, 
          "with any of the elements of minVarFuncDict, minVar function is unlikely to perform correctly"))
    } else {
      fadd = minVarFuncDict[il][[1]]
      fadd = gsub("idx", li, fadd) # pattern, replace by, input
      flist = c(flist, fadd)
    }
  }
  fstring = paste("{x = x+1/100000; c(", paste(flist, collapse = " , "), ")*x }")

  minVar = eval(bquote(function(x) .(parse(text =  fstring)[[1]])))
}


# Find ecoregions of the parks
  if (is.null(parks$ecoID)) {

    s01 = system.time(paeco <- over(SpatialPoints(coordinates(parks),
        proj4string = CRS(proj4string(ecoregions))), ecoregions)) # 'over' from library(sp)
    leids = which(is.na(paeco[ecoID]))

    print(s01[[1]])
    print(paste("will do second attempt to find ecoregions of ", length(leids),"parks"))
    pps = matrix(NA, nrow = length(leids), ncol = 10) 

    if(length(leids)>0){ #begin_a
    for (il in 1:length(leids)) {
      ptt = spsample(parks[leids[il],], 100, "regular")
      ptt = SpatialPointsDataFrame(ptt, data = data.frame(id = rep(il, length(ptt))))

      	if (il == 1)	{

        pts = ptt
      			} else pts = rbind(pts, ptt)

    				}

    s1 = system.time(pp <- over(pts, ecoregions)[,ecoID])[[1]]
    pts$ecoID =  pp

    for (il in 1:length(leids)) {
      pp = pts$ecoID[pts$id == il]      
      pp = pp[!is.na(pp)]
      ppd = c(table(pp))
      pp = unique(pp)
      print(paste(il, paste(pp, collapse = " ")))
      pps[il,1:(length(pp)+1)] = c(leids[il], pp)
      if (length(ppd) > 0) pps[il,6:(length(pp)+5)] = ppd
      if (sum(pps[il,6:10], na.rm = TRUE) > 10) paeco[leids[il], ecoID] = pps[il, which.max(pps[il,6:10])+1]
				}
    pps = cbind(pps, rowSums(pps[,6:10],na.rm = TRUE))

} # end_a
 # Assign ecoregion
    parks$ecoID = paeco[,ecoID]
  
  }


# FUNCTION DEFINITION ******
# Function to call the ecohri function
mecohri = function(ecoreg, ecoregions, ecoregs, ecoID,...) {
  removeTmpFiles(h=6)
  write(paste("Entering ecoreg", ecoreg, "with", Sys.getpid()), file = "mecohri.txt", append = TRUE)
  eids <- which(ecoregions@data[,ecoID] %in% ecoreg) # changed @javier
  if (length(eids) > 0) {
    leco <- ecoregions[eids,] # changed @javier
    print(dim(leco))
    ress <- try(ecohri(leco, ecoID = ecoID, ...)) # name changed to ress by @javier (26.07.13)
    #print(ress) # added by Javier (26.07.13)
    if (is(ress, "try-error")) ress = -999
  } else ress = NULL
  ids = which(ecoreg %in% ecoregs)
  
  	# write anyway results to csv file!
	write.table(ress$hriRes, file = paste(ecoreg,"_hriRes.csv",sep=''),sep=',', row.names=FALSE,col.names = FALSE)
	write.table(ress$hriRes2, file = paste(ecoreg,"_hriRes2.csv",sep=''),sep=',', row.names=FALSE,col.names = FALSE)
	write.table(ress$hriInRes, file = paste(ecoreg,"_hriInRes.csv",sep=''),sep=',', row.names=FALSE,col.names = FALSE)
	write.table(ress$hriInRes2, file = paste(ecoreg,"_hriInRes2.csv",sep=''),sep=',', row.names=FALSE,col.names = FALSE)
	write.table(ress$errors, file = "errors.csv",sep=',',append=TRUE, row.names=FALSE)
	
  if (!is.null(ress) && length(ress)>1 && dim(ress$hriRes)[1]!=0) { # changed by Javier (26.07.13) and updated (29.07.13) # ress!=-999
    print(paste("**** mecohri ", ecoreg, "tiles", attr(ress,"tr")$n, "parks", sum(!is.na(ress$hriRes[,1])), "****"))
    write(paste("Done ecoreg", ecoreg, "with", Sys.getpid(), 
                "tiles", attr(ress,"tr")$n, "parks", sum(!is.na(ress$hriRes[,1]))), file = "mecohri2.txt", append = TRUE)
	
  } else {
    print(paste("**** mecohri ", ecoreg, "NULL?" , "****"))
    write(paste("Done ecoreg", ecoreg, "with", Sys.getpid(), "NULL"), file = "mecohri2.txt", append = TRUE)
  }  
  ress
}
# FUNCTION DEFINITION ******


# FUNCTION DEFINITION ******
# Compute areas and sort the ecoregions according to their size 
getAreaPolygons = function(x) { 
  holes = unlist(lapply(x@Polygons, function(x) x@hole))
  areas = unlist(lapply(x@Polygons, function(x) x@area))
  area = ifelse(holes, -1, 1) * areas
  sum(area)
}
# FUNCTION DEFINITION ******


ecoregs = unique(ecoregions@data[,ecoID])
ecoregions$area =  sapply(ecoregions@polygons, getAreaPolygons)
ecoregs2 = data.frame(ecoregs = ecoregs, area = sapply(ecoregs, 
                                                       FUN = function(x) sum(ecoregions$area[ecoregions@data[,ecoID] == x])))
ecoregs2[-827,]->ecoregs3 
ecoregs_tmp = ecoregs3[order(ecoregs3$area, decreasing = TRUE),] # modified by @javier (26.07.13)
ecoregs_tmp2 = ecoregs_tmp[-ecoregs_tmp$ecoregs<0,] # added by @javier (26.07.13)

ecoregs = ecoregs_tmp2[c(2,4,10),] #changed by @javier
parks$area = sapply(parks@polygons, getAreaPolygons)/1e6


# Set up for parallelization if necessary
#if (nclus > 1 & length(ecoregions) > 3) { 
 nclus <- 3
  print('nclus 6; launching mecohri!') # @javier 29.07.13
  library(parallel)
  cl = makeCluster(nclus, outfile = "")
  clusterEvalQ(cl, library(eHab))
#  closeAllConnections()

# Using clusterApplyLB, as it seems to be the only parallel function that properly
# balances the workload on the processors, i.e., new tasks are distributed before
# all nodes have returned.
  s3 = system.time(results <- clusterApplyLB(cl, ecoregs$ecoregs, fun = mecohri, # call to mecohri @javier
      ecoregions = ecoregions, ecoregs = ecoregs$ecoregs, parks = parks, 
      indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
      pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
      ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID))


 stopCluster(cl)

ecoregs = ecoregs_tmp2[c(1,3,5:9,11:25),] #changed by @javier
parks$area = sapply(parks@polygons, getAreaPolygons)/1e6


# Set up for parallelization if necessary
#if (nclus > 1 & length(ecoregions) > 3) { 
 nclus <- 6
  print('nclus 6; launching mecohri!') # @javier 29.07.13
  library(parallel)
  cl = makeCluster(nclus, outfile = "")
  clusterEvalQ(cl, library(eHab))
#  closeAllConnections()

# Using clusterApplyLB, as it seems to be the only parallel function that properly
# balances the workload on the processors, i.e., new tasks are distributed before
# all nodes have returned.
  s3 = system.time(results <- clusterApplyLB(cl, ecoregs$ecoregs, fun = mecohri, # call to mecohri @javier
      ecoregions = ecoregions, ecoregs = ecoregs$ecoregs, parks = parks, 
      indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
      pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
      ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID))


 stopCluster(cl)
 #} #else {
  # print('nclus = 1; launching mecohri!') # @javier 29.07.13
  # We are not doing parallel processing
  # s4 = system.time(results <- sapply(ecoregs$ecoregs, FUN = mecohri, # call to mecohri @javier
                                           # ecoregions = ecoregions, ecoregs = ecoregs$ecoregs, parks = parks, 
                                           # indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
                                           # pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
                                           # ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID))
# }

ecoregs = ecoregs_tmp2[26:825,] #changed by @javier

nclus <- 8
  print('nclus 10; launching mecohri!') # @javier 29.07.13
  library(parallel)
  cl = makeCluster(nclus, outfile = "")
  clusterEvalQ(cl, library(eHab))
#  closeAllConnections()

# Using clusterApplyLB, as it seems to be the only parallel function that properly
# balances the workload on the processors, i.e., new tasks are distributed before
# all nodes have returned.
  s3 = system.time(results <- clusterApplyLB(cl, ecoregs$ecoregs, fun = mecohri, # call to mecohri @javier
      ecoregions = ecoregions, ecoregs = ecoregs$ecoregs, parks = parks, 
      indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
      pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
      ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID))


 stopCluster(cl)

# This is just to be able to check individual ecoregions with the ecohri function @jon; It skips the paralell and other funtions in the mecohri @javier
#if (FALSE) { 
#  s4 = system.time(results <- ecohri(ecoregions[ecoregions$ECO_ID %in% ecoregs$ecoregs,],  
#                                    parks = parks, 
#                                   indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
#                                   pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
#                                   ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID))
#}

# Reorganize results
#  first = TRUE
#  hriRes = NULL
#  hriInRes = NULL
#  hriRes2 = NULL
#  hriInRes2 = NULL
#  errors = list()
#  for (ieco in 1:length(results)) {
#    ress <<- results[[ieco]]
#    if (!is.null(ress) & !is(ress, "try-error")) {
#      hriRes <<- rbind(hriRes, ress$hriRes) # keep in wkspc @javier 02.08.13
#      hriRes2 <<- rbind(hriRes2, ress$hriRes2) # keep in wkspc @javier 02.08.13
#      hriInRes <<- rbind(hriInRes, ress$hriInRes) # keep in wkspc @javier 02.08.13
#      hriInRes2 <<- rbind(hriInRes2, ress$hriInRes2) # keep in wkspc @javier 02.08.13
#      errors <<- c(errors, ress$errors) # keep in wkspc @javier 02.08.13
#    }
#  }

print(s4)

system('cat header.txt *_hriRes.csv | grep -v - > hriRes.csv')
system('cat header.txt *_hriRes2.csv | grep -v - > hriRes2.csv')
system('cat header.txt *_hriInRes.csv | grep -v - > hriInRes.csv')
system('cat header.txt *_hriInRes2.csv | grep -v - > hriInRes2.csv')
print('DONE')
#write.csv(hriRes, file = "hriRes.csv")
#write.csv(hriRes2, file = "hriRes2.csv")
#write.csv(hriInRes, file = "hriInRes.csv")
#write.csv(hriInRes2, file = "hriInRes2.csv")
#write.csv(errors, file = "errors.csv")





