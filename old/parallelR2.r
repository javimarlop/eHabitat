library(eHab)
library(rgeos)
library(raster)
library(sp)
library(rgdal)

removeTmpFiles(h=6)

parkdir = "parks/"
ecodir = "ecoregions/"
tiffdir = "tiffs/"
pngdir = "nopng"      # nopng gives no png output, otherwise they are dumped in this directory
ecoBuffer = 250000    # The size of the ecobuffer outside the ecoregion of the park - must be relative to the projection of the data
parkdir = "parks/"
ecodir = "ecoregions/"
wdpaid = "WDPA_ID"     # The column name giving the id of the parks/polygons
ecoID = "eco_id"      # The column name giving the ecoregion id in the ecoregions shapefile
inddir = "Variables"  # The directory of the tiff files with the indicator variables
minVar = 1

### Run it only when changing parks, ecoregions or indicators!
# if (!exists("indicators")) {
  # tiffs = list.files(inddir,pattern = ".tif", full.names = TRUE)
  # indicators = stack(tiffs)

# if (!exists("ecoregions")) {
  # shpfiles = list.files(ecodir, pattern = "shp")
  # shpfiles = strsplit(shpfiles, ".shp")
  # if (length(shpfiles) > 1) {
    # warning(paste("several shapefiles found in directory", ecodir, "using", shpfiles[[1]]))
  # }
  # shpfile = shpfiles[[1]]
  # ecoregions = readOGR(ecodir,shpfile) 
# }
# if (!exists("parks")) {
  # shpfiles = list.files(parkdir, pattern = "shp")
  # shpfiles = strsplit(shpfiles, ".shp")
  # if (length(shpfiles) > 1) {
    # warning(paste("several shapefiles found in directory", parkdir, "using", shpfiles[[1]]))
  # }
  # shpfile = shpfiles[[1]]
  # parks = readOGR(parkdir,shpfile) 
# }  

#gIntersects(ecoregions,parks,byid=T)->ep # matrix

### Run it only when changing parks, ecoregions or indicators!

# checking that all proj4strings are equal, transforming if not. 
# Also modifying the first one to make sure to have a valid string.
pstring = CRSargs(CRS(projection(indicators)))
if (!pstring == proj4string(indicators)) proj4string(indicators) = CRS(pstring)
if (!pstring == proj4string(parks)) parks = spTransform(parks, CRS(pstring))
if (!pstring == proj4string(ecoregions)) ecoregions = spTransform(ecoregions, CRS(pstring))

options(rasterMaxMemory = 3e22, chunksize=1e+7) # changed from 3e11 by @javier # testing maximum memory available


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

########

#1:14458->base # original
#base<-dim(ep)[2]
1:20->base # testing

require(rgeos)
require(rgdal)

# data.frame(matrix(NA,nrow=1,ncol=3))->results
# names(results)<-c('wdpa_id','eco_id','percentage')
# write.table(results,'results.csv',row.names=F,sep=';')

	nclus = 10
	library(parallel)
	cl = makeCluster(nclus, outfile = "")
#  clusterEvalQ(cl, library(eHab))
	out<-clusterApplyLB(cl, base, fun = ecoperpa2, ep = ep, parks = parks, ecoregions = ecoregions, results = results)
	stopCluster(cl)

####

  cl = makeCluster(nclus, outfile = "")
  clusterEvalQ(cl, library(eHab))
#  closeAllConnections()

# Using clusterApplyLB, as it seems to be the only parallel function that properly
# balances the workload on the processors, i.e., new tasks are distributed before
# all nodes have returned.
  s3 = results <- clusterApplyLB(cl, ecoregs$ecoregs, fun = mecohri2, 
      ecoregions = ecoregions, ecoregs = ecoregs$ecoregs, parks = parks, 
      indicators = indicators, pvals = seq(0.05, 1, 0.05), tiffdir = tiffdir, 
      pngdir = pngdir, hriRes = NULL, hriRes2 = NULL, minVar = minVar, 
      ecoBuffer = ecoBuffer, wdpaid = wdpaid, ecoID = ecoID)


 stopCluster(cl)

####
# ecoperpa2<-function(e, ep, parks, ecoregions, results){

# require(rgeos)
# require(rgdal)

# r<-0
# #data.frame(matrix(NA,nrow=1,ncol=3))->results
# #names(results)<<-c('wdpa_id','eco_id','percentage')

# length(which(ep[,e]))->lp

# if(lp>0){

# parks$WDPA_ID[which(ep[,e])]->pas

# for (p in 1:lp){

# r<-r+1

# # Example showing the parks related to the 3rd ecoregion
# #plot(ecoregions[3,])
# #plot(parks[ep[,3],],add=T,col='red')

# try(gIntersection(ecoregions[e,],parks[which(ep[,e])[p],])->ppee)

# #plot(p1e3,col='red')
# #plot(parks[which(ep[,3])[1],],add=T)
# #plot(ecoregions[3,],add=T)

# try(gArea(parks[which(ep[,e])[p],])->app)
# try(gArea(ppee)->appee)
# try(per_ep<-100*appee/app)

# try(results[r,]<-c(pas[p],ecoregions$eco_id[e],per_ep)) # no double assignment!
# print(paste('row:',r))
# # add a try for ecoregions without parks


# }
# print(paste('PA',p))
# }
# write.table(results,'results.csv',row.names=F,sep=';',col.names=F,append=T)
# #print('END')
# print(results)
# }

