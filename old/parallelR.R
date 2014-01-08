
library(rgeos)
library(rgdal)

parkdir = "parks/"
ecodir = "ecoregions/"

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

gIntersects(ecoregions,parks,byid=T)->ep # matrix

########

#1:14458->base # original
1:5->base # testing

require(rgeos)
require(rgdal)

data.frame(matrix(NA,nrow=1,ncol=3))->results
names(results)<-c('wdpa_id','eco_id','percentage')
write.table(results,'results.csv',row.names=F,sep=';')

	nclus = 10
	library(parallel)
	cl = makeCluster(nclus, outfile = "")
#  clusterEvalQ(cl, library(eHab))
	out<-clusterApplyLB(cl, base, fun = ecoperpa2, ep = ep, parks = parks, ecoregions = ecoregions, results = results)
	stopCluster(cl)

####

ecoperpa2<-function(e, ep, parks, ecoregions, results){

require(rgeos)
require(rgdal)

r<-0
#data.frame(matrix(NA,nrow=1,ncol=3))->results
#names(results)<<-c('wdpa_id','eco_id','percentage')

length(which(ep[,e]))->lp

if(lp>0){

parks$WDPA_ID[which(ep[,e])]->pas

for (p in 1:lp){

r<-r+1

# Example showing the parks related to the 3rd ecoregion
#plot(ecoregions[3,])
#plot(parks[ep[,3],],add=T,col='red')

try(gIntersection(ecoregions[e,],parks[which(ep[,e])[p],])->ppee)

#plot(p1e3,col='red')
#plot(parks[which(ep[,3])[1],],add=T)
#plot(ecoregions[3,],add=T)

try(gArea(parks[which(ep[,e])[p],])->app)
try(gArea(ppee)->appee)
try(per_ep<-100*appee/app)

try(results[r,]<-c(pas[p],ecoregions$eco_id[e],per_ep)) # no double assignment!
print(paste('row:',r))
# add a try for ecoregions without parks


}
print(paste('PA',p))
}
write.table(results,'results.csv',row.names=F,sep=';',col.names=F,append=T)
#print('END')
print(results)
}

