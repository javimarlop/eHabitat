library(eHab)
if (FALSE) {
options(error = recover)
library(rgdal)
mamms = readOGR("E:/data/Species/Mammterr","MAMMTERR")
wdpa = readOGR("E:/data/parks/wdpa","wdpat")
twpa = wdpa[wdpa$name %in% c("Serengeti", "Selous", "Kruger", "Kafue", "Salonga", "Dja", "Minkebe", "Kigosi")  ,]
rm(wdpa)
gc()
continents = readOGR("e:/data/shapes/continents", "continent")
africa = continents[continents$CONTINENT == "Africa",]


ares = list()
for (ip in c(1:6,8:9)) {
#pa = twdpa[twdpa$name == "Serengeti",]
pa = twpa[ip,]

papts = spsample(pa,100,"regular")
mpa = over(mamms, papts)

maml = mamms[which(!is.na(mpa)),]


epratio = raster("F:/WorldClim/Current/10m/LRM_Bioclim/wc_10m_epratio.tif")
bio = raster("F:/WorldClim/Current/10m/LRM_Bioclim/wc_10m_bio.tif")
prec = raster("F:/WorldClim/Current/10m/LRM_Bioclim/wc_10m_prec.tif")
indicators = stack(epratio, bio, prec)
indicators = crop(indicators, extent(africa))    
indispdf = rasterToPoints(indicators, spatial = TRUE)
gridded(indispdf) = TRUE
proj4string(indispdf) = CRS(proj4string(maml))

dirr = "F:/WorldClim/Forecast/10m/LRM_Bioclim/"
indis = list("Current" = indispdf)
for (year in c("2020", "2050", "2080")) {
  epratio = raster(paste(dirr,"wc_10m_HADCM3_A2a_", year, "_epratio.tif", sep = ""))
  bio = raster(paste(dirr, "wc_10m_HADCM3_A2a_", year, "_bio.tif", sep = ""))
  prec = raster(paste(dirr,"wc_10m_HADCM3_A2a_", year, "_prec.tif", sep = ""))
  indicators = stack(epratio, bio, prec)
  indicators = crop(indicators, extent(africa))    
  indis[[year]] = rasterToPoints(indicators, spatial = TRUE)
  gridded(indis[[year]]) = TRUE
  proj4string(indis[[year]]) = CRS(proj4string(maml))
}
}




mamo = maml
maml = maml[maml$SHAPE_area < median(maml$SHAPE_area),]
spplot(pa, "sqkm", col.regions = bpy.colors(), lwd = 2, xlim = bbox(maml[,1]), ylim = bbox(maml[,2]),
   panel = function(x,y, ...){
    panel.polygonsplot(x,y, ...)
    sp.polygons(maml, col="red",fill=3, alpha = 0.02)
    })
   
#hrep = whri(pa, geometry = maml, indicators = indispdf) 

hrires = list()
whrires = list()
for (ifor in c("Current", "2020", "2050", "2080")) {
  hrires[[ifor]] = hri(habitat = pa, indicators = indis[["Current"]], 
             indicators2 = indis[[ifor]], forecast = TRUE)
  print(paste("done hres",pa$name, ip, ifor, attr(hrires[[ifor]], "hri")))
  whrires[[ifor]] = hri(habitat = maml, indicators = indis[["Current"]], 
             indicators2 = indis[[ifor]], forecast = TRUE)
  print(paste("done whres",pa$name, ip, ifor, attr(whrires[[ifor]], "hri")))
}


hres = indis[["Current"]]
whres = indis[["Current"]]
for (ifor in c("Current", "2020", "2050", "2080")) {
  hres@data = cbind(hres@data, hrires[[ifor]]$pHab)
  print(paste("done hres",ifor))
  whres@data = cbind(whres@data, whrires[[ifor]]$pHab)
  print(paste("done whres",ifor))
}

names(hres)[4:7] = c("Current", "y2020", "y2050", "y2080")
names(whres)[4:7] = c("Current", "y2020", "y2050", "y2080")
ares[[ip]] = list(hres = hres, whres = whres, hrires = hrires, whrires = whrires)



setwd("e:/text/semcon/spatialEcology/figs")
#ip = 6
hres = ares[[ip]]$hres
pdf(paste("hrinn",ip,".pdf", sep = ""))
print(spplot(hres,c( "y2050","y2080","Current", "y2020"), col.regions = bpy.colors(),
#    xlim = bbox(maml[,1]), ylim = bbox(maml[,2]),
   panel = function(x,y, ...){
    panel.gridplot(x,y, ...)
    sp.polygons(pa, col="red",fill=3, lwd = 3)
    }))
dev.off()



# ip = 6
whres = ares[[ip]]$whres
pdf(paste("whrinn",ip,".pdf", sep = ""))
print(spplot(whres, c( "y2050","y2080","Current", "y2020"), col.regions = bpy.colors(),
#    xlim = bbox(maml[,1]), ylim = bbox(maml[,2]),
   panel = function(x,y, ...){
    panel.gridplot(x,y, ...)
    sp.polygons(pa, col="red",fill=3, lwd = 3)
    }))
dev.off()     



for (ip in c(1:6,8:9)) {
  for (ifor in 1:4) {
    hres = ares[[ip]]$hrires[[ifor]]
    inHA = attr(hres, "inHA")
    whres = ares[[ip]]$whrires[[ifor]]
    hrepl = sum(hres$pHab > 0.5)/length(attr(hres, "inHA"))
    hwrepl = sum(whres$pHab > 0.5)/length(attr(hres, "inHA"))
    print(paste(ip, twpa$name[ip], length(attr(hres, "inHA")), 
          length(attr(whres, "inHA")),round(hrepl,3), round(hwrepl,3),
          round(mean(hres$pHab[inHA]), 3), round(mean(whres$pHab[inHA]),3),
          round(min(hres$pHab[inHA]), 3), round(min(whres$pHab[inHA]),3)))
  }
}
}