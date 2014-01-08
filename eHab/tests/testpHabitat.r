library(eHab)
options(error = recover)
data(eHab)
summary(pHabitat(indicators,protectedArea))
#debug(predict.gstat)
system.time(hrep <- mhri(indicators,protectedArea, pval = 0.95, nsim = 4, nboot = 4))

#jar <- paste(system.file(package="dismo"), "/java/maxent.jar", sep='')
#  if (!file.exists(jar)) {
#	  system(paste("cp", "E:\\R\\Misc\\dismo\\dismoold\\inst\\java\\maxent.jar", 
#            paste(system.file(package="dismo"), "/java/", sep = "")))	
#  }

hrep = hri(indicators,protectedArea)
attr(hrep,"hri")
summary(hrep)
hrep = hri(indicators,protectedArea)
summary(hrep)
hrep = hri(indicators[,1],protectedArea)
summary(hrep)
	
attr(hrep,"mean")
attr(hrep,"cov")
s1 = system.time(hrep <- hri(indicators,protectedArea, pval = 0.95))
if (interactive()) s1
s2 = system.time(hrep2 <- hri(indicators,protectedArea, pval = 0.95, method = "maxent"))
if (interactive()) s2
s3 = system.time(hrep3 <- hri(indicators,protectedArea, pval = 0.95, method = "bioclim"))
if (interactive()) s3
s4 = system.time(hrep4 <- hri(indicators,protectedArea, pval = 0.95, method = "domain"))
if (interactive()) s4

hreps = hri(indicators,protectedArea, forecast = TRUE)


meanPoly = attr(hrep,"mean")
covPoly = attr(hrep,"cov")
hrep1 = hri(indicators2 = indicators,habitat = protectedArea, 
    covPoly = covPoly, meanPoly = meanPoly, pval = 0.95)
summary(hrep1)

polyData = attr(hrep,"polyData")
hrep2 = hri(indicators2 = indicators,habitat = protectedArea, polyData = polyData, pval = 0.95)
hrep3 = hri(indicators = indicators,habitat = protectedArea, polyData = polyData, pval = 0.95)
hrep4 = hri(indicators2 = indicators,habitat = protectedArea, polyData = polyData, pval = 0.95, forecast = TRUE)
all.equal(hrep@data, hrep1@data)
all.equal(hrep1@data, hrep2@data)
all.equal(hrep1@data, hrep3@data)
all.equal(hrep1@data, hrep4@data)


system.time(hrep <- mhri(indicators[,2:4],protectedArea, pval = 0.95, nsim = 4, nboot = 4))



#############################
if (FALSE) {
set.seed(1)
hrep = mhri(indicators,protectedArea, pval = 0.95, nsim = 4)
summary(hrep)
attr(hrep,"mean")
attr(hrep,"cov")



findPatches(hrep)

}


system.time(hrep <- mhri(indicators,protectedArea, pval = 0.95, nsim = 4, nboot = 4))
