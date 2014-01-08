mecohri = function(ecoreg, ecoregions, ecoregs, ecoID,...) {
  removeTmpFiles(h=12)
  write(paste("Entering ecoreg", ecoreg, "with", Sys.getpid()), file = "mecohri.txt", append = TRUE)
  eids <- which(ecoregions@data[,ecoID] %in% ecoreg) # changed @javier
  if (length(eids) > 0) {
    leco <- ecoregions[eids,] # changed @javier
    print(dim(leco))
    ress <- try(ecohri(leco, ecoID = ecoID, ...)) # name changed to ress by @javier (26.07.13)
    print(ress) # added by Javier (26.07.13)
    if (is(ress, "try-error")) ress = -999
  } else ress = NULL
  ids = which(ecoreg %in% ecoregs)
  if (!is.null(ress) && length(ress)>1 && dim(ress$hriRes)[1]!=0) { # changed by Javier (26.07.13) and updated (29.07.13) # ress!=-999
    print(paste("**** mecohri ", ecoreg, "tiles", attr(ress,"tr")$n, "parks", sum(!is.na(ress$hriRes[,1])), "****"))
    write(paste("Done ecoreg", ecoreg, "with", Sys.getpid(), 
                "tiles", attr(ress,"tr")$n, "parks", sum(!is.na(ress$hriRes[,1]))), file = "mecohri2.txt", append = TRUE)

	write.table(ress$hriRes, file = "hriRes.csv",sep=',', append=TRUE, row.names=FALSE,col.names = ifelse(exists("hriRes.csv"), FALSE, TRUE))
	write.table(ress$hriRes2, file = "hriRes2.csv",sep=',',append=TRUE, row.names=FALSE,col.names = ifelse(exists("hriRes2.csv"), FALSE, TRUE))
	write.table(ress$hriInRes, file = "hriInRes.csv",sep=',',append=TRUE, row.names=FALSE,col.names = ifelse(exists("hriInRes.csv"), FALSE, TRUE))
	write.table(ress$hriInRes2, file = "hriInRes2.csv",sep=',',append=TRUE, row.names=FALSE,col.names = ifelse(exists("hriInRes2.csv"), FALSE, TRUE))
	write.table(ress$errors, file = "errors.csv",sep=',',append=TRUE, row.names=FALSE)

  } else {
    print(paste("**** mecohri ", ecoreg, "NULL" , "****"))
    write(paste("Done ecoreg", ecoreg, "with", Sys.getpid(), "NULL"), file = "mecohri2.txt", append = TRUE)
  }  
  ress
}