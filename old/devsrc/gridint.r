gridint = function(spgrid, newpoints, unc = FALSE, model) {
  if (unc) {
    return
  
  
  
  } else {
    gridpar = gridparameters(spgrid)
    dx = gridpar$cellsize[1]
    dy = gridpar$cellsize[2]
    nx = gridpar$cells.dim[1]
    ny = gridpar$cells.dim[2]
    x0 = gridpar$cellcentre.offset[1]
    y0 = gridpar$cellcentre.offset[2]
    dat = as.matrix(spgrid)
    
    nnew = dim(newpoints)[1]
    pNew = cbind(newpoints,rep(0,dim(newpoints)[1]))
#    nx ,ny, nnew, pNew, x0, y0, dx, dy, dat
    bilinres = .Fortran("bilin", nx, ny, nnew, pNew, x0, y0, dx, dy, dat)  
#    bilinres = .Fortran("bilin", nx ,ny, nnew, pNew, x0, y0, dx, dy, dat)  
    bilinres[[4]]
}
}