
KappaV <-
function(shp1.path, shp2.path,
  shp1.fieldID = "ID", shp2.fieldID = shp1.fieldID,
  shp1.fieldOS = "OS", shp2.fieldOS = shp1.fieldOS, plot = FALSE) {
  shp1       = readShapeSpatial(shp1.path)
  shp2       = readShapeSpatial(shp2.path)
  shp1.fieldID2 <- paste(shp1.fieldID, "2", sep="")
  shp2.fieldID2 <-paste(shp1.fieldID, "2", sep="")
  ov <- over(shp1, shp2, returnList=T)

  shp1.lev <- sort(unique(shp1@data[, shp1.fieldOS]))
  shp2.lev <- sort(unique(shp2@data[, shp2.fieldOS]))

  if (plot) {
    layout(matrix(1:2, ncol=2))
    plot(shp1, col=topo.colors(length(unique(shp1@data[, shp1.fieldOS])))[shp1@data[, shp1.fieldOS]+1])
    plot(shp2, col=topo.colors(length(unique(shp2@data[, shp1.fieldOS])))[shp2@data[, shp2.fieldOS]+1])
    layout(matrix(1))
  }
  
  res <- matrix(0, nrow=length(shp1.lev), ncol=length(shp2.lev),
                dimnames=list(shp1.lev, shp2.lev))

  n1 <- length(shp1@polygons)
  n2 <- length(shp2@polygons)

  shp1@data[, shp1.fieldID2] <- 1:n1
  shp2@data[, shp2.fieldID2] <- 1:n2

  for (i in 1:n1){
    xi <- SpatialPolygons(list(shp1@polygons[[i]]))
    pi <- as(xi, "gpc.poly")
    #cat("\n", i, " ")
    for (j in seq(along=ov[[i]])) {
      yj     <- SpatialPolygons(list(shp2@polygons[[ov[[i]][j]]]))
      pj     <- as(yj, "gpc.poly")    
      int.ij <- area.poly(intersect(pi, pj))
      #int.ij <- gpclib::area.poly(gpclib::intersect(pi, pj)) # previously used with gpclib
      if (int.ij != 0) {
        #cat("*")
        ri <- which(shp1.lev == shp1@data[shp1@data[ ,shp1.fieldID2]==i, shp1.fieldOS])
        ci <- which(shp2.lev == shp2@data[shp2@data[ ,shp2.fieldID2]==ov[[i]][j], shp2.fieldOS])
        res[ri, ci] <- res[ri, ci] + int.ij
      #} else {(cat("."))}
      } else {}
    }
  }
  #cat("\n")
  return(list(confusion.matrix=res, kappa.v=Kappa(res)))
}
