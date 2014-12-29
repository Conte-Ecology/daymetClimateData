# catchmentsShapefileDaymet - a shapefile of the catchments polygons in the Daymet projection

tileCatchmentsShapefile <- function(catchmentsShapefileDaymet){
  # This function only works for shapefiles with units that are lat/lon. Make sure the shapefile is in the daymet coordinate system before running this function on it.
  
  require(raster)
  require(maptools)
  
  # Get the bounds of the shapefile for tiling
  EXT <- extent(catchmentsShapefileDaymet)
  
  xmin <- floor(EXT@xmin)
  xmax <- ceiling(EXT@xmax)
  ymin <- floor(EXT@ymin)
  ymax <- ceiling(EXT@ymax)

  # Get the centroids of the polygons for tiling
  centroids <- data.frame(coordinates(catchmentsShapefileDaymet))
  names(centroids) <- c('X', 'Y')
  
  # Create tile bounds (by 2 degrees) for dividing catchments. Add 1 to ensure a last tile is not missed if there is an odd difference.
  Xs <- seq(from = xmin, to = xmax+1, by = 2)
  Ys <- seq(from = ymin, to = ymax+1, by = 2)
  
  # Get number of X and Y tiles for looping
  XL <- length(Xs) - 1
  YL <- length(Ys) - 1
  
  # Cerate storage for shapefiles
  catchmentsTiles <- list()
  
  
  tileNum <- 0
  # Loop through latitude (Y) bounds of tiles
  for( i in 1:YL ){
    
    # Y bounds
    YBmin <- Ys[i]
    YBmax <- Ys[i+1]
    
    # Loop through longitude (X) bounds of tiles
    for (j in 1:XL){
      
      # X bounds
      XBmin <- Xs[j]
      XBmax <- Xs[j+1]
  
      # Select the catchments that fall within the current tile boundaries
      rows <- which(centroids$X >= XBmin & centroids$X < XBmax & centroids$Y >= YBmin & centroids$Y < YBmax)
      
      # If catchments fall into the tile then store them as a separate shapefile
      if( length(rows) > 0 ){
        tileNum <- tileNum +1
        
        catchmentsTiles[[tileNum]] <- catchmentsShapefileDaymet[rows,]
      }
      
      # Print status
      print(paste0( "Tiling shapefile: ", j +(XL*(i-1)), " of ", XL*YL, " iterations complete.") )
    }
  }
  
  # Save memory
  rm(catchmentsShapefileDaymet);gc()

  return(catchmentsTiles)
}