determineSpatialRelationships <- function(catchmentsShapefile, zoneField, daymetDirectory, variables, years){
    
  # Tiles large catchments into smaller spatial objects for processing
  # Compares shapefile with spatial daymet information
  # Determines the Daymet mosaic bounds for the shapefile
  # Determines the array positions in the mosaic subset for each polygon
    
  # Packages used
  require(dplyr)
  require(ncdf4) 
  require(rgdal)
  
  source('C:/KPONEIL/GitHub/projects/daymetClimateData/scripts/functions/tileCatchmentsShapefile.R')
  
  # Transform the shapefile into the coordinate system so the units are in lat/lon
  # This makes the shapefile comparable to the coordinates provided by Daymet in WGS
  transformShapefile <- spTransform(catchmentsShapefile, 
                                    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"), 
                                    class = "SpatialPolygonsDataFrame")
  
  # Tiles the catchment shapefile into 2 degree sections to prevent memory issues when pulling netCDF data
  tiledShapefile <- tileCatchmentsShapefile(transformShapefile)
  
  # Connect with the sample netCDF to 
  spatialNCDF <- nc_open(file.path(daymetDirectory, paste0( variables[1], '_', years[1], '.nc4') ) )  #netcdf     
  
  # Latitude coordinates in full mosaic
  mosaicLat = ncvar_get ( nc    = spatialNCDF, 
                          varid = "lat", 
                          start = c(1,1), 
                          count = c(spatialNCDF$var$lat$varsize[1], 
                                    spatialNCDF$var$lat$varsize[2]) )
  
  # Longitude coordinates in full mosaic
  mosaicLon = ncvar_get ( nc    = spatialNCDF, 
                          varid = "lon", 
                          start = c(1,1), 
                          count = c(spatialNCDF$var$lon$varsize[1], 
                                    spatialNCDF$var$lon$varsize[2]) )
  
  # Output objects
  mosaicIndeces <- list() # Change to ncdfIndeces
  tileIndeces  <- list() # Change to tileIndeces
  
  
  # Loop through tiles to determine the index information for the tiles
  for(t in 1:length(tiledShapefile) ){
    
    tile <- tiledShapefile[[t]]
    
    # Print status
    print(paste("Beginning spatial analysis of tile", t, "of", length(tiledShapefile), ", which has", nrow(tile), "catchments."))
    
    # Calculate the indeces of the mosaic netcdf for the tile
    # -------------------------------------------------------
    # Get the extent of the shapefile
    tileExtent <- extent(tile)
    
    # Determine the array positions of coordinates within the shapefile extent for spatial referencing
    arrayIndeces <- which( mosaicLat >= tileExtent@ymin & 
                             mosaicLat <= tileExtent@ymax & 
                             mosaicLon >= tileExtent@xmin & 
                             mosaicLon <= tileExtent@xmax, 
                           arr.ind = T)
    
    # Corners of the box in the array
    minRow <- min(arrayIndeces[,1])
    maxRow <- max(arrayIndeces[,1])
    minCol <- min(arrayIndeces[,2])
    maxCol <- max(arrayIndeces[,2])
    
    # Number of rows and columns
    countx = maxRow - minRow + 1 
    county = maxCol - minCol + 1 
    
    # Save the spatial references for use in other netCDF files
    mosaicIndeces[[t]] <- data.frame(tile = t, 
                                     minRow = minRow, 
                                     maxRow = maxRow, 
                                     minCol = minCol, 
                                     maxCol = maxCol, 
                                     countx = countx, 
                                     county = county)
    
    # Create a spatial object of the tile Daymet points
    # -------------------------------------------------
    # Latitude coordinates in the current tile
    tileLat = ncvar_get( nc    = spatialNCDF, 
                         varid = "lat",     
                         start = c(minRow, minCol),    
                         count = c(countx,county) )
    
    # Longitude coordinates in the current tile
    tileLon = ncvar_get( nc    = spatialNCDF, 
                         varid = "lon",     
                         start = c(minRow, minCol),    
                         count = c(countx,county) )
    
    # Close connection
    nc_close(spatialNCDF)
    
    # Convert coordinate to a spatial object
    tilePoints <- as.data.frame( cbind( as.vector(tileLon), as.vector(tileLat) ) )
    colnames(tilePoints) <- c("Longitude", "Latitude")
    tilePointsSpPts   <- SpatialPoints(tilePoints, proj4string = CRS(proj4string(transformShapefile)))
    
    # Daymet points get assigned to polygons using 2 methods:
    # 1. A spatial join is used to determine all points fall inside of each polygon. These will eventually get spatially averaged
    # 2. Polygons that are too small to overlap any points get assigned the point that is nearest to their centroid
    
    
    # Spatial Join
    # ------------
    # Overlay points on catchment shapefile
    overPoints <- over(tilePointsSpPts, tile)
    
    spatialJoin <- cbind(overPoints, tilePoints)
    
    polygonPoints <- spatialJoin[which(!is.na(spatialJoin[,zoneField])),]
    
    # Nearest Point
    # -------------
    # Determine the polygons without assigned points
    missingPolygons <- tile[!tile@data[,zoneField] %in% polygonPoints[,zoneField],]
    
    missingCentroids <- data.frame(missingPolygons@data[,zoneField], coordinates(missingPolygons) )
    names(missingCentroids) <- c(zoneField, 'LON', 'LAT')
    
    missingIDs <- missingPolygons@data[,zoneField]
    
    nearPoints <- as.data.frame(matrix(nrow = length(missingIDs), ncol = 3))
    names(nearPoints) <- c(zoneField, 'Longitude', 'Latitude')
    
    # Assign the nearest point to polygons missing point assignments
    for ( i in seq_along(missingIDs) ) {
      
      tempLat <- missingCentroids$LAT[missingCentroids[,zoneField] == missingIDs[i]]
      tempLon <- missingCentroids$LON[missingCentroids[,zoneField] == missingIDs[i]]
      
      distances <- spDistsN1(as.matrix(tilePoints), c(tempLon, tempLat), longlat = TRUE)
      minDist <- min(distances)
      distPos <- which(distances == minDist)[1]
      
      nearPoints[i ,zoneField] <- missingIDs[i] 
      nearPoints$Longitude[i]  <- tilePoints[distPos, 1]
      nearPoints$Latitude[i]   <- tilePoints[distPos, 2]
    }
    
    # Join the two versions of point assignments
    finalPoints <- rbind(polygonPoints[,c(zoneField, 'Longitude', 'Latitude')], nearPoints)
    
    
    # Determine the position in the tile array
    # ----------------------------------------
    finalPoints$subRow <- NA; finalPoints$subCol <- NA
    
    for ( m in 1:nrow(finalPoints) ){
      
      # Find th position in the array of the variable
      position <- which(tileLon == finalPoints$Longitude[m] & tileLat == finalPoints$Latitude[m], arr.in = TRUE)
      
      finalPoints$subRow[m] <- as.numeric(position[,1])
      finalPoints$subCol[m] <- as.numeric(position[,2])
    }
    
    tileIndeces[[t]] <- finalPoints
    names(tileIndeces)[t] <- paste0("Tile", t)
    
    print(paste("Tile", t, "complete."))
  }
  
  # Prepare output
  # --------------
  mosaicIndeces <- as.data.frame(rbind_all(mosaicIndeces))
  
  output <- list(mosaicIndeces, tileIndeces)
  names(output) <- c('mosaicIndeces', 'tileIndeces')
  
  return(output)
}