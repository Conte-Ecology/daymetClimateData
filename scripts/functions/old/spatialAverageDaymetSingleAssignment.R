# NetCDF - The NetCDF connection to the file
# variable - character string, the name of the variable
# catchmentsShapefile - a shapefile of the catchments polygon


spatialAverageDaymet <- function(NetCDF, variable, catchmentsShapefile){
  
  require(raster)
  require(ncdf4)
  require(sp)
  require(rgdal)
  require(tcltk)
  
  
  catchmentsShapefile <- tileList[[1]]
  
  NetCDF <- nc_open('F:/KPONEIL/SourceData/climate/DAYMET/raw/prcp_1980.nc4')    #netcdf     
  
  
  # Index the portion of the NetCDF covered by the shapefile
  # --------------------------------------------------------
  # Read in variables
  lat = ncvar_get ( nc=NetCDF, varid="lat", start = c(1,1), count = c(NetCDF$var$lat$varsize[1], NetCDF$var$lat$varsize[2]) )
  lon = ncvar_get ( nc=NetCDF, varid="lon", start = c(1,1), count = c(NetCDF$var$lon$varsize[1], NetCDF$var$lon$varsize[2]) )
  
  ## Get the extent of the shapefile
  EXT <- extent(catchmentsShapefile)
  
  # Positions in the array of coordinates within the shapefile extent
  matInd <- which(lat >= EXT@ymin & lat <= EXT@ymax & lon >= EXT@xmin & lon <= EXT@xmax, arr.ind = T)
  
  # Corners of the box in the array
  minRow <- min(matInd[,1])
  maxRow <- max(matInd[,1])
  minCol <- min(matInd[,2])
  maxCol <- max(matInd[,2])
  
  # Number of rows and columns
  countx = maxRow - minRow + 1 
  county = maxCol - minCol + 1 
  
  # Remove the full NetCDF lat/lon coordinates. These will be replaced with the indexed set of coords
  rm(lat,lon)
  
  lat = ncvar_get( nc = NetCDF, varid="lat",     start = c(minRow, minCol),    count = c(countx,county) )
  lon = ncvar_get( nc = NetCDF, varid="lon",     start = c(minRow, minCol),    count = c(countx,county) )
  
  
  
  
  
  
  masterCoords <- cbind( as.vector(lon), as.vector(lat))
  colnames(masterCoords) <- c("Longitude", "Latitude")
  masterCoords       <- as.data.frame(masterCoords)
  masterCoordsMatrix <- as.matrix(masterCoords)
  masterCoordSpPts   <- SpatialPoints(masterCoords, proj4string = CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
  
  # Generate list of catchment centroids
  centroids <- data.frame(catchmentsShapefile@data$FEATUREID, coordinates(catchmentsShapefile) )
  names(centroids) <- c('FEATUREID', 'LON', 'LAT')
  
  # Overlay points on catchment shapefile
  overPoints <- over(masterCoordSpPts, catchmentsShapefile)
  
  #save(overPoints, file = 'C:/KPONEIL/workspace/overPointsAllTile1.RData')
  
  
  # Generate a list of the points and which 
  pointsInside <- overPoints[which(!is.na(overPoints$FEATUREID)),]
  
  # Add the coordinates of the points inside the catchments
  ptsInCoords <- as.data.frame(masterCoordSpPts[as.numeric(row.names(pointsInside))])
  
  pointsInside$Longitude <- ptsInCoords$x
  pointsInside$Latitude  <- ptsInCoords$y
  
  #######
  for ( i in seq_along(fids) ){
    
    setTkProgressBar(progressBar, i, label=paste("Spatial averaging ", round(i/length(fids)*100, 0), "% complete"))
    
    inside <- pointsInside[which(pointsInside$FEATUREID == fids[i]),c('Longitude', 'Latitude')]
    
    #If no point falls within the catchment, find the nearest one:
    #-------------------------------------------------------------
    if(nrow(inside) == 0 ){
      
      tempLat <- centroids$LAT[centroids$FEATUREID == fids[i]]
      tempLon <- centroids$LON[centroids$FEATUREID == fids[i]]
      
      distances <- spDistsN1(masterCoordsMatrix, c(tempLon, tempLat), longlat = TRUE)
      minDist <- min(distances)
      distpos <- which(distances == minDist)[1]
      
      nearLon  <- masterCoords[distpos, 1]
      nearLat  <- masterCoords[distpos, 2]
      
      inside[1,1] <- nearLon
      inside[1,2] <- nearLat
    }
    
    # Index data for points within the watershed
    # ------------------------------------------
    for ( m in 1:nrow(inside) ){
      
      # Find th position in the array of the variable
      position <- which(lon == inside$Longitude[m] & lat == inside$Latitude[m], arr.in = TRUE)
      # Pull those values
      indVar <- var[as.numeric(position[,1]), as.numeric(position[,2]), 1:365]
      # Join them in an array
      if (m == 1) {tempVar <- data.frame(dOY, indVar)} else(tempVar <- cbind(tempVar, indVar))
    }
    
    # Which columns to average
    toBeAvg <- which(names(tempVar) == 'indVar')
    # Take the average across all points
    ifelse( length(toBeAvg) > 1, R <- rowMeans(tempVar[,toBeAvg], na.rm = TRUE, dims = 1),  R <- tempVar[,2] )
    
    # Add means to data.frame
    #varMeans[,paste(fids[i])] <- R
    
    
    dfRows <- c( (1 + length(dOY)*(i-1)) : (length(dOY) + length(dOY)*(i-1)) )
    
    varMeans$FEATUREID[dfRows] <- fids[i]
    varMeans[dfRows, variable] <- R 
  }
  close(progressBar)
  
  
  
  ################################################################################
  
  ## Get the extent of the shapefile
  EXT <- extent(catchmentsShapefile)
  
  # Positions in the array of coordinates within the shapefile extent
  matInd <- which(lat >= EXT@ymin & lat <= EXT@ymax & lon >= EXT@xmin & lon <= EXT@xmax, arr.ind = T)
  
  # Corners of the box in the array
  minRow <- min(matInd[,1])
  maxRow <- max(matInd[,1])
  minCol <- min(matInd[,2])
  maxCol <- max(matInd[,2])
  
  # Number of rows and columns
  countx = maxRow - minRow + 1 
  county = maxCol - minCol + 1 
  
  # Remove the full NetCDF lat/lon coordinates. These will be replaced with the indexed set of coords
  rm(lat,lon)
  
  # Read the variables for the subsetted netcdf
  #print("Reading NetCDF values...")
  lat = ncvar_get( nc = NetCDF, varid="lat",     start = c(minRow, minCol),    count = c(countx,county) )
  lon = ncvar_get( nc = NetCDF, varid="lon",     start = c(minRow, minCol),    count = c(countx,county) )
  dOY = ncvar_get( nc = NetCDF, varid="yearday", start = 1,                    count = NetCDF$var$yearday$varsize  )
  var = ncvar_get( nc = NetCDF, varid=variable,  start = c(minRow, minCol, 1), count = c(countx,county,length(dOY)) )
  
  # Correction for Daymet doy which starts at 0.
  dOY <- dOY + 1 
  
  # Get the year for output
  YEAR <- ncatt_get( NetCDF, varid = 0, attname="start_year")$value
  
  # Index and replace missval
  # -------------------------
  for( h in 1:length(NetCDF$var) ){
    if ( NetCDF$var[[h]]$name == variable ) {varIndex <- h}
  }
  missingValue <- NetCDF$var[[varIndex]]$missval
  
  # Replace
  var <- replace(var, var == missingValue, NA)
  
  # Prep the points for indexing
  # ----------------------------
  
  # Join coordinate lists
  masterCoords <- cbind( as.vector(lon), as.vector(lat))
  colnames(masterCoords) <- c("Longitude", "Latitude")
  masterCoords       <- as.data.frame(masterCoords)
  masterCoordsMatrix <- as.matrix(masterCoords)
  masterCoordSpPts   <- SpatialPoints(masterCoords, proj4string = CRS(proj4string(catchmentsShapefile)))

  
  # Generate list of catchment centroids
  centroids <- data.frame(catchmentsShapefile@data$FEATUREID, coordinates(catchmentsShapefile) )
  names(centroids) <- c('FEATUREID', 'LON', 'LAT')
  
  # Overlay points on catchment shapefile
  overPoints <- over(masterCoordSpPts, catchmentsShapefile)
  
  # Generate a list of the points and which 
  pointsInside <- overPoints[which(!is.na(overPoints$FEATUREID)),]
  
  # Add the coordinates of the points inside the catchments
  ptsInCoords <- as.data.frame(masterCoordSpPts[as.numeric(row.names(pointsInside))])
  
  pointsInside$Longitude <- ptsInCoords$x
  pointsInside$Latitude  <- ptsInCoords$y
  
  # List of FeatureIDs to average over
  fids <- catchmentsShapefile@data$FEATUREID
  
  # Replace with means
  varMeans <- as.data.frame(matrix(nrow = length(fids)*length(dOY), ncol = 4))
  names(varMeans) <- c('FEATUREID', 'Year', 'DayOfYear', variable)
  
  varMeans$Year <- YEAR
  varMeans$DayOfYear <- rep(dOY, length(fids))
  
  progressBar <- tkProgressBar(title = "Progress Bar", min = 0, max = length(fids), width = 300)
  
  for ( i in seq_along(fids) ){
    
    setTkProgressBar(progressBar, i, label=paste("Spatial averaging ", round(i/length(fids)*100, 0), "% complete"))
    
    inside <- pointsInside[which(pointsInside$FEATUREID == fids[i]),c('Longitude', 'Latitude')]
    
    #If no point falls within the catchment, find the nearest one:
    #-------------------------------------------------------------
    if(nrow(inside) == 0 ){
      
      tempLat <- centroids$LAT[centroids$FEATUREID == fids[i]]
      tempLon <- centroids$LON[centroids$FEATUREID == fids[i]]
      
      distances <- spDistsN1(masterCoordsMatrix, c(tempLon, tempLat), longlat = TRUE)
      minDist <- min(distances)
      distpos <- which(distances == minDist)[1]
      
      nearLon  <- masterCoords[distpos, 1]
      nearLat  <- masterCoords[distpos, 2]
      
      inside[1,1] <- nearLon
      inside[1,2] <- nearLat
    }
    
    # Index data for points within the watershed
    # ------------------------------------------
    for ( m in 1:nrow(inside) ){
      
      # Find th position in the array of the variable
      position <- which(lon == inside$Longitude[m] & lat == inside$Latitude[m], arr.in = TRUE)
      # Pull those values
      indVar <- var[as.numeric(position[,1]), as.numeric(position[,2]), 1:365]
      # Join them in an array
      if (m == 1) {tempVar <- data.frame(dOY, indVar)} else(tempVar <- cbind(tempVar, indVar))
    }
    
    # Which columns to average
    toBeAvg <- which(names(tempVar) == 'indVar')
    # Take the average across all points
    ifelse( length(toBeAvg) > 1, R <- rowMeans(tempVar[,toBeAvg], na.rm = TRUE, dims = 1),  R <- tempVar[,2] )
    
    # Add means to data.frame
    #varMeans[,paste(fids[i])] <- R
    
    
    dfRows <- c( (1 + length(dOY)*(i-1)) : (length(dOY) + length(dOY)*(i-1)) )
    
    varMeans$FEATUREID[dfRows] <- fids[i]
    varMeans[dfRows, variable] <- R 
  }
  close(progressBar)
  
  # Return the records
  return(varMeans)
}




