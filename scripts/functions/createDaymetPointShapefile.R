# This function generates a SpatialPoints object of all of the points that fall within the extent of a given shapefile

# boundaryShapefile - a SpatialPolygonsDataFrame object of the boundary around the Daymet points
# NetCDF - the netcdf to access
# includeVariable - TRUE/FALSE of whether or not to include the variable values in the attribute table
# dOY - if includeVariable is TRUE, the day of year specified to pull the variable from



createDaymetPointShapefile <- function(boundaryShapefile, NetCDF, includeVariable, dOY){
  
  require(ncdf4)
  require(raster)
  require(sp)
    
  if( proj4string(boundaryShapefile) != "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"){
    stop("Error in 'createDaymetPointShapefile' function: The spatial relationship between the shapefile and Daymet grid cannot be determined. The shapefile must be in WGS.")
  }
  
  # Read in variables
  lat = ncvar_get ( nc=NetCDF, varid="lat", start = c(1,1), count = c(NetCDF$var$lat$varsize[1], NetCDF$var$lat$varsize[2]) )
  lon = ncvar_get ( nc=NetCDF, varid="lon", start = c(1,1), count = c(NetCDF$var$lon$varsize[1], NetCDF$var$lon$varsize[2]) )
    
  # Convert NetCDF info to spatial data for overlay
  # -----------------------------------------------
  # Get the extent of the shapefile
  EXT <- extent(boundaryShapefile)
  
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
  
  # Number of rows and columns (Minimum of 2 to make sure the variable is read in as an array for uniformity.)
  countX = max(maxRow - minRow + 1, 2)
  countY = max(maxCol - minCol + 1, 2)
  
  # Remove the full NetCDF lat/lon coordinates. These will be replaced with the indexed set of coords
  rm(lat,lon)
  
  # Read the variables for the subsetted netcdf
  lat = ncvar_get( nc = NetCDF, varid="lat",     start = c(minRow, minCol),    count = c(countX,countY) )
  lon = ncvar_get( nc = NetCDF, varid="lon",     start = c(minRow, minCol),    count = c(countX,countY) )
  if( includeVariable == T ) {
    
    # Get variable name
    variable <- names(NetCDF$var)[names(NetCDF$var) %in% c("tmax", "tmin", "prcp", "dayl", "srad", "vp", "swe")]
    
    var = ncvar_get( nc = NetCDF, varid=variable,  start = c(minRow, minCol, dOY), count = c(countX,countY,1) )
    
    # Index and replace missval
    # -------------------------
    for( h in 1:length(NetCDF$var) ){
      if ( NetCDF$var[[h]]$name == variable ) {varIndex <- h}
    }
    missingValue <- NetCDF$var[[varIndex]]$missval
    
    # Replace with NA (for averaging purposes)
    var <- replace(var, var == missingValue, -9999)
    var <- replace(var, is.na(var),          -9999)
    
  }
  
  daymetCoordinates <- cbind( as.vector(lon), as.vector(lat))
  colnames(daymetCoordinates) <- c("Longitude", "Latitude")
  daymetCoordinates       <- as.data.frame(daymetCoordinates)

  # Make it a spatial points object
  # ----------------------------------
  if( includeVariable == T ){
    daymetCoordinates$variable <- as.vector(var)
    daymetCoordinatesSpPts   <- SpatialPointsDataFrame(daymetCoordinates[,1:2],  daymetCoordinates, proj4string = CRS(proj4string(boundaryShapefile))) 
  } else( daymetCoordinatesSpPts   <- SpatialPoints(daymetCoordinates,  proj4string = CRS(proj4string(boundaryShapefile))))

  
  # Return the point spatial object
  return(daymetCoordinatesSpPts)
}