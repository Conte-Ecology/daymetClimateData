spatialAverageFromMosaic <- function(daymetMosaicFilePath, spatialIndeces, zoneField){
  
  require(lubridate)
  require(ncdf)
  
  mosaicIndeces <- spatialIndeces$mosaicIndeces
  shapefileIndeces <- spatialIndeces$shapefileIndeces
  
  NetCDF <- nc_open( file.path(daymetMosaicFilePath) )
  
  # Get the year and variable directly from netcdf file
  year <- ncatt_get( NetCDF, varid = 0, attname="start_year")$value
  variable <- names(NetCDF$var)[names(NetCDF$var) %in% c("tmax", "tmin", "prcp", "dayl", "srad", "vp", "swe")]
    
  #Read "Day of Year" and make it start with 1
  dOY <- ncvar_get( nc = NetCDF, varid="yearday", start = 1, count = NetCDF$var$yearday$varsize  )
  dOY <- dOY + 1
  
  # Read variable
  shapeVar <- ncvar_get( nc   = NetCDF, 
                         varid = variable,  
                         start = c(mosaicIndeces$minRow[t], mosaicIndeces$minCol[t], 1), 
                         count = c(mosaicIndeces$countx[t], mosaicIndeces$county[t], length(dOY)) )
    
    
  # Index and replace missing value with NA
  # ---------------------------------------
  # Determine the currently assigned missing value
  for( h in 1:length(NetCDF$var) ){
    if ( NetCDF$var[[h]]$name == variable ) {varIndex <- h}
  }
  missingValue <- NetCDF$var[[varIndex]]$missval
    
  # Replace
  shapeVar <- replace(shapeVar, shapeVar == missingValue, NA)
    
  # Close the NetCDF connection
  nc_close(NetCDF)
    
  # Spatially average Daymet records for zones
  # ------------------------------------------
  # Setup storage dataframe (one year-long row for each point assigned to a zone)
  varPoints <- data.frame(matrix(nrow = nrow(shapefileIndeces), ncol = 366))
  names(varPoints) <- c(zoneField, dOY)
  varPoints[,zoneField] <- shapefileIndeces[,zoneField]
    
  # Loop through days, pulling records from variable array
  for (sel in dOY){
    selection <- as.matrix(data.frame(shapefileIndeces[c('subRow', 'subCol')], sel) )
    varPoints[,sel+1] <- shapeVar[selection]
  }
    
  # Add duplicate column Because dplyr struggles with variable column names
  varPoints$ZONE <- varPoints[,zoneField]
  varMeans <-  group_by(varPoints, ZONE) %>% 
    summarise_each(funs(mean))%>%
    group_by() %>%
    dplyr::select( -ZONE )
    
    
  # Melt means into output format
  varMelt <- melt(varMeans, id=c(zoneField))
  names(varMelt) <- c(zoneField, 'DayOfYear', variable)
  varMelt$DayOfYear <- as.numeric(varMelt$DayOfYear)
  varMelt$Year <- year
  varMelt$Date <- paste(parse_date_time(paste0(varMelt$Year,"-", varMelt$DayOfYear), "y-j", tz = "EST") )
  
  output <- varMelt[,c(zoneField, variable, "Date")]
  
  return(output)
}