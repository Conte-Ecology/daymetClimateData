convertDaymetToSQLiteDatabase <- function(databaseConnection, tableName, years, variables, daymetDirectory, catchmentsShapefile, daymetProj4string){
  
  # Start time
  S <- proc.time()[3]
  
  require(RSQLite)
  require(RSQLite.extfuns)
  require(dplyr)
  require(lubridate)
  
  catchmentsProj <- transformShapefile(shapefile = catchmentsShapefile, 
                                       newProjection = daymetProj4string)
  
  # Tile the catchments shapefile if it is large
  tileList <- tileCatchmentsShapefile(catchmentsProj); 
  
  # Print status
  print( "Shapefile tiling complete, begining spatial averaging.")
  
  # Save memory
  rm(catchmentsShapefile);rm(catchmentsProj);gc()
  
  for ( k in 1:length(tileList)){
    
    subCatchments <- tileList[[k]]
    
    for (year in years){
      
      annualRecord <- NULL
      for (var in variables){
        
        # Open the NetCDF file
        NCDF <- nc_open(file.path(daymetDirectory, paste0(var, '_',year , '.nc4')))    #netcdf           
        
        currentRecord <- spatialAverageDaymet(NetCDF = NCDF,
                                              variable = var,
                                              catchmentsShapefile = subCatchments)
        
        nc_close(NCDF)
        
        # Store the variables records for the current year
        if( is.null(annualRecord)){
          annualRecord <- currentRecord
        } else( annualRecord <- left_join(annualRecord, currentRecord, by = c('FEATUREID', 'Year', 'DayOfYear')))
        print( paste0(var, " complete for ", year, " in tile #", k, " of ", length(tileList), ".") )
      }# End variables loop
      
      # Create a "Date" column
      annualRecord$Date <- paste(parse_date_time(paste0(annualRecord$Year,"-", annualRecord$DayOfYear), "y-j", tz = "EST") )
      
      # Return the desired columns 
      annualRecord <- annualRecord[,c('FEATUREID', 'Date', variables)]
      
      
      # Status update with time elapsed
      print(paste0("Writing results into database for ", year, " in tile #", k, " of ", length(tileList), ". Elapsed time: ", (proc.time()[3] - S)/3600, " hours."))
      
      # If the table doesn't exist in the database, create it. If if does, append to it.
      if(!tableName %in% dbListTables(databaseConnection) ){
        dbWriteTable(conn = databaseConnection, 
                     name = tableName, 
                     value = annualRecord, 
                     append = FALSE, 
                     row.names = FALSE)
      } else( dbWriteTable(conn = databaseConnection, 
                           name = tableName, 
                           value = annualRecord, 
                           append = TRUE,  
                           row.names = FALSE) )
      
    }# End year loop
  }# End variable loop
}# End function