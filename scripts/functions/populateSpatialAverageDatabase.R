populateSpatialAverageDatabase <- function(zonesShapefile, zoneField, daymetDirectory, variables, years, databaseFilePath, tableName){
  
  # Load required packages
  require(dplyr)
  require(RSQLite)
  require(RSQLite.extfuns)
  
  # Load functions
  source('C:/KPONEIL/GitHub/projects/daymetClimateData/scripts/functions/tileShapefile.R')
  source('C:/KPONEIL/GitHub/projects/daymetClimateData/scripts/functions/determineSpatialRelationshipsSingleShapefile.R')
  source('C:/KPONEIL/GitHub/projects/daymetClimateData/scripts/functions/spatialAverageFromMosaic.R')
  
  # Create and connect to the database
  # ----------------------------------
  # If the database does not exist then create it, otherwise establish connection with existing database for uploading.
  if ( !file.exists( databaseFilePath ) ) { 
    src_sqlite(databaseFilePath, create = T) 
  } else ( print("Database already exists. Processed records will be added to existing database."))
  
  database <- dbConnect(SQLite(), databaseFilePath)
  
  # Spatial processing
  # ------------------
  # Checks shapefile projection
  if (proj4string(zonesShapefile) != "+proj=lcc +ellps=WGS84 +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0"){
    stop("Error: The correct projection has not been defined for the shapefile. Shapefile should be read in as Lambert (+proj=lcc +ellps=WGS84 +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0).")
  }
  
  # Transform the shapefile into the coordinate system so the units are in lat/lon. This makes the shapefile comparable to the coordinates provided by Daymet NetCDFs in WGS.
  transformShapefile <- spTransform(zonesShapefile,
                                    CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"),
                                    class = "SpatialPolygonsDataFrame")
  
  # If the shapefile is larger than 2 degrees by 2 degrees then it gets tiled.
  EXT <- extent(transformShapefile)
  
  # If the shapefile fits into one tile then skip tiling, else tile shapefile.
  if( abs(EXT@xmin - EXT@xmax) <= 2 & abs(EXT@ymax - EXT@ymin) <= 2){
    print("Shapefile is not large enough to warrant tiling. Original structure is preserved.")
    tiledShapefile <- list(transformShapefile)
  }else{
    # Tiles the shapefile into 2 degree sections to prevent memory issues when pulling netCDF data
    tiledShapefile <- tileShapefile(shapefile = transformShapefile) 
  }
  
  # Record averaging by tile, year, and variable
  # --------------------------------------------
  for( t in length(tiledShapefile) ){
    
    # Subsets the Daymet mosaic netcdf file based on the provided shapefile
    spatialIndeces <- determineSpatialRelationships(zonesShapefile = tiledShapefile[[t]],
                                                    zoneField = zoneField,
                                                    exampleDaymetFile = file.path(daymetDirectory, paste0(variables[1], '_', years[1], '.nc4') ) )
    for( year in years ){
      
      S <- proc.time()[3]
      
      annualRecord <- NULL
      
      for( variable in variables ){
          
        print(paste0("Spatially averaging '", variable, "' records in ", year, " for shapefile ", t, " of ", length(tiledShapefile), "." ))
        
        ncdfPath <- file.path(daymetDirectory, paste0(variable, '_', year, '.nc4') )
           
        record <- spatialAverageFromMosaic(daymetMosaicFilePath = ncdfPath, 
                                           spatialIndeces = spatialIndeces, 
                                           zoneField = zoneField)  
        
        # Store the variables records for the current year
        if(is.null(annualRecord)){
          annualRecord <- record
        } else( annualRecord <- left_join(annualRecord, record, by = c(zoneField, 'Date')))
      }# End variables loop
      
      # Upload to database
      # ------------------
      # If the table doesn't exist in the database, create it. If if does, append to it.
      if(!tableName %in% dbListTables(database) ){
        dbWriteTable(conn = database, 
                     name = tableName, 
                     value = annualRecord, 
                     append = FALSE, 
                     row.names = FALSE)
      } else( dbWriteTable(conn = database, 
                           name = tableName, 
                           value = annualRecord, 
                           append = TRUE,  
                           row.names = FALSE) )
      
      # Status update with time elapsed
      print(paste0("Done writing results into database for ", year, " in tile #", t, " of ", length(tiledShapefile), ". Elapsed time: ", (proc.time()[3] - S)/60, " minutes."))
      
    }# End years loop
  }# End tile loop
}# End function
