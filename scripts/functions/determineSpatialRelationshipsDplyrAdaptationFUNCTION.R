
rm(list = ls())

require(sp)
require(ncdf4)
require(maptools)
require(raster)
require(rgdal)
require(lubridate)
library(RSQLite)
library(RSQLite.extfuns)
library(dplyr)
library(reshape2)

#source('C:/KPONEIL/GitHub/personal/dataProcessing/scripts/functions/tileCatchmentsShapefile.R')
#source('C:/KPONEIL/GitHub/personal/dataProcessing/scripts/functions/determineSpatialRelationships.R')



proj4.Lambert <- "+proj=lcc +ellps=WGS84 +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0"




# Temporal range
START_YEAR <- 1980
END_YEAR <- 2013

YEARS <- seq(from = START_YEAR, to = END_YEAR, by = 1)

# Variables
VARIABLES <- c('tmax', 'tmin', 'prcp', 'dayl', 'srad', 'vp', 'swe')

DAYMET_DIRECTORY <- 'F:/KPONEIL/SourceData/climate/DAYMET/raw'

tableName <- "climateRecord"

CATCHMENTS <- readShapePoly('C:/KPONEIL/GitHub/projects/daymetClimateData/gisFiles/shapefiles/NENYHRD_AllCatchments_DaymetPrj.shp', proj4string=CRS(proj4.Lambert))




DATABASE_NAME <- 'DaymetByNENYHRDCatchments'


setwd('C:/KPONEIL/GitHub/projects/daymetClimateData/gisFiles')

# If the database does not exist then create one
if ( !file.exists(DATABASE_NAME) ) {src_sqlite(DATABASE_NAME, create = T)}

# Connect to the database
database <- dbConnect(SQLite(), DATABASE_NAME)
databaseConnection <- database


#if( year == years[1] & variable == variables[1])

# Make this a function....

# -------------------------------------------------------------------------------------------
# Determine the spatial relationships between catchments, daymet points, and netCDF positions
# -------------------------------------------------------------------------------------------

daymetDirectory
catchmentShapefile
zoneField
variables
years
workingDirectory
zoneField


require(dplyr)
require(lubridate)
require(ncdf4)
require(reshape2)
require(RSQLite)
require(RSQLite.extfuns)


# Create and connect to the database
# ----------------------------------
databasePath <- file.path(DAYMET_DIRECTORY, DATABASE_NAME)

# If the database does not exist then create one
if ( !file.exists( databasePath ) ) 
{ src_sqlite(databasePath, create = T) }

# Connect to the database
database <- dbConnect(SQLite(), databasePath)




catsId <- CATCHMENTS[1:10,]
catsId@data$ZONE_FIELD <- catsId@data[,ZONE_FIELD]




allYears <- seq(from = START_YEAR, to = END_YEAR, by = 1)


spatialIndeces <- determineSpatialRelationships(catchmentsShapefile = catsId, 
                                                  zoneField = ZONE_FIELD,
                                                  daymetDirectory = DAYMET_DIRECTORY,
                                                  variables = VARIABLES,
                                                  years = allYears)

save(spatialIndeces, file = file.path('C:/KPONEIL/workspace', paste(DATABASE_NAME, '_spatialIndeces.RData') ) )






populateDatabase <- function(spatialIndeces, zoneField, daymetDirectory, variables, years, databaseConnection)

mosaicIndeces <- spatialIndeces$mosaicIndeces
tileIndeces <- spatialIndeces$tileIndeces


# Loop through tiles of catchments
for( t in 1:length(tileIndeces) ){
   
  curTile <- tileIndeces[[t]]
  
  cats <- unique(curTile[,zoneField])
  
  # Status update
  print(paste("Beginning variable indexing for tile", t, "of", length(tileIndeces), ", which has", length(cats), "catchments."))
  
  for ( year in years ){
    
    S <- proc.time()[3]
    
    annualRecord <- NULL
    
    for( variable in variables ){
      
      #print(variable)
      
      NetCDF <- nc_open( file.path(daymetDirectory, paste0(variable, '_', year, '.nc4') ) )
      
      tileVar = ncvar_get( nc    = NetCDF, 
                           varid = variable,  
                           start = c(mosaicIndeces$minRow[t], mosaicIndeces$minCol[t], 1), 
                           count = c(mosaicIndeces$countx[t], mosaicIndeces$county[t], 365) )
      
      # Index and replace missval
      # -------------------------
      for( h in 1:length(NetCDF$var) ){
        if ( NetCDF$var[[h]]$name == variable ) {varIndex <- h}
      }
      missingValue <- NetCDF$var[[varIndex]]$missval
      
      # Replace
      tileVar <- replace(tileVar, tileVar == missingValue, NA)
      # -------------------------
      
      # Close the connection
      nc_close(NetCDF)
      
      # Setup storage dataframe
      varPoints <- data.frame(matrix(nrow = nrow(curTile), ncol = 366))
      names(varPoints) <- c(zoneField, 1:365)
      varPoints[,zoneField] <- curTile[,zoneField]
      
      # Loop through days, pulling records from variable array
      for (sel in 1:365){
        select <- as.matrix(data.frame(curTile[c('subRow', 'subCol')], sel) )
        varPoints[,sel+1] <- tileVar[select]
      }
      
      
      
      # Average records by catchment
      # Add duplicate column Because dplyr struggles with variable column names
      varPoints$ZONE <- varPoints[,zoneField]
      varMeans <-  group_by(varPoints, ZONE) %>% 
                    summarise_each(funs(mean))%>%
                    group_by() %>%
                    select( -ZONE )
                  
      # Melt means into output format
      varMelt <- melt(varMeans, id=c(zoneField))
      names(varMelt) <- c(zoneField, 'DayOfYear', variable)
      varMelt$DayOfYear <- as.numeric(varMelt$DayOfYear)
      varMelt$Year <- year
    
      # Store the variables records for the current year
      if( is.null(annualRecord)){
        annualRecord <- varMelt
      } else( annualRecord <- left_join(annualRecord, varMelt, by = c(zoneField, 'Year', 'DayOfYear')))

    }# End variables loop
    
    # Create a "Date" column
    annualRecord$Date <- paste(parse_date_time(paste0(annualRecord$Year,"-", annualRecord$DayOfYear), "y-j", tz = "EST") )
    
    # Return the desired columns 
    annualRecord <- annualRecord[,c(zoneField, 'Date', variables)]    
    
    # Upload to database
    # ------------------
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
    
    # Status update with time elapsed
    print(paste0("Done writing results into database for ", year, " in tile #", t, " of ", length(tileIndeces), ". Elapsed time: ", (proc.time()[3] - S)/60, " minutes."))
    
  }# End years loop 
}# End tile loop 









#..........................................................
#varMeans <- data.frame(matrix(nrow = length(cats), ncol = 366))
#names(varMeans) <- c('FEATUREID', 1:365)
#
#S <- proc.time()[3]
#for( k in 1:length(cats) ){
#  
#  print(k)
#  curCatch <- curTile[curTile$FEATUREID == cats[k], ]
#  
#  for( f in 1:nrow(curCatch) ){
#      
#    newRow <- min(which(is.na(varMeans$FEATUREID)))
#    
#    varMeans$FEATUREID[newRow] <- cats[k]
#    
#    varMeans[newRow,2:366] <- tileVar[curCatch$subRow[f], curCatch$subCol[f], 1:365]
#
#  }
#}

#E <- proc.time()[3]
