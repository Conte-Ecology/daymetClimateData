
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

source('C:/KPONEIL/GitHub/personal/dataProcessing/scripts/functions/tileCatchmentsShapefile.R')


proj4.WGS  <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

proj4.Lambert <- "+proj=lcc +ellps=WGS84 +datum=WGS84 +lat_1=25 +lat_2=60 +lat_0=42.5 +lon_0=-100 +x_0=0 +y_0=0"




# Temporal range
START_YEAR <- 1980
END_YEAR <- 2013

years <- seq(from = START_YEAR, to = END_YEAR, by = 1)

# Variables
variables <- c('tmax', 'tmin', 'prcp', 'dayl', 'srad', 'vp', 'swe')

daymetDirectory <- 'F:/KPONEIL/SourceData/climate/DAYMET/raw'

tableName <- "climateRecord"

catchmentsShapefile <- readShapePoly('C:/KPONEIL/GitHub/projects/daymetClimateData/gisFiles/shapefiles/NENYHRD_AllCatchments_DaymetPrj.shp', proj4string=CRS(proj4.Lambert))




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
transformShapefile <- spTransform(catchmentsShapefile, CRS(proj4.WGS), class = "SpatialPolygonsDataFrame")

catchTiles <- tileCatchmentsShapefile(transformShapefile)

spatialNCDF <- nc_open(file.path(daymetDirectory, paste0( variables[1], '_', years[1], '.nc4') ) )  #netcdf     

lat = ncvar_get ( nc=spatialNCDF, varid="lat", start = c(1,1), count = c(spatialNCDF$var$lat$varsize[1], spatialNCDF$var$lat$varsize[2]) )
lon = ncvar_get ( nc=spatialNCDF, varid="lon", start = c(1,1), count = c(spatialNCDF$var$lon$varsize[1], spatialNCDF$var$lon$varsize[2]) )

tileIndeces <- list() # Change to ncdfIndeces
tilePoints <- list() # Change to tileIndeces

for(t in 1:length(catchTiles) ){
  
  tile <- catchTiles[[t]]
  
  print(paste("Beginning spatial analysis of tile", t, "of", length(catchTiles), ", which has", nrow(tile), "catchments."))
  
  ## Get the extent of the shapefile
  EXT <- extent(tile)
  
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
  
  tileIndeces[[t]] <- data.frame(tile = t, minRow = minRow, maxRow = maxRow, minCol = minCol, maxCol = maxCol, countx = countx, county = county)

  trimLat = ncvar_get( nc = spatialNCDF, varid="lat",     start = c(minRow, minCol),    count = c(countx,county) )
  trimLon = ncvar_get( nc = spatialNCDF, varid="lon",     start = c(minRow, minCol),    count = c(countx,county) )
  
  masterCoords <- as.data.frame( cbind( as.vector(trimLon), as.vector(trimLat) ) )
  colnames(masterCoords) <- c("Longitude", "Latitude")
  masterCoordsMatrix <- as.matrix(masterCoords)
  masterCoordSpPts   <- SpatialPoints(masterCoords, proj4string = CRS(proj4string(transformShapefile)))
  
  # Overlay points on catchment shapefile
  overPoints <- over(masterCoordSpPts, tile)
  
  spatJoin <- cbind(overPoints, masterCoords)
  
  join <- spatJoin[which(!is.na(spatJoin$FEATUREID)),] #or... select(spatJoin, FEATUREID != NA)
  
  missingCats <- tile[!tile@data$FEATUREID %in% join$FEATUREID,]
  
  centroids <- data.frame(missingCats@data$FEATUREID, coordinates(missingCats) )
  names(centroids) <- c('FEATUREID', 'LON', 'LAT')
  
  fids <- missingCats@data$FEATUREID
  
  miss <- as.data.frame(matrix(nrow = length(fids), ncol = 3))
  names(miss) <- c('FEATUREID', 'Longitude', 'Latitude')

  for ( i in seq_along(fids) ) {
    
    tempLat <- centroids$LAT[centroids$FEATUREID == fids[i]]
    tempLon <- centroids$LON[centroids$FEATUREID == fids[i]]
      
    distances <- spDistsN1(masterCoordsMatrix, c(tempLon, tempLat), longlat = TRUE)
    minDist <- min(distances)
    distpos <- which(distances == minDist)[1]
        
    miss$FEATUREID[i] <- fids[i] 
    miss$Longitude[i] <- masterCoords[distpos, 1]
    miss$Latitude[i]  <- masterCoords[distpos, 2]
  }

  catchPoints <- rbind(join[,c('FEATUREID', 'Longitude', 'Latitude')], miss)
  
  catchPoints$subRow <- NA; catchPoints$subCol <- NA
  
  for ( m in 1:nrow(catchPoints) ){
      
    # Find th position in the array of the variable
    position <- which(trimLon == catchPoints$Longitude[m] & trimLat == catchPoints$Latitude[m], arr.in = TRUE)
  
    catchPoints$subRow[m] <- as.numeric(position[,1])
    catchPoints$subCol[m] <- as.numeric(position[,2])
  }

  tilePoints[[t]] <- catchPoints
  names(tilePoints)[t] <- paste0("Tile", t)
  
  print(paste("Tile", t, "complete."))
}


tileIndeces <- rbind_all(tileIndeces)

nc_close(spatialNCDF)


# End spatial function....
# ---------------------------------------------------------------------------------------------------------------
rm(catchmentsShapefile, transformShapefile, lat, lon, trimLat, trimLon)# This changes with the previous step becoming a function

save(tileIndeces, tilePoints, file = file.path('C:/KPONEIL/workspace', paste(DATABASE_NAME, '_tileInfo.RData') ) )



# Loop through tiles of catchments
mainbeg <- proc.time()[3]
for( t in 1:length(tilePoints) ){
   
  curTile <- tilePoints[[t]]
  
  cats <- unique(curTile$FEATUREID)
  
  # Status update
  print(paste("Beginning variable indexing for tile", t, "of", length(tilePoints), ", which has", length(cats), "catchments."))
  
  for ( year in years ){
    
    S <- proc.time()[3]
    
    annualRecord <- NULL
    
    for( variable in variables ){
      
      print(variable)
      
      #NetCDF <- nc_open( file.path( daymetDirectory, paste0(variable, '_', year, '.nc4') ) )    #netcdf     
      NetCDF <- nc_open( file.path('F:/KPONEIL/SourceData/climate/DAYMET/raw', paste0(variable, '_', year, '.nc4') ) )
      
      tileVar = ncvar_get( nc    = NetCDF, 
                           varid = variable,  
                           start = c(tileIndeces$minRow[t], tileIndeces$minCol[t], 1), 
                           count = c(tileIndeces$countx[t], tileIndeces$county[t], 365) )
      
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
      names(varPoints) <- c('FEATUREID', 1:365)
      varPoints$FEATUREID <- curTile$FEATUREID
      
      # Loop through days, pulling records from variable array
      for (sel in 1:365){
        select <- as.matrix(data.frame(curTile[c('subRow', 'subCol')], sel) )
        varPoints[,sel+1] <- tileVar[select]
      }
      
      # Average records by catchment
      varMeans <- varPoints %>% group_by(FEATUREID) %>% summarise_each(funs(mean))
      
      # Melt means into output format
      varMelt <- melt(varMeans, id=c("FEATUREID"))
      names(varMelt) <- c('FEATUREID', 'DayOfYear', variable)
      varMelt$DayOfYear <- as.numeric(varMelt$DayOfYear)
      varMelt$Year <- year
    
      # Store the variables records for the current year
      if( is.null(annualRecord)){
        annualRecord <- varMelt
      } else( annualRecord <- left_join(annualRecord, varMelt, by = c('FEATUREID', 'Year', 'DayOfYear')))

    }# End variables loop
    
    # Create a "Date" column
    annualRecord$Date <- paste(parse_date_time(paste0(annualRecord$Year,"-", annualRecord$DayOfYear), "y-j", tz = "EST") )
    
    # Return the desired columns 
    annualRecord <- annualRecord[,c('FEATUREID', 'Date', variables)]    
    
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
    print(paste0("Done writing results into database for ", year, " in tile #", t, " of ", length(tilePoints), ". Elapsed time: ", (proc.time()[3] - S)/60, " minutes."))
    
  }# End years loop 
}# End tile loop 
mainend <- proc.time()[3]

print(paste0("Total processing time: ", (mainend - mainbeg)/3600))
# 31.53 hours for NENY HRD catchments








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
