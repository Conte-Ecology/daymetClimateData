# databaseFilePath - the file path of the SQLite database
# tableName - character string, the name of the table in the database
# startDate - a character string of the start date of the record to pull (Format: "yyyy-mm-dd")
# endDate  - a character string of the end date of the record to pull (Format: "yyyy-mm-dd")
# featureIDs - a vector of numeric values, FEATUREIDs to pull from the database
# variables - vector of character strings, the variables to download as named by Daymet (e.g. "prcp")

returnRecordFromDaymetDB <- function(databaseFilePath, tableName, startDate, endDate, featureIDs, variables){
  
  require(dplyr)
  require(lubridate)
  
  # Connect to the database
  dBase <- src_sqlite(databaseFilePath, create = F)
  
  # Set the range to pull from the database
  dateRange <- paste(seq(from = ymd(startDate, tz = "EST"),
                           to = ymd(endDate, tz = "EST"),
                           by = "day"))
          
  # Create iteration storage
  tempRecord <- list()
  
  # Loop through variable tables in the database
  for ( i in seq_along(featureIDs) ){
    
    # Select the current variable table
    sourceTable <- tbl(dBase, sql(paste0("SELECT * FROM ", tableName)))
      
    # Pull the FeatureIDs and dates
    tempDF <- collect(filter(sourceTable, 
                                      FEATUREID == featureIDs[i], 
                                      Date %in% dateRange))
    tempRecord[[i]] <- tempDF[,c('FEATUREID', 'Date', variables)]
  }
    
  # Join all catchment records
  outRecord <- rbind_all(tempRecord)

  # Return the record
  return(outRecord)
}