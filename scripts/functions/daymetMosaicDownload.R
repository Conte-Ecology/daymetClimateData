
# years - numeric, a vector of years to download
# variables - vector of character strings, the variables to download as named by Daymet (e.g. "prcp")
# destinationFolder - A character string of the file path to the folder to download the files to
# retryFailedDownloads - TRUE/FALSE A check of whether the file was downloaded successfully. If TRUE, the check will be performed and the file redownloaded if the original is corrpupt.

daymetMosaicDownload <- function(years, variables, destinationFolder, retryFailedDownloads){
  
  # Default to retrying corrupted files
  if( is.null(retryFailedDownloads) ) { retryFailedDownloads = TRUE }

  require(ncdf4)
  
  # Loop through years
  for ( year in years ){
    
    # Loop through variables
    for( var in variables){
            
      #THREDDS server address
      address <- paste0('http://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1219/', year, '/', var, '_', year, '.nc4')
      
      # Name the output file
      outFile <- file.path(destinationFolder, paste0(var, '_', year, '.nc4'))
            
      # Alert the user if the file already exists in the target directory. Also checks for file corruption.
      if(file.exists(outFile)){
        print(paste0("File '", outFile ,"' already exists in download directory."));

        # Alerts the user if the existing file is corrupt
        if(is(try(nc_open(outFile), silent=T),"try-error")){
          print(paste0("Existing file: '", outFile, "' is corrupt. Suggest deleting this file."))
        }# end corruption check
      }# end file.exists check
      
      # If the file doesn't exist, download it
      while(!file.exists(outFile)){
        
        # Start time
        beg <- proc.time()[3]
        
        # Download the file
        download.file(url      = address, 
                      destfile = outFile, 
                      quiet    = FALSE, 
                      mode     = 'wb')
            
        # Print download time
        print(paste0("Download took ", (proc.time()[3] - beg)/3600, " hours.") )
        
        # Test to see if the file downloaded correctly. If not, retry download.
        fileTest <- try(nc_open(outFile), silent=T)
        if(is(fileTest,"try-error") & retryFailedDownloads) {
          file.remove(outFile);print("File corrupt. Removing and redownloading....")
        } else{nc_close(fileTest); rm(fileTest)}
      }# end while 
    }# end variable loop
  }# end year loop
}# end function
