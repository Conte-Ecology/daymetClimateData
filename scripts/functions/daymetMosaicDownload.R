# startYear - numeric, the first year in sequence to begin downloading files for
# endYear - numeric, the last year in sequence to begin downloading files for
# variables - vector of character strings, the variables to download as named by Daymet (e.g. "prcp")
# destinationFolder - A character string of the file path to the folder to download the files to
# retryFailedDownloads - TRUE/FALSE A check of whether the file was downloaded successfully. If TRUE, the check will be performed and the file redownloaded if the original is corrpupt.

daymetMosaicDownload <- function(startYear, endYear, variables, destinationFolder, retryFailedDownloads){
  
  # Need to include a section that checks if the downloaded file is not corrupt and can be openned
  require(ncdf4)
  
  # Loop through years
  for ( year in seq(from = startYear, to = endYear, by = 1) ){
    
    # Loop through variables
    for( var in variables){
            
      #THREDDS server address
      address <- paste0('http://thredds.daac.ornl.gov/thredds/fileServer/ornldaac/1219/', year, '/', var, '_', year, '.nc4')
      
      # Name the output file
      outFile <- file.path(destinationFolder, paste0(var, '_', year, '.nc4'))
            
      # Don't download if the file already exists
      if(file.exists(outFile)){
        print(paste0("File '", outFile ,"' already exists in download directory."));
        #existingFileTest <- try(nc_open(outFile), silent=T);
        if(is(try(nc_open(outFile), silent=T),"try-error")) {print(paste0("Existing file: '", outFile, "' is corrupt. Suggest deleting this file."))}
      }# end if
      
      # If the file doesn't exist, download it
      while(!file.exists(outFile)){
        
        # Start time
        beg <- proc.time()[3]
        
        # Download the file
        download.file(url = address, destfile = outFile, quiet = FALSE, mode = 'wb')
          
        # Time download
        runTime <- (proc.time()[3] - beg)/3600
          
        # Print download time
        print(paste0("Download took ", runTime, " hours.") )
        
        # Test to see if the file downloaded correctly. If not, re-attempt.
        fileTest <- try(nc_open(outFile), silent=T)
        if(is(fileTest,"try-error") & retryFailedDownloads) {file.remove(outFile);print("File corrupt. Removing and redownloading....")} else {nc_close(fileTest); rm(fileTest)}
      }# end while 
      
    }# end variable loop
  }# end year loop
}# end function
