
library(ncdf4)

# Temporal range
START_YEAR <- 1980
END_YEAR <- 2013

# Variables
variables <- c('tmax', 'tmin', 'prcp', 'dayl', 'srad', 'vp', 'swe')

daymetDirectory <- 'F:/KPONEIL/SourceData/climate/DAYMET/raw'


years <- seq(from = 1980, to = 2013, by = 1)      

for (year in years){
  for (var in variables){
    
    NetCDF <- nc_open(file.path(daymetDirectory, paste0(var, '_',year , '.nc4')))    #netcdf     
    
    lat = ncvar_get ( nc=NetCDF, varid="lat", start = c(1000,1400), count = c(1, 1) )
    lon = ncvar_get ( nc=NetCDF, varid="lon", start = c(1000,1400), count = c(1, 1) )
    
    print(paste(lat, lon))
  }
}