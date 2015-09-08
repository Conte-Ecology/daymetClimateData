library(rgdal)
library(ncdf4)

source("C:/KPONEIL/GitHub/projects/daymetClimateData/scripts/functions/createDaymetPointShapefile.R")


bound <- readOGR(dsn = "C:/KPONEIL/workspace", layer = "nenyNHRD25kmbuf_daymet")

boundSP <- spTransform(bound, CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0") )

NCDF <- nc_open('F:/KPONEIL/SourceData/climate/DAYMET/raw/tmax_1985.nc4')    #netcdf     

test <- createDaymetPointShapefile(boundaryShapefile = boundSP,
                                       NetCDF = NCDF, 
                                       includeVariable = T, 
                                       dOY = 307)



setwd("C:/KPONEIL/workspace")
writeOGR(test,  
         ".", 
         layer = "dayPtNov03_1985tmax", 
         driver = "ESRI Shapefile")
