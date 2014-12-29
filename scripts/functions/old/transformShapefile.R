# shapefile - the shapefile to transform
# newProjection - the proj4string of the coordinate system/projection that the shapefile is getting projected into

transformShapefile <- function(shapefile, newProjection){
  # It is important to note that this simple function requires in depth metadata on the transformation between the two projections.
  # For example, the "+towgs84=0,0,0" portion of the proj4string of a shapefile indicates the transformation parameters. (See: http://proj.maptools.org/gen_parms.html)
  # If this information is not known, it is best to transform the shapefile in a GIS program ahead of time.
  
  require(sp)
  
  # Transform the shapefile to be in the Daymet projection if necessary
  if(!proj4string(shapefile) == newProjection){
    
    print("Spatially transforming shapefile...")
    shapefilePROJ <- spTransform(shapefile, newProjection, class = "SpatialPolygonsDataFrame")
  } else(shapefilePROJ <- shapefile )
  
  return(shapefilePROJ)
}