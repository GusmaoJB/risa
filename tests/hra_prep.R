# Prepare a list of rasters and or shp to hra.
# takes as input SpatRaster or sf data.frame objects.
# If area is not porvided, it will be derived as a convex hull from the stressor blob
# All rasters and shp will be reprojected after the CRS of the area of interest
# All rasters and shp will be cropped after the area of interest
# The output is a list ready for hra() function
