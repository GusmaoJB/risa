library(risa)
library(terra)
library(sf)

# For even more complex cases, we have to use hra()

# Area of interest
aio <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/area_of_interest/area_of_interest.shp")

# Species
pont <- terra::rast("/home/jojo/Documents/pontal_projects/risa_example_maps/species_kernel_maps/species1/raster/raster.tif")
sota <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/species_kernel_maps/species2/shp/shp.shp")

# Stressors
gill <- terra::rast("/home/jojo/Documents/pontal_projects/risa_example_maps/stressor_kernel_maps/stressor1/raster/raster.tif")
traw <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/stressor_kernel_maps/stressor1/shp/shp.shp")

# criteria
crit2 <- read.csv("/home/jojo/Documents/pontal_projects/risa_example_shp/criteria.csv")

# quick plot
terra::plot(pont)
terra::plot(sota)

# Prepare for HRA
pont <- convert_to_decimal_degrees(pont)
pont
sota
gill
traw

standard_maps <- byra_prep3(list(pont, sota), list(gill, traw), quiet=FALSE)
terra::plot(standard_maps$species_kernel_maps$species_1)
terra::plot(standard_maps$stressor_kernel_maps$stressor_1)
terra::plot(standard_maps$area_of_interest)
terra::plot(standard_maps$overlap_maps$species_1$stressor_1)

################################################################################
# When we already have density maps for species/habitats
