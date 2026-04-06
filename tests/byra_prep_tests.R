library(risa)
library(terra)
library(sf)

# For even more complex cases, we have to use hra()

# Area of interest
aio <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/area_of_interest/area_of_interest.shp")

# Species
pont <- terra::rast("/home/jojo/pontoporia.tif")
sota <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/species_kernel_maps/species2/shp/shp.shp")

# Stressors
gill <- terra::rast("/home/jojo/Documents/pontal_projects/risa_example_maps/stressor_kernel_maps/stressor1/raster/raster.tif")
traw <- read_sf("/home/jojo/Documents/pontal_projects/risa_example_maps/stressor_kernel_maps/stressor2/shp/shp.shp")

# criteria
crit2 <- read.csv("/home/jojo/Documents/pontal_projects/risa_example_shp/criteria.csv")

# quick plot
terra::plot(pont)
terra::plot(sota)
terra::plot(gill)
terra::plot(traw)

# Prepare for HRA
pont
sota
gill
traw

standard_maps <- byra_prep(list(pont, sota), list(gill, traw), quiet=FALSE)

terra::plot(standard_maps$species_kernel_maps$species_1)
terra::plot(standard_maps$species_kernel_maps$species_2)
terra::plot(standard_maps$stressor_kernel_maps$stressor_1)
terra::plot(standard_maps$stressor_kernel_maps$stressor_2)
terra::plot(standard_maps$overlap_maps$species_1$stressor_1)
terra::plot(standard_maps$overlap_maps$species_1$stressor_2)
terra::plot(standard_maps$overlap_maps$species_2$stressor_1)
terra::plot(standard_maps$overlap_maps$species_2$stressor_2)


# The output is not keeping the names of the original objects!!!
################################################################################
