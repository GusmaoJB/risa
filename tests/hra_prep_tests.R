# Setting my working directory
setwd("/home/jojo/Documents/pontal_projects/risa2/tests")

# Loading packages
library(risa)
library(terra)

# Loading data
species_df <- read.csv("species_distr.csv")
stressor_df <- read.csv("stressor_distr.csv")
criteria <- read.csv("criteria.csv")

# Running ByRA
byra <- quick_byra(species_df, stressor_df, criteria)

r1 <- byra$kde_maps$species_kernel_maps$Dolphin$raster
r2 <- convert_to_decimal_degrees(byra$kde_maps$species_kernel_maps$Turtle$shp)

r3 <- convert_to_decimal_degrees(byra$kde_maps$stressor_kernel_maps$Gillnet$raster)
r4 <- byra$kde_maps$stressor_kernel_maps$Trawling$shp

r1
r2
r3
r4

spp_list <- list(r1, r2)
str_list <- list(r3, r4)

names(spp_list)

test1 <- byra_prep3(spp_list, str_list, reclass_cat = TRUE)

terra::plot(test1$species_distributions$species_1)
terra::plot(test1$species_distributions$species_2)
terra::plot(test1$species_kernel_maps$species_1)
terra::plot(test1$species_kernel_maps$species_2)
terra::plot(test1$stressor_kernel_maps$stressor_1)
terra::plot(test1$stressor_kernel_maps$stressor_2)
par(mfrow=c(2,2))

terra::plot(test1$overlap_maps$species_1$stressor_1)
terra::plot(test1$overlap_maps$species_1$stressor_2)
terra::plot(test1$overlap_maps$species_2$stressor_1)
terra::plot(test1$overlap_maps$species_2$stressor_2)
plot(test1$area_of_interest, col="red")


test1$overlap_maps$species_2$stressor_2
test1$area_of_interest
terra::plot(test1$overlap_maps$species_2$stressor_2)
lines(test1$area_of_interest, col="red", lwd=2)

test1$area_of_interest

plot(sf::st_geometry(test1$area_of_interest), border = "red")  # set extent from AOI
terra::plot(test1$overlap_maps$species_1$stressor_1, add=TRUE)
terra::plot(test1$overlap_maps$species_1$stressor_2, add=TRUE)
terra::plot(test1$overlap_maps$species_2$stressor_1, add=TRUE)
terra::plot(test1$overlap_maps$species_2$stressor_2, add=TRUE)
terra::plot(test1$species_kernel_maps$species_1, add=TRUE)
terra::plot(test1$species_kernel_maps$species_2, add=TRUE)
terra::plot(test1$stressor_kernel_maps$stressor_1, add=TRUE)
terra::plot(test1$stressor_kernel_maps$stressor_2, add=TRUE)

