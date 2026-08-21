library(risa)
library(sf)
library(terra)
library(dplyr)
library(tidyr)
library(readxl)
library(ggplot2)
library(patchwork)

setwd("/home/jojo/Documents/pontal_projects/risa2/tests/")

spp_df <- read.csv("data/species_input.csv")
str_df <- read.csv("data/stressor_input.csv")
criteria <- read.csv("data/criteria.csv")

#Input
sea_lion_df <- readxl::read_excel("/home/jojo/Documents/pontal_projects/risa2/tests/data2/sea_lion_df.xlsx")
seal_rast <- terra::rast("/home/jojo/Documents/pontal_projects/risa2/tests/data2/seal_hotspots.tif")
gillnet_sf <- sf::read_sf("/home/jojo/Documents/pontal_projects/risa2/tests/data2/gillnet_polygons/gillnet_polygons.shp")
criteria2 <- read.csv("/home/jojo/Documents/pontal_projects/risa2/tests/data2/criteria2.csv")

# Treatment
# First we need the path to the folder where we keep the kml files
folder_path <- paste(getwd(), "data2/driftnet_polygons", sep = "/")

# Then we recover the file names in that folder
kml_files <- list.files(
  folder_path,
  pattern = "\\.kml$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Then we load all files as sf vectors and store them in a list
kml_list <- lapply(kml_files, function(f) {
  x <- sf::st_read(f, quiet = TRUE)
  x$source_file <- basename(f)
  return(x)
})

# Then we bind all sf vectors together into a single object
driftnet_sf <- dplyr::bind_rows(kml_list)

# First we need to make a reclassificaiton matrix
re_mat <- reclass_matrix(seal_rast,
                         n_classes = 3,
                         exclude_lowest = FALSE)

# Then use the reclassification matrix to reclassify our raster
seal_rast_rec <- terra::classify(
  seal_rast,
  re_mat,
  include.lowest = TRUE
)

# Plot it
terra::plot(seal_rast_rec, col=hcl.colors(3, rev=TRUE))

# KDE
sea_lion_kde <- get_class_kernel(sea_lion_df)

# Plot
terra::plot(sea_lion_kde, col=hcl.colors(3, rev=TRUE))

# Convert sf to SpatVector
v_gillnet <- terra::vect(gillnet_sf)
v_driftnet <- terra::vect(driftnet_sf)

# Add value 1 of 'Rating' to each polygon
v_gillnet$Rating <- 1
v_driftnet$Rating <- 1

# Reproject to a metric CRS before defining resolution in meters
# Example: WGS 84 / UTM zone 30N (EPSG:32630)
v_gillnet_m <- terra::project(v_gillnet, "EPSG:32630")
v_driftnet_m <- terra::project(v_driftnet, "EPSG:32630")

# Create a template raster (adjust res according to scale!)
template_gillnet <- terra::rast(
  terra::ext(v_gillnet_m),
  res = 500,          # pixel size in meters
  crs = terra::crs(v_gillnet_m))

template_driftnet <- terra::rast(
  terra::ext(v_driftnet_m),
  res = 500,
  crs = terra::crs(v_driftnet_m))

# Rasterize and sum overlapping polygons
gillnet_overlap_r <- terra::rasterize(
  v_gillnet_m,
  template_gillnet,
  field = "Rating",
  fun = "sum",
  background = NA)

driftnet_overlap_r <- terra::rasterize(
  v_driftnet_m,
  template_driftnet,
  field = "Rating",
  fun = "sum",
  background = NA)

# Plot result
par(mfrow = c(1,2))
terra::plot(gillnet_overlap_r,
            main = "Gillnet overlaps",
            col=hcl.colors(5, rev=TRUE))
terra::plot(driftnet_overlap_r,
            main = "Driftnet overlaps",
            col=hcl.colors(5, rev=TRUE))

# Build reclassification matrices for each stressor
gillnet_re_mat <- reclass_matrix(gillnet_overlap_r,
                                 exclude_lowest = FALSE)
driftnet_re_mat <- reclass_matrix(driftnet_overlap_r,
                                  exclude_lowest = FALSE)

# Reclassify maps
gillnet_hotspots <- terra::classify(gillnet_overlap_r,
                                    gillnet_re_mat,
                                    include.lowest = TRUE)
driftnet_hotspots <- terra::classify(driftnet_overlap_r,
                                     driftnet_re_mat,
                                     include.lowest = TRUE)

# Plot reclassified maps
par(mfrow = c(1,2))
terra::plot(gillnet_hotspots,
            main = "Gillnet hotspots",
            col=hcl.colors(3, rev=TRUE))
terra::plot(driftnet_hotspots,
            main = "Driftnet hotspots",
            col=hcl.colors(3, rev=TRUE))

# Preparing input lists
species_list <- list(Sea.lion = sea_lion_kde,
                     Seal = seal_rast_rec)

stressor_list <- list(Gillnet = gillnet_hotspots,
                      Driftnet = driftnet_hotspots)

# Harmonize maps for ByRA
byra_input <- byra_prep(species_list,
                        stressor_list,
                        quiet = FALSE)

byra3 <- quick_byra(byra_input,
                    criteria = criteria2,
                    equation = "euclidean",
                    quiet = FALSE,
                    return_crs = "4326")

hra(byra_input, )

byra_input

