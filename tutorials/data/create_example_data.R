# Creating example data for exaple 2
library(sf)
library(dplyr)
library(terra)

# Seal
set.seed(555)
seal_df <- data.frame(long = rnorm(100, 0, 0.15),
                      lat = rnorm(100, 0, 0.2), species = "Seal")
seal_kde <- get_class_kernel(seal_df, continuous = TRUE, output_min = 0, n_classes = 1)
# terra::writeRaster(seal_kde, "/home/jojo/Documents/pontal_projects/risa2/tests/data/seal_hotspots.tif", overwrite = TRUE)

# Sea lion
set.seed(444)
groups <- sample(c(1:21), size = 90, replace = TRUE, prob = seq(from = 100, to = 1, length.out = 21))
sea_lion_df <- data.frame(long = rnorm(90, 0, 0.2),
                          lat = rnorm(90, 0, 0.16), 
                          species = "Sea lion",
                          group_size = groups)
# write.csv(sea_lion_df, "/home/jojo/Documents/pontal_projects/risa2/tests/data/sea_lion_df.csv", row.names = FALSE)

# Gillnet
set.seed(666)
str_df2 <- rbind(data.frame(long = rnorm(100, 0, 0.18),
                            lat = rnorm(100, 0, 0.2), fisherman = "fisherman1"),
                 data.frame(long = rnorm(50, 0.2, 0.2),
                            lat = rnorm(50, 0, 0.17), fisherman = "fisherman2"),
                 data.frame(long = rnorm(70, 0, 0.1),
                            lat = rnorm(70, -0.1, 0.15), fisherman = "fisherman3"),
                 data.frame(long = rnorm(30, -0.2, 0.1),
                            lat = rnorm(30, -0.1, 0.15), fisherman = "fisherman4"),
                 data.frame(long = rnorm(40, 0.1, 0.1),
                            lat = rnorm(40, 0.2, 0.15), fisherman = "fisherman5"))

str_sf <- st_as_sf(
  str_df2,
  coords = c("long", "lat"),
  crs = 4326
)

hulls_sf <- str_sf %>%
  group_by(fisherman) %>%
  summarise(
    geometry = st_convex_hull(st_union(geometry)),
    .groups = "drop"
  )

#st_write(hulls_sf, "/home/jojo/Documents/pontal_projects/risa2/tests/data2/gillnet_polygons.shp", delete_dsn = TRUE)

# Driftnet
set.seed(231)
str_df3 <- rbind(data.frame(long = rnorm(100, 0, 0.18),
                            lat = rnorm(100, 0, 0.2), fisherman = "fisherman1"),
                 data.frame(long = rnorm(50, 0.2, 0.2),
                            lat = rnorm(50, 0, 0.17), fisherman = "fisherman2"),
                 data.frame(long = rnorm(30, -0.2, 0.1),
                            lat = rnorm(30, -0.1, 0.15), fisherman = "fisherman3"),
                 data.frame(long = rnorm(40, 0.1, 0.1),
                            lat = rnorm(40, 0.2, 0.15), fisherman = "fisherman4"))

str_sf2 <- st_as_sf(
  str_df3,
  coords = c("long", "lat"),
  crs = 4326
)

hulls_sf2 <- str_sf2 %>%
  group_by(fisherman) %>%
  summarise(
    geometry = st_convex_hull(st_union(geometry)),
    .groups = "drop"
  )

hulls_sf2 <- st_transform(hulls_sf2, 4326)

hulls_list <- split(hulls_sf2, hulls_sf2$fisherman)
plot(hulls_list$fisherman4)

setwd("/home/jojo/Documents/pontal_projects/risa2/tests/data2")

# lapply(names(hulls_list), function(nm) {
#   st_write(hulls_list[[nm]], paste0(nm, ".kml"), 
#            driver = "KML", delete_dsn = TRUE)})

# Load Sea lion data from an Excel file
sea_lion_df <- readxl::read_excel("sea_lion_df.xlsx")

# Load Seal data from a .tif raster
seal_rast <- terra::rast("seal_hotspots.tif")

# Load Gillnet data from a shapefile
gillnet_sf <- sf::read_sf("gillnet_polygons/gillnet_polygons.shp")

# Load criteria data
criteria2 <- read.csv("criteria2.csv")


# Load Driftnet data from multiple KML files

# First we need the path to the folder where we keep the kml files
folder_path <- paste(getwd(), "driftnet_polygons", sep = "/")

# Then we recover the file names in that folder
kml_files <- list.files(
  folder_path,
  pattern = "\\.kml$",
  full.names = TRUE,
  ignore.case = TRUE
)

# Then we load all files as sf vectors and store then in a list
kml_list <- lapply(kml_files, function(f) {
  x <- st_read(f, quiet = TRUE)
  x$source_file <- basename(f)
  return(x)
})

# Then we bind all sf vectors together into a single object
driftnet_sf <- bind_rows(kml_list)


# Prepare data
re_mat <- reclass_matrix(seal_rast, n_classes = 3, exclude_lowest = FALSE)

r_valid <- terra::classify(
  seal_rast,
  re_mat,
  include.lowest = TRUE
)

terra::plot(r_valid)


sea_lion_kde <- get_class_kernel(sea_lion_df, group_size = "group_size", radius = 5000)
terra::plot(sea_lion_kde)

# Convert sf to SpatVector
v_gillnet <- terra::vect(gillnet_sf)
v_driftnet <- terra::vect(driftnet_sf)

# Add value 1 of 'Rating' to each polygon
v_gillnet$Rating <- 1
v_driftnet$Rating <- 1

# Reproject to a metric CRS before defining resolution in meters
# Example: SIRGAS 2000 / UTM 22S
v_gillnet_m <- terra::project(v_gillnet, "EPSG:31982")
v_driftnet_m <- terra::project(v_driftnet, "EPSG:31982")

# Create a template raster (adjust res according to the scale of your data!)
v_gillnet_r <- terra::rast(
  terra::ext(v_gillnet_m),
  res = 500,          # pixel size in meters
  crs = terra::crs(v_gillnet_m))

v_driftnet_r <- terra::rast(
  terra::ext(v_driftnet_m),
  res = 500,          
  crs = terra::crs(v_driftnet_m))

# Rasterize and sum overlapping polygons
gillnet_overlap_r <- terra::rasterize(
  v_gillnet_m,
  v_gillnet_r,
  field = "Rating",
  fun = "sum",
  background = NA)

driftnet_overlap_r <- terra::rasterize(
  v_driftnet_m,
  v_driftnet_r,
  field = "Rating",
  fun = "sum",
  background = NA)

# Plot result
par(mfrow = c(1,2))
terra::plot(gillnet_overlap_r, main = "Gillnet overlaps")
terra::plot(driftnet_overlap_r, main = "Driftnet overlaps")
dev.off()

# Build reclassification matrices
gillnet_re_mat <- reclass_matrix(gillnet_overlap_r, exclude_lowest = FALSE)
driftnet_re_mat <- reclass_matrix(driftnet_overlap_r, exclude_lowest = FALSE)

# Reclassify maps
gillnet_hotspots <- terra::classify(gillnet_overlap_r, gillnet_re_mat, include.lowest = TRUE)
driftnet_hotspots <- terra::classify(driftnet_overlap_r, driftnet_re_mat, include.lowest = TRUE)

# Plot reclassified maps
par(mfrow = c(1,2))
terra::plot(gillnet_hotspots, main = "Gillnet hotspots")
terra::plot(driftnet_hotspots, main = "Driftnet hotspots")


# Build reclassification matrices
gillnet_re_mat <- reclass_matrix(gillnet_overlap_r, exclude_lowest = FALSE)
driftnet_re_mat <- reclass_matrix(driftnet_overlap_r, exclude_lowest = FALSE)

# Reclassify maps
gillnet_hotspots <- terra::classify(gillnet_overlap_r, gillnet_re_mat, include.lowest = TRUE)
driftnet_hotspots <- terra::classify(driftnet_overlap_r, driftnet_re_mat, include.lowest = TRUE)

# Plot reclassified maps
par(mfrow = c(1,2))
terra::plot(gillnet_hotspots, main = "Gillnet hotspots")
terra::plot(driftnet_hotspots, main = "Driftnet hotspots")
dev.off()

# Preparing input lists
species_list <- list(`Sea lion` = sea_lion_kde,
                     Seal = seal_rast_rec)

stressor_list <- list(Gillnet = gillnet_hotspots,
                      Driftnet = driftnet_hotspots)

# Harmonize maps for ByRA
library(ggplot2)
byra_input <- byra_prep(species_list, stressor_list, quiet = FALSE)
class(byra_input)
risaplot2(byra_input)






