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

# This time, we ask for decimal degrees by setting return_crs = "4326"
byra4 <- quick_byra(byra_input,
                    criteria = criteria2,
                    equation = "euclidean",
                    return_crs = "4326",
                    quiet = FALSE)

# Loading the packages (install if you don't have them)
library(geodata)
library(ggspatial)

# Load world map data
world_map <- geodata::world(resolution = 3, path = "data/")

# We don't need all the countries, so we extract only what we need
ghana <- world_map[world_map$NAME_0 %in% c("Ghana"), ]

# Convert it to an sf vector, so it is easier to plot with ggplot2
ghana_sf <- st_as_sf(ghana)

# Making ggplot-type risk maps with risaplot
risk_plots <- risaplot(byra4)
byra4$area_of_interest

library(geodata)
ghana_adm <- geodata::gadm(country = "GHA", level = 2, path = tempdir())  # level varies by country
ghana_adm_sf <- st_as_sf(ghana_adm)

risk_plots[[6]] +
  geom_sf(data=ghana_sf, fill="gray", col="black") +
  geom_sf(data=ghana_adm_sf, fill="#00000000", col="gray50") +
  coord_sf(xlim=c(-1.15,1),
           ylim=c(4,6)) +
  annotation_scale(location = "br", height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "br",
                         pad_y=unit(1, "cm"),
                         width = unit(0.6, "cm"),
                         height = unit(0.8, "cm"))


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

# This time, we ask for decimal degrees by setting return_crs = "4326"
byra4 <- quick_byra(byra_input,
                    criteria = criteria2,
                    equation = "euclidean",
                    return_crs = "4326",
                    quiet = FALSE)

# Loading the packages (install if you don't have them)
library(geodata)
library(ggspatial)

# Load world map data
world_map <- geodata::world(resolution = 3, path = "data/")

# We don't need all the countries, so we extract only what we need
ghana <- world_map[world_map$NAME_0 %in% c("Ghana"), ]

# Convert it to an sf vector, so it is easier to plot with ggplot2
ghana_sf <- st_as_sf(ghana)

# Making ggplot-type risk maps with risaplot
risk_plots <- risaplot(byra4)
byra4$area_of_interest

library(geodata)
ghana_adm <- geodata::gadm(country = "GHA", level = 2, path = tempdir())  # level varies by country
ghana_adm_sf <- st_as_sf(ghana_adm)

risk_plots[[6]] +
  geom_sf(data=ghana_sf, fill="gray", col="black") +
  geom_sf(data=ghana_adm_sf, fill="#00000000", col="gray50") +
  coord_sf(xlim=c(-1.15,1),
           ylim=c(4,6)) +
  annotation_scale(location = "br", height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "br",
                         pad_y=unit(1, "cm"),
                         width = unit(0.6, "cm"),
                         height = unit(0.8, "cm"))

# First, let's extract what we need
seal_total <- byra4$Seal$total_raw
sealion_total <- byra4$Sea.lion$total_raw
ecosys_total <- byra4$ecosys_risk_raw
aio_sf <- byra4$area_of_interest

# Since risk maps are SpatRaster objects, so we need to convert them to data.frame
seal_total_df <- as.data.frame(seal_total, xy=TRUE, na.rm=TRUE)
sealion_total_df <- as.data.frame(sealion_total, xy=TRUE, na.rm=TRUE)
ecosys_total_df <- as.data.frame(ecosys_total, xy=TRUE, na.rm=TRUE)

# This time, we will define a different palete for species and ecosystem maps
spp_palette <- c("#DCE95B", "#32AE7C", "#255668")
eco_palette <- c("#FFD99F", "#F7606D", "#7D1D67")

# Since we will standardize the aesthetics of the plots, we can define a template
template <- ggplot() +
  geom_sf(data=ghana_sf, fill="gray", col="black") +
  geom_sf(data=ghana_adm_sf, fill="#00000000", col="gray50") +
  coord_sf(xlim=c(-1.15,1),
           ylim=c(4,6)) +
  annotation_scale(location = "br",
                   height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "br",
                         pad_y=unit(1, "cm"),
                         width = unit(0.6, "cm"),
                         height = unit(0.8, "cm")) +
  theme_bw()

# The layers are SpatRasters
sp1 <- template +
  geom_tile(data=seal_total_df, aes(x=x, y=y, fill=Rating)) +
  scale_color_gradientn(colors=spp_palette) +
  theme(panel.grid = element_blank(),
        legend.position = "none")


# Since we will standardize the aesthetics of the plots, we can define a template
template <- ggplot() +
  geom_sf(data=ghana_sf, fill="gray", col="black") +
  geom_sf(data=ghana_adm_sf, fill="#00000000", col="gray50") +
  coord_sf(xlim=c(-1.15,1),
           ylim=c(4,6)) +
  annotation_scale(location = "br",
                   height = unit(0.2, "cm")) +
  annotation_north_arrow(location = "br",
                         pad_y=unit(1, "cm"),
                         width = unit(0.6, "cm"),
                         height = unit(0.8, "cm")) +
  theme_bw() +
  labs(x = "Longitude", y="Latitude")

sp1 <- template +
  geom_tile(data=seal_total_df, aes(x=x, y=y, fill=Rating)) +
  scale_fill_gradientn(colors=spp_palette) + # set our color palette
  theme(panel.grid = element_blank(),
        legend.position = "none") + # we will omit the color guides in this one
  ggtitle("a) Seal (total risk)") + # Includes a title
  annotation_scale(location = "br",
                   height = unit(0.2, "cm")) +
  theme_bw()  +
  labs(x = "Longitude", y="Latitude")

# Sea lion
sp2 <- template +
  geom_tile(data=sealion_total_df, aes(x=x, y=y, fill=Rating)) +
  scale_fill_gradientn(colors=spp_palette) +
  theme(panel.grid = element_blank()) +
  ggtitle("b) Sea lion (total risk)")

# All pennipeds combined
tot <- template +
  geom_tile(data=ecosys_total_df, aes(x=x, y=y, fill=lyr.1)) +
  scale_fill_gradientn(colors=eco_palette) +
  theme(panel.grid = element_blank()) +
  ggtitle("c) Pennipeds (combined total risk)") +
  annotation_north_arrow(location = "br", # North arrow
                         pad_y=unit(1, "cm"),
                         width = unit(0.6, "cm"),
                         height = unit(0.8, "cm")) +
  annotation_scale(location = "br",
                   height = unit(0.2, "cm")) +
  theme_bw()  +
  labs(x = "Longitude", y="Latitude")

# Plot all using the operators "|" (side-by-side) and "/" (one-over-another)

(sp1 | sp2) / tot









