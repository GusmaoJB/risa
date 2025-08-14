# Create test data
vec1 <- df_to_shp(data.frame(long = c(1,2,2,4), lat = c(4,4,2,2)))
vec2 <- df_to_shp(data.frame(long = c(2,5,4,6), lat = c(4,4,2,2)))
vec_list <- list(vec1, vec2)
vec_list

df <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
vec <- df_to_shp(df)             # silent
vec2 <- df_to_shp(df, quiet=FALSE) # verbose
vec2

df <- data.frame(
  long = c(1,2,2,4, 2,5,4,6),
  lat  = c(4,4,2,2, 4,4,2,2),
  species = rep(c("sp1","sp2"), each = 4)
)
vec_list <- df_to_list(df, group = "species")
vec_list

# Not a list
list_depth_base(42)
# Empty list
list_depth_base(list())
# Nested list
nested <- list(a = list(b = list(c = 1)))
list_depth_base(nested)

# Create test data
coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
# Gussing coordinates
guess_crs(coords)

# Create test data
coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
coords_vec <- df_to_shp(coords)

# Transform vector to a metric projection
transform_to_metric(coords_vec)

df <- data.frame(long = c(277000,389000,389000,611000),
                 lat  = c(442000,442000,221000,221000))
sf_obj <- sf::st_as_sf(df, coords = names(df), crs = 32631, remove = FALSE)
convert_to_decimal_degrees(sf_obj)

# Create test data
coords <- data.frame(long = c(1,2,2,4), lat = c(4,4,2,2))
coords_vec <- df_to_shp(coords)
coords_vec_metric <- transform_to_metric(coords_vec)$shape
# Create area of interest (default: convex hull)
aoi <- create_area(coords)
# Plot results
plot(sf::st_geometry(aoi), border = "blue", axes=TRUE)
plot(sf::st_geometry(coords_vec_metric), add = TRUE, col = "red")

# Creating example data
# Example
r <- terra::rast(nrows=5, ncols=5, vals=c(0.01,2,3,4,5))
reclass_matrix(r)

# Creating test data
df <- data.frame(long = rnorm(120, 0, 10), lat = rnorm(120, 0, 10))
# Generating reclassified Kernel densities estimates (3 classes)
kde <- get_class_kernel(df)
# Plot KDE map
plot(kde)

species  <- data.frame(long = rnorm(80, 0, 10),  lat = rnorm(80, 0, 10))
stressor <- data.frame(long = rnorm(100, 0, 10), lat = rnorm(100, 0, 10))
kde_spe <- get_class_kernel(species,  output_layer_type = "raster")
kde_str <- get_class_kernel(stressor, output_layer_type = "raster")
# Ensure same grid (nearest to preserve class labels)
kde_str <- terra::project(kde_str, kde_spe, method = "near")
overlap <- get_overlap_kernel(kde_spe, kde_str, method = "product", output_layer_type = "raster")
terra::plot(overlap)
overlap


# Creating test data
# Creating test data
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "sp1"),
                           data.frame(long = rnorm(60, 0, 10),
                           lat = rnorm(60, 0, 10), species = "sp2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 10),
                            lat = rnorm(100, 0, 10), stressor = "trawling"),
                                              data.frame(long = rnorm(50, 0, 10),
                                                  lat = rnorm(100, 0, 5), stressor = "gillnet"))
                                                      # Create kernel maps of species and stressor distributions and overlap maps
                                                      risa_maps <- risa_prep(spp_df, str_df)
# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)
# Species and Stressor distributions
dev.off()
par(mfrow = c(2,2))
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 1")
plot(risa_maps$species_kernel_maps$sp1$shp, add = TRUE)
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species 2")
plot(risa_maps$species_kernel_maps$sp2$shp, add = TRUE)
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Gillnet")
plot(risa_maps$stressor_kernel_maps$gillnet$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Trawling")
plot(risa_maps$stressor_kernel_maps$trawling$shp, add = TRUE, col=c("lightgreen", "green", "darkgreen"))
# Overlap maps
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Gillnet")
plot(risa_maps$overlap_maps$sp1$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Gillnet")
plot(risa_maps$overlap_maps$sp2$gillnet$shp, add = TRUE, col=c("yellow", "orange", "red"))
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species1 vs Trawling")
plot(risa_maps$overlap_maps$sp1$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))
plot(sf::st_geometry(risa_maps$area_of_interest), border = "blue", axes=TRUE, main="Species2 vs Trawling")
plot(risa_maps$overlap_maps$sp2$trawling$shp, add = TRUE, col=c("yellow", "orange", "red"))






