############### Test Field ###############
setwd("/home/jojo/Documents/pontal_projects/risa")
library(risa)
library(sf)

# Loading real data
avist_data <- read.csv("/home/jojo/Documents/pontal_projects/megacost/avistagens_atualizadas.csv")
sg_data <- avist_data[avist_data$sp == "sg", c("Easting", "Northing", "tam_max", "id_grupo", "data", "area_code", "best_transect")]
sg_data$tam_max <- as.numeric(sg_data$tam_max)
sg_data_clean <- aggregate(. ~ data + area_code + best_transect, data=sg_data[,-4], FUN=mean)
sg_data_clean$tam_max <- round(sg_data_clean$tam_max)
sg_data_cleaner <- sg_data_clean[,c("Easting", "Northing", "tam_max")]

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 1),
                           lat = rnorm(80, 0, 1), species = "species1"),
                data.frame(long = rnorm(60, 0, 1),
                           lat = rnorm(60, 0, 1), species = "species2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 1.5),
                           lat = rnorm(100, 0, 3), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 3),
                           lat = rnorm(100, 0, 1.5), stressor = "stressor2"))

test_data1 <- spp_df[spp_df$species == "species1",]
test_data2 <- test_data1[1:30,]
test_data3 <- test_data1[1:10,]
test_data4 <- test_data1[1:5,]

# Function to plot KDE and points
plot_kernel_points <- function(data,
                               x_col = "long",
                               y_col = "lat",
                               input_crs = "EPSG:4326",
                               point_col = "red",
                               point_pch = 16,
                               point_cex = 0.8,
                               ...) {

  kde_raster <- get_class_kernel2(
    data,
    input_crs = input_crs,
    output_layer_type = "raster",
    ...
  )

  pts <- terra::vect(
    data,
    geom = c(x_col, y_col),
    crs = input_crs
  )

  pts_proj <- terra::project(
    pts,
    terra::crs(kde_raster)
  )

  terra::plot(kde_raster)

  terra::plot(
    pts_proj,
    add = TRUE,
    pch = point_pch,
    col = point_col,
    cex = point_cex
  )

  invisible(list(
    raster = kde_raster,
    points = pts_proj
  ))
}

head(sg_data_cleaner)
test <- get_class_kernel2(sg_data_cleaner, input_crs = "EPSG:32722", output_layer_type = "raster")
test
plot_kernel_points(data = sg_data_cleaner, x_col = "Easting",
                   y_col = "Northing", input_crs = "EPSG:32722")

head(sg_data_cleaner)
test <- get_class_kernel2(sg_data_cleaner, input_crs = "EPSG:32722")
plot(test)

library(ggplot2)
ggplot() +
  geom_point(data = spp_df[spp_df$species == "species1",],
             aes(x=long, y=lat)) +
  geom_point(data = spp_df[1:5,],
             aes(x=long, y=lat), col="red")


  plot(rbind.data.frame(spp_df[spp_df$species == "species1",1:2],
                        spp_df[1:5,1:2]))



#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
crit_list <- criteria_reshape(df)

# Risa prep
input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])

#export
#risa::export_maps(input_maps, out_dir = "/home/jojo/Documents/pontal_projects/risa_example_maps")

# test
sp_sf <- input_maps$species_kernel_maps$species1$shp
plot(sp_sf)
kde_raster <- terra::rasterize(vect(sp_sf), rast(sp_sf, resolution=1000))
terra::plot(kde_raster)

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)
# export to compare with invest
# export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")

# Selecting spatially explicit criteria ratings
# Note that the rasters in the stressors's list are named after the respective attribute in the criteria table.
rast_list <- list(
   species1 = list(
     stressor1 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor1$raster),
     stressor2 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species1$stressor2$raster)
   ),
   species2 = list(
     stressor1 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor1$raster),
     stressor2 = list(
       intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
       `likelihood of interaction`=risa_maps$overlap_maps$species2$stressor2$raster)
   )
 )

# Species' distributions rasters
spp_dist <- list(species1 = risa_maps$species_distributions$species1$raster,
                 species2 = risa_maps$species_distributions$species2$raster)



# single value for all stressors, linear decay within 10 km
res <- hra(rast_list[[1]], spp_dist[[1]], crit_list[[1]],
          equation = "multiplicative",
          r_max = 3, n_overlap = 2)

# per-stressor buffers (named)
res_2_2 <- hra(rast_list, spp_dist, crit_list,
             equation = "euclidean",
             decay = "linear",
             r_max = 3, n_overlap = 2,
             buffer_m = c(stressor1 = 500000, stressor2 = 1000000))

spp_df

output <- quick_byra(spp_df, str_df, df)
output$summary_stats


risaplot(output$kde_maps)
risaplot(output)

input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])

test2 <- quick_byra2(risa_maps, criteria = df)
risaplot(test2)

r_stack <- terra::rast(input_maps$species_distributions)

input_maps$species_distributions$species1

raster_list <- reshape_risa_maps(input_maps, crit_names)
species_distr <- input_maps$species_distributions





"area_of_interest" %in% names(byra_test)

terra::plot(byra_test$species1$stressor1$E_criteria)
terra::plot(byra_test$species1$stressor1$C_criteria)

byra_test$summary_stats
terra::plot(byra_test$ecosys_risk_raw)
terra::plot(byra_test$ecosys_risk_classified)



