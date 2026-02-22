############### Test Field ###############
setwd("/home/jojo/Documents/pontal_projects/risa")
library(risa)

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 3),
                           lat = rnorm(80, 0, 3), species = "species1"),
                data.frame(long = rnorm(60, 0, 3),
                           lat = rnorm(60, 0, 3), species = "species2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 1.5),
                           lat = rnorm(100, 0, 3), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 3),
                           lat = rnorm(100, 0, 1.5), stressor = "stressor2"))


#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
crit_list <- criteria_reshape(df)

# Risa prep
input_maps <- risa_prep(spp_df, str_df)
input_maps_single <- risa_prep(spp_df[,-3], str_df[,-3])

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



