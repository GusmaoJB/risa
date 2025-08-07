############### Test Field ###############
library(risa)

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "species1"),
                data.frame(long = rnorm(60, 0, 10),
                           lat = rnorm(60, 0, 10), species = "species2"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           lat = rnorm(100, 0, 10), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 10),
                           lat = rnorm(100, 0, 5), stressor = "stressor2"))

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)

export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")

#Load example data
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)

#Reshape criteria table
crit_list <- criteria_reshape(df)

crit_df <- crit_list[[1]]
crit_df

# Selecting spatially explicit criteria ratings
rast_list <- list(
  species1 = list(
    stressor1 = list(
      intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
      `likelihood of interaction` = risa_maps$overlap_maps$species1$stressor1$raster),
    stressor2 = list(
      intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
      `likelihood of interaction` = risa_maps$overlap_maps$species1$stressor2$raster)
  ),
  species2 = list(
    stressor1 = list(
      intensity = risa_maps$stressor_kernel_maps$stressor1$raster,
      `likelihood of interaction` = risa_maps$overlap_maps$species2$stressor1$raster),
    stressor2 = list(
      intensity = risa_maps$stressor_kernel_maps$stressor2$raster,
      `likelihood of interaction` = risa_maps$overlap_maps$species2$stressor2$raster)
  )
)


# Species distribution
sp_dist <- risa_maps$species_distributions$species1$raster


test <- habitat_risk_assessment(rast_list, sp_dist, crit_df, equation = "multiplicative")

terra::plot(risa_maps$species_kernel_maps$species1$raster)
terra::plot(risa_maps$stressor_kernel_maps$stressor1$raster)
terra::plot(risa_maps$stressor_kernel_maps$stressor2$raster)

terra::plot(risa_maps$overlap_maps$species1$stressor1$raster)
terra::plot(risa_maps$overlap_maps$species1$stressor2$raster)

terra::plot(sqrt(((test$stressor1$E_criteria-1)^2 + (test$stressor1$C_criteria-1)^2)))


terra::plot(terra::mask(
  test$stressor1$C_criteria))


terra::plot(test$stressor1$C_criteria)

terra::plot(test$stressor1$Risk_map_raw)
terra::plot(test$stressor2$Risk_map_raw)

terra::plot(test$stressor1$Risk_map)
terra::plot(test$stressor2$Risk_map)

get_stats(test)
test$total
