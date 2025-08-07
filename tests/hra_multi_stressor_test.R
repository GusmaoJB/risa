############### Tesst Field ###############

# Creating test data
set.seed(12)
spp_df <- rbind(data.frame(long = rnorm(80, 0, 10),
                           lat = rnorm(80, 0, 10), species = "species1"))
str_df <- rbind(data.frame(long = rnorm(100, 0, 5),
                           lat = rnorm(100, 0, 10), stressor = "stressor1"),
                data.frame(long = rnorm(50, 0, 10),
                           lat = rnorm(100, 0, 5), stressor = "stressor2"))

# Create kernel maps of species and stressor distributions and overlap maps
risa_maps <- risa_prep(spp_df, str_df)

# export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")

#Load example data
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)

#Reshape criteria table
crit_list <- criteria_reshape(df)

test <- habitat_risk_assessment(risa_maps, crit_list[[1]], equation = "euclidean")

terra::plot(risa_maps$species_kernel_maps$species1$raster)
terra::plot(risa_maps$stressor_kernel_maps$stressor1$raster)
terra::plot(risa_maps$stressor_kernel_maps$stressor2$raster)

terra::plot(risa_maps$overlap_maps$species1$stressor1$raster)
terra::plot(risa_maps$overlap_maps$species1$stressor2$raster)

terra::plot(sqrt(((test$stressor1$E_criteria-1)^2 + (test$stressor1$C_criteria-1)^2)))


terra::plot(test$stressor1$E_criteria)
terra::plot(test$stressor1$C_criteria)

terra::plot(test$stressor1$Risk_map_raw)
terra::plot(test$stressor2$Risk_map_raw)

terra::plot(test$stressor1$Risk_map)
terra::plot(test$stressor2$Risk_map)
