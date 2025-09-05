############### Test Field ###############

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
# export to compare with invest
# export_maps(risa_maps, "C:/Users/gusma/Documents/research/test_hra/maps")

#Load example data
#path <- "C:/Users/gusma/Documents/research/test_hra/1sp_2stressors.csv"
path <- system.file("extdata", "multi_species_criteria.csv", package = "risa")
df <- read.csv(path)
#Reshape criteria table
crit_list <- criteria_reshape(df)

crit_list
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

terra::plot(res$stressor1$Risk_map)
terra::plot(res$stressor2$Risk_map)
terra::plot(res$total_raw)
terra::plot(res$total)


# per-stressor buffers (named)
res4 <- hra2(rast_list, spp_dist, crit_list,
             equation = "euclidean",
             decay = "linear",
             r_max = 3, n_overlap = 2,
             buffer_m = c(stressor1 = 500000, stressor2 = 1000000))

terra::plot(res4$species1$stressor1$E_criteria)
terra::plot(res4$species1$stressor1$C_criteria)

terra::plot(res4$species1$stressor2$E_criteria)
terra::plot(res4$species1$stressor2$C_criteria)

terra::plot(res4$species1$stressor1$Risk_map)
terra::plot(res4$species1$stressor2$Risk_map)

terra::plot(res4$species1$stressor1$Risk_map_raw)
terra::plot(res4$species1$stressor2$Risk_map_raw)

terra::plot(res4$species1$total_raw)
terra::plot(res4$species2$total_raw)

terra::plot(res4$species1$total)
terra::plot(res4$species2$total)

terra::plot(res4$ecosys_risk_classified)

